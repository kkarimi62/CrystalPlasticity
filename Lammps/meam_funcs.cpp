/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Sebastian Hütter (OvGU)
------------------------------------------------------------------------- */

#include "meam.h"
#include "math_special.h"
#include <cmath>
#include <cassert>
#include <math.h>
#include <stdio.h>

using namespace LAMMPS_NS;


//-----------------------------------------------------------------------------
// Compute G(gamma) based on selection flag ibar:
//   0 => G = sqrt(1+gamma)
//   1 => G = exp(gamma/2)
//   2 => not implemented
//   3 => G = 2/(1+exp(-gamma))
//   4 => G = sqrt(1+gamma)
//  -5 => G = +-sqrt(abs(1+gamma))
//
double
MEAM::G_gam(const double gamma, const int ibar, int& errorflag) const
{
  double gsmooth_switchpoint;

  switch (ibar) {
    case 0:
    case 4:
      gsmooth_switchpoint = -gsmooth_factor / (gsmooth_factor + 1);
      if (gamma < gsmooth_switchpoint) {
        //         e.g. gsmooth_factor is 99, {:
        //         gsmooth_switchpoint = -0.99
        //         G = 0.01*(-0.99/gamma)**99
        double G = 1 / (gsmooth_factor + 1) * pow((gsmooth_switchpoint / gamma), gsmooth_factor);
        return sqrt(G);
      } else {
        return sqrt(1.0 + gamma);
      }
    case 1:
      return MathSpecial::fm_exp(gamma / 2.0);
    case 3:
      return 2.0 / (1.0 + MathSpecial::fm_exp(-gamma));
    case -5:
      if ((1.0 + gamma) >= 0) {
        return sqrt(1.0 + gamma);
      } else {
        return -sqrt(-1.0 - gamma);
      }
  }
  errorflag = 1;
  return 0.0;
}

//-----------------------------------------------------------------------------
// Compute G(gamma and dG(gamma) based on selection flag ibar:
//   0 => G = sqrt(1+gamma)
//   1 => G = exp(gamma/2)
//   2 => not implemented
//   3 => G = 2/(1+exp(-gamma))
//   4 => G = sqrt(1+gamma)
//  -5 => G = +-sqrt(abs(1+gamma))
//
double
MEAM::dG_gam(const double gamma, const int ibar, double& dG, double& ddG) const
{
  double gsmooth_switchpoint;
  double G;

  switch (ibar) {
    case 0:
    case 4:
      gsmooth_switchpoint = -gsmooth_factor / (gsmooth_factor + 1);
      if (gamma < gsmooth_switchpoint) {
        //         e.g. gsmooth_factor is 99, {:
        //         gsmooth_switchpoint = -0.99
        //         G = 0.01*(-0.99/gamma)**99
        G = 1 / (gsmooth_factor + 1) * pow((gsmooth_switchpoint / gamma), gsmooth_factor);
        G = sqrt(G);
        dG = -gsmooth_factor * G / (2.0 * gamma);
        ddG = 0.5 * gsmooth_factor * ( G / gamma - dG ) / gamma;
        return G;
      } else {
        G = sqrt(1.0 + gamma);
        dG = 1.0 / (2.0 * G);
        ddG =  dG * (-1.0 / (2.0 * G * G ));
        return G;
      }
    case 1:
      G = MathSpecial::fm_exp(gamma / 2.0);
      dG = G / 2.0;
      ddG = dG / 2.0;
      return G;
    case 3:
      G = 2.0 / (1.0 + MathSpecial::fm_exp(-gamma));
      dG = G * (2.0 - G) / 2;
      ddG = dG * (1.0 - G);
      return G;
    case -5:
      if ((1.0 + gamma) >= 0) {
        G = sqrt(1.0 + gamma);
        dG = 1.0 / (2.0 * G);
        ddG =  dG * (-1.0 / (2.0 * G * G ));
        return G;
      } else {
        G = -sqrt(-1.0 - gamma);
        dG = -1.0 / (2.0 * G);
        ddG =  dG * (1.0 / (2.0 * G * G ));
        return G;
      }
  }
  dG = 1.0;
  ddG = 0.0; 
  return 0.0;
}

//-----------------------------------------------------------------------------
// Compute ZBL potential
//
double
MEAM::zbl(const double r, const int z1, const int z2)
{
  int i;
  const double c[] = { 0.028171, 0.28022, 0.50986, 0.18175 };
  const double d[] = { 0.20162, 0.40290, 0.94229, 3.1998 };
  const double azero = 0.4685;
  const double cc = 14.3997;
  double a, x;
  // azero = (9pi^2/128)^1/3 (0.529) Angstroms
  a = azero / (pow(z1, 0.23) + pow(z2, 0.23));
  double result = 0.0;
  x = r / a;
  for (i = 0; i <= 3; i++) {
    result = result + c[i] * MathSpecial::fm_exp(-d[i] * x);
  }
  if (r > 0.0)
    result = result * z1 * z2 / r * cc;
  return result;
}

//-----------------------------------------------------------------------------
// Compute embedding function F(rhobar) and derivative F'(rhobar), eqn I.5
//
double
MEAM::embedding(const double A, const double Ec, const double rhobar, double& dF, double& ddF ) const
{
  const double AEc = A * Ec;

  if (rhobar > 0.0) {
      const double lrb = log(rhobar);
      dF = AEc * (1.0 + lrb);
      ddF = AEc * (1.0/rhobar); //--- 2nd deriv.
      return AEc * rhobar * lrb;
  } else {
    if (this->emb_lin_neg == 0) {
      dF = 0.0;
      ddF = 0.0;
      return 0.0;
    } else {
      dF = - AEc;
      ddF = 0.0;
      return - AEc * rhobar;
    }
  }
}

//-----------------------------------------------------------------------------
// Compute Rose energy function, I.16
//
double
MEAM::erose(const double r, const double re, const double alpha, const double Ec, const double repuls,
            const double attrac, const int form)
{
  double astar, a3;
  double result = 0.0;

  if (r > 0.0) {
    astar = alpha * (r / re - 1.0);
    a3 = 0.0;
    if (astar >= 0)
      a3 = attrac;
    else if (astar < 0)
      a3 = repuls;

    if (form == 1)
      result = -Ec * (1 + astar + (-attrac + repuls / r) * MathSpecial::cube(astar)) * MathSpecial::fm_exp(-astar);
    else if (form == 2)
      result = -Ec * (1 + astar + a3 * MathSpecial::cube(astar)) * MathSpecial::fm_exp(-astar);
    else
      result = -Ec * (1 + astar + a3 * MathSpecial::cube(astar) / (r / re)) * MathSpecial::fm_exp(-astar);
  }
  return result;
}

//-----------------------------------------------------------------------------
// Shape factors for various configurations
//
void
MEAM::get_shpfcn(const lattice_t latt, const double sthe, const double cthe, double (&s)[3])
{
  switch (latt) {
    case FCC:
    case BCC:
    case B1:
    case B2:
      s[0] = 0.0;
      s[1] = 0.0;
      s[2] = 0.0;
      break;
    case HCP:
      s[0] = 0.0;
      s[1] = 0.0;
      s[2] = 1.0 / 3.0;
      break;
    case CH4: // CH4 actually needs shape factor for diamond for C, dimer for H
    case DIA:
    case DIA3:
      s[0] = 0.0;
      s[1] = 0.0;
      s[2] = 32.0 / 9.0;
      break;
    case DIM:
      s[0] = 1.0;
      s[1] = 2.0 / 3.0;
      //        s(3) = 1.d0 // this should be 0.4 unless (1-legendre) is multiplied in the density calc.
      s[2] = 0.40; // this is (1-legendre) where legendre = 0.6 in dynamo is accounted.
      break;
    case LIN: //linear, theta being 180
      s[0] = 0.0;
      s[1] = 8.0 / 3.0; // 4*(co**4 + si**4 - 1.0/3.0) in zig become 4*(1-1/3)
      s[2] = 0.0;
      break;
    case ZIG: //zig-zag
    case TRI: //trimer e.g. H2O
      s[0] = 4.0*pow(cthe,2);
      s[1] = 4.0*(pow(cthe,4) + pow(sthe,4) - 1.0/3.0);
      s[2] = 4.0*(pow(cthe,2) * (3*pow(sthe,4) + pow(cthe,4)));
      s[2] = s[2] - 0.6*s[0]; //legend in dyn, 0.6 is default value.
      break;
    default:
      s[0] = 0.0;
      //        call error('Lattice not defined in get_shpfcn.')
  }
}

//-----------------------------------------------------------------------------
// Number of first neighbors for reference structure
//
int
MEAM::get_Zij(const lattice_t latt)
{
  switch (latt) {
    case FCC:
      return 12;
    case BCC:
      return 8;
    case HCP:
      return 12;
    case DIA:
    case DIA3:
      return 4;
    case DIM:
      return 1;
    case B1:
      return 6;
    case C11:
      return 10;
    case L12:
      return 12;
    case B2:
      return 8;
    case CH4: // DYNAMO currently implemented this way while it needs two Z values, 4 and 1
      return 4;
    case LIN:
    case ZIG:
    case TRI:
      return 2;
      //        call error('Lattice not defined in get_Zij.')
  }
  return 0;
}

//-----------------------------------------------------------------------------
// Number of second neighbors for the reference structure
//   a = distance ratio R1/R2 (a2nn in dynamo)
//   numscr = number of atoms that screen the 2NN bond
//   S = second neighbor screening function (xfac, a part of b2nn in dynamo)
//
int
MEAM::get_Zij2(const lattice_t latt, const double cmin, const double cmax,
                const double stheta, double& a, double& S)
{

  double C, x, sijk;
  int Zij2 = 0, numscr = 0;

  switch (latt) {

  case FCC:
    Zij2 = 6;
    a = sqrt(2.0);
    numscr = 4;
    break;

  case BCC:
    Zij2 = 6;
    a = 2.0 / sqrt(3.0);
    numscr = 4;
    break;

  case HCP:
    Zij2 = 6;
    a = sqrt(2.0);
    numscr = 4;
    break;

  case B1:
    Zij2 = 12;
    a = sqrt(2.0);
    numscr = 2;
    break;

  case DIA: // 2NN
    Zij2 = 12;
    a = sqrt(8.0 / 3.0);
    numscr = 1;
    if (cmin < 0.500001) {
        //          call error('can not do 2NN MEAM for dia')
    }
    break;

  case DIA3: // 3NN
    Zij2 = 12;
    a = sqrt(11.0 / 3.0);
    numscr = 4;
    if (cmin < 0.500001) {
        //          call error('can not do 2NN MEAM for dia')
    }
    break;

  case CH4: //does not have 2nn structure so it returns 0
  case LIN: //line
  case DIM:
    //        this really shouldn't be allowed; make sure screening is zero
    a = 1.0;
    S = 0.0;
    return 0;

  case TRI: //TRI
    Zij2 = 1;
    a = 2.0*stheta;
    numscr=2;
    break;

  case ZIG: //zig-zag
    Zij2 = 2;
    a = 2.0*stheta;
    numscr=1;
    break;

  case L12:
    Zij2 = 6;
    a = sqrt(2.0);
    numscr = 4;
    break;

  case B2:
    Zij2 = 6;
    a = 2.0 / sqrt(3.0);
    numscr = 4;
    break;
  case C11:
    // unsupported lattice flag C11 in get_Zij
    break;
  default:
    // unknown lattic flag in get Zij
    //        call error('Lattice not defined in get_Zij.')
    break;
  }

  // Compute screening for each first neighbor
  if (latt==DIA3){ // special case for 3NN diamond structure
    C = 1.0;
  } else {
    C = 4.0 / (a * a) - 1.0;
  }
  x = (C - cmin) / (cmax - cmin);
  sijk = fcut(x);
  // There are numscr first neighbors screening the second neighbors
  S = MathSpecial::powint(sijk, numscr);
  return Zij2;
}

int
MEAM::get_Zij2_b2nn(const lattice_t latt, const double cmin, const double cmax, double &S)
{

  double x, sijk, C;
  int numscr = 0, Zij2 = 0;
  switch (latt) {
  case ZIG: //zig-zag for b11s and b22s term
  case TRI: //trimer for b11s
    Zij2 = 2;
    numscr=1;
    break;
  default:
    // unknown lattic flag in get Zij
    //        call error('Lattice not defined in get_Zij.')
    break;
  }
  C = 1.0;
  x = (C - cmin) / (cmax - cmin);
  sijk = fcut(x);
  S = MathSpecial::powint(sijk, numscr);
  return Zij2;
}

//------------------------------------------------------------
//-------- compute higher order derivatives of density
//------------------------------------------------------------
double
MEAM::Get_ddrho1drdr(int i, //--- deriv of 4.30(a) wrt rij
                    double rij, double sij, 
                    double rhoa1j, double drhoa1j, double ddrhoa1j,
                    double arg1i1,
                    double arg1i1_d
                    ){
        double rij2 = rij * rij;
        double a1 = 2 * sij / rij;
        double da1 = -2.0 * sij / rij2;
        double drho1dr1 = a1 * (drhoa1j - rhoa1j / rij) * arg1i1; //--- 4.30(a)
   
         double ddrho1ddr1 = da1 * (drhoa1j - rhoa1j / rij) * arg1i1 + 
                              a1 * (ddrhoa1j - (drhoa1j/rij-rhoa1j/rij2)) * arg1i1 +
                              a1 * (drhoa1j - rhoa1j / rij) * arg1i1_d;
        
//        double ddrho1ddr1 = a1 * ( ( - 0.5 * drho1dr1 / sij ) + ( ddrhoa1j - drhoa1j / rij + rhoa1j / rij2 ) * arg1i1 + (drhoa1j - rhoa1j / rij) * arg1i1_d ); 
        return ddrho1ddr1;
}

void
MEAM::Get_ddrho1drmdr(int i,
                       double rij, double sij, double* delij,
                       double rhoa1j, double drhoa1j,
                       double** arho1, double* darho1dri,
                       double* ddrho1drmdr1 //--- modify 
                    ){
        double a1;
        int m;
        a1 = 2.0 * sij / rij;
        double da1 = -a1/rij;
        for (m = 0; m < 3; m++) {
          ddrho1drmdr1[m] =  da1 * rhoa1j * arho1[i][m] + a1 * drhoa1j * arho1[i][m] + a1 * rhoa1j * darho1dri[m];
//          ddrho1drmdr1[m] = a1 * ( ( - rhoa1j * arho1[i][m] / rij ) + ( rhoa1j * darho1dri[m] ) + ( arho1[i][m] * drhoa1j ) );
        }
}
//---------------------------------------------------------------------------
void
MEAM::Get_ddrho1drmdrn(int i, //--- deriv of Eq. (4.30c) wrt rn
                       double rij, double sij,
                       double rhoa1j,
                       double* ddrho1drmdrn1 //--- modify 
                    ){
        double a1;
        int m, n, nv2;
        a1 = 2.0 * sij / rij;
        nv2 = 0;
        for (m = 0; m < 3; m++) {
         for (n = m; n < 3; n++) {
          ddrho1drmdrn1[ nv2 ] = a1 * ( rhoa1j * ( ( m == n ) ? 1.0 : 0.0 ) * sij / rij ) * rhoa1j;
          nv2++;
         }
        }
}


double
MEAM::Get_ddrho2drdr(int i, 
                    double rij, double sij, 
                    double rhoa2j, double drhoa2j, double ddrhoa2j,
                    double* arho2b,
                    double arg1i2,
                    double arg1i2_d
                    ){
        double rij2 = rij * rij;
        double a2 = 2 * sij / rij2;
        double darho2b = drhoa2j * sij;
        double argg_1 = (drhoa2j - 2 * rhoa2j / rij);
        double argg_2 = arg1i2;
        double argg_3 = 2.0 / 3.0 * arho2b[i] * drhoa2j * sij;
        double d_argg_1 = ddrhoa2j-2*(-rhoa2j/rij+drhoa2j)/rij;
        double d_argg_2 = arg1i2_d;
        double d_argg_3 = 2.0 / 3.0 * (darho2b * drhoa2j+arho2b[i]*ddrhoa2j) * sij;
        double ddrho2ddr1 = (-2 * a2 / rij) * argg_1 * argg_2 + a2 * (d_argg_1*argg_2+argg_1*d_argg_2)-d_argg_3;
        
        return ddrho2ddr1;
}

void
MEAM::Get_ddrho2drmdr(int i,
                       double rij, double sij, double* delij,
                       double rhoa2j, double drhoa2j,
                       double** arho2,
                       double* darho2dr,
                       double* ddrho2drmdr1 //--- modify 
                    ){
         double rij2 = rij * rij;
         double a2 = 4 * sij / rij2;
         double drho2drm1[3];
         int m, n;
        for (m = 0; m < 3; m++) {
          drho2drm1[m] = 0.0;
          ddrho2drmdr1[m] = 0.0;
          for (n = 0; n < 3; n++) {
            drho2drm1[m] += arho2[i][this->vind2D[m][n]] * delij[n]; //--- 4.30(f): arho2 is Y_{2i\sigma\alpha}
            ddrho2drmdr1[m] += darho2dr[this->vind2D[m][n]] * delij[n]; //this->vind2D[m][n]??????????
          }
          //--- d^2rho/drm/dr
          ddrho2drmdr1[m] *= rhoa2j;
          ddrho2drmdr1[m] += drho2drm1[m] * (drhoa2j-2*rhoa2j/rij); 
          ddrho2drmdr1[m] *= a2;
        }
}        
//-------------------
void
MEAM::Get_ddrho2drmdrn(int i,
                       double rij, double sij, double* delij,
                       double rhoa2j,
                       double** arho2,
                       double* ddrho2drmdrn1 //--- modify 
                    ){
        double drho2drm1[3];
        int m, n, nv2;
        double arg;
        double rij2 = rij * rij;
        double a2 = 4 * sij / rij2;
        nv2 = 0.0;
        for (m = 0; m < 3; m++) {
//          drho2drm1[m] = 0.0;
//          ddrho2drmdr1[m] = 0.0;
          for (n = 0; n < 3; n++) {
//             arg = rhoa2j * sij * ( delij[ m ] * delij[ n ] / rij2 + (m == n ? 1 : 0) );
//             ddrho2drmdrn1[ nv2 ] = a2 * rhoa2j * ( arg + arho2[i][this->vind2D[n][m]] ); //???
             ddrho2drmdrn1[ nv2 ] = (4/rij2)*rhoa2j*rhoa2j*sij*sij*(m == n ? 1 : 0);
            nv2++;
          }
        }
}     
//-------------------
double
MEAM::Get_ddrho3drdr( double rij, double sij, 
                    double rhoa3j, double drhoa3j, double ddrhoa3j,
                    double arg1i3, double arg3i3,
                    double arg1i3_d, double arg3i3_d
                    ){
        double rij2 = rij * rij;
        double rij3 = rij * rij2;
        double a3 = 2 * sij / rij3;
        double a3a = 6.0 / 5.0 * sij / rij;

        //--- 2nd deriv. wrt r
        double argg_1 = drhoa3j - 3 * rhoa3j / rij;
        double argg_2 = arg1i3;
        double argg_3 = drhoa3j - rhoa3j / rij;
        double argg_4 = arg3i3;
        double argg_1_d = ddrhoa3j - 3.0 * (- rhoa3j / rij + drhoa3j ) / rij;
        double argg_2_d = arg1i3_d;
        double argg_3_d = ddrhoa3j - ( drhoa3j - rhoa3j / rij ) / rij;
        double argg_4_d = arg3i3_d;
        double ddrho3ddr1 = a3  * (-3.0*argg_1*argg_2/rij+argg_1_d*argg_2+argg_1*argg_2_d)-
                     a3a * (-argg_3*argg_4/rij+argg_3_d*argg_4+argg_3*argg_4_d); //--- deriv. 4.30(g) wrt r
        return ddrho3ddr1;
}
//-------------------------------
void
MEAM::Get_ddrho3drmdr( int i,
                       double rij, double sij, double* delij,
                       double rhoa3j, 
                       double drhoa3j, 
                       double* darho3dr, double* darho3bdr,
                       double* ddrho3drmdr1 //--- modify 
                     ){
        double drho3drm1[3];
        int m, n, p, nv2;
        double arg;
        double rij2 = rij * rij;
        double rij3 = rij * rij2;
        double a3 = 6 * sij / rij3;
        double a3a = 6 * sij / (5 * rij);
        for (m = 0; m < 3; m++) {
          drho3drm1[m] = 0.0;
          ddrho3drmdr1[m] = 0.0;
          nv2 = 0;
          for (n = 0; n < 3; n++) {
            for (p = n; p < 3; p++) {
              arg = delij[n] * delij[p] * this->v2D[nv2];
              drho3drm1[m] = drho3drm1[m] + arho3[i][this->vind3D[m][n][p]] * arg; //--- 4.30(i)
              //
              ddrho3drmdr1[m] += darho3dr[this->vind3D[m][n][p]] * arg; //--- deriv. 4.30(i) wrt. r
              //
              nv2 = nv2 + 1;
            }
          }
          //
          ddrho3drmdr1[m] *= rhoa3j;  //--- deriv. 4.30(i) wrt. r
          ddrho3drmdr1[m] += drho3drm1[m] * (-3*rhoa3j/rij+drhoa3j);
          ddrho3drmdr1[m] *= a3;
          ddrho3drmdr1[m] += -a3a * ((-rhoa3j/rij+drhoa3j)*arho3b[i][m]+rhoa3j*darho3bdr[m]);
//           //
//           ddrho3drmdr2[m] *= rhoa3i;  //--- deriv. 4.30(i) wrt. r (atom j)
//           ddrho3drmdr2[m] += drho3drm2[m] * (-3*rhoa3i/rij+drhoa3i);
//           ddrho3drmdr2[m] *= a3;
//           ddrho3drmdr2[m] += a3a * ((-rhoa3i/rij+drhoa3i)*arho3b[j][m]+rhoa3i*darho3bdr[j][m]) //--- negative sign???
          //
        }
}
//-------------------------------
void
MEAM::Get_ddrho3drmdrn( int i,
                       double rij, double sij, double* delij,
                       double rhoa3j, 
                       double* ddrho3drmdrn1 //--- modify 
                     ){
        double drho3drm1[3];
        int m, n, p, nv2;
        double arg1;
        double rij2 = rij * rij;
        double rij3 = rij * rij2;
        double a3 = 6 * sij / rij3;
        double a3a = 6 * sij / (5 * rij);
   
//         nv3 = 0;
//         for (m = 0; m < 3; m++) {
//           for (n = m; n < 3; n++) {
//             arg1 = 0.0;
//             for (p = n; p < 3; p++) {
// //              arg = delij[n] * delij[p] * this->v2D[nv2]; //--- ?????????
//               arg1 +=  ( arho3[i][this->vind3D[m][n][p]] + arho3[i][this->vind3D[m][p][n]] ) * delij[p];
//             }
//             arg1 +=  rhoa3j * sij * ((m == n?1 : 0) * rij2 + 2.0 * delij[m]*delij[n]) / rij;
//             ddrho3drmdrn1[nv3] = rhoa3j * (a3 * arg1 -  a3a * rhoa3j * sij * ( m == n?1 : 0 ) / rij );
//             nv3++;
//           }
         for (m = 0; m < 3; m++) {
          for (n = m; n < 3; n++) {
            ddrho3drmdrn1[nv2] = rhoa3j * (a3*(rhoa3j*sij*rij*( m == n ? 1 : 0 )+4*(rhoa3j*sij/rij)*(delij[m]*delij[n]))
                                           -a3a*(rhoa3j*sij*( m == n ? 1 : 0 )/rij));
            nv2++;
          }
         }
           //
 //--- negative sign???
          //
//        }
}
//-------------------------------

double
MEAM::Get_ddrhodrdr(  int i, int elti,
                      double* shpi,
                      double t1i, double t2i, double t3i,
                      double dt1dr1, double dt2dr1, double dt3dr1,
                      double ddt1drdr1, double ddt2drdr1, double ddt3drdr1,
                      double* rho0, double* rho1, double* rho2, double* rho3, 
                      double drho0dr1, double drho1dr1, double drho2dr1, double drho3dr1, 
                      double ddrho0drdr1, double ddrho1drdr1, double ddrho2drdr1, double ddrho3drdr1,
                      double drhodr1
                     ){

        double dgamdr = (shpi[0] * dt1dr1 + shpi[1] * dt2dr1 + shpi[2] * dt3dr1) / (Zarray[i] * Zarray[i]); //--- d\Gamma^{ref}/dr: deriv of Eq. (4.6)
        double ddgamdrdr = (shpi[0] * ddt1drdr1 + shpi[1] * ddt2drdr1 + shpi[2] * ddt3drdr1) / (Zarray[i] * Zarray[i]); //--- d^2\Gamma^{ref}/drdr: 2nd deriv of Eq. (4.6)
        double rho_bkgd = rho_bkgd_array[ i ];
        double drho_bkgd_dr = this->rho0_meam[elti] * Zarray[ i ] * dGbar_array[ i ] * dgamdr; //--- deriv of Eq. (4.5) wrt. r
        double ddrho_bkgd_drdr = this->rho0_meam[elti] * Zarray[ i ] * ( ddGbar_array[i] * dgamdr * dgamdr + dGbar_array[ i ] * ddgamdrdr ); //--- 2nd deriv of Eq. (4.5) wrt. r
        double dgammadr =  (dt1dr1 * rho1[i] + t1i * drho1dr1 + //--- deriv. of Eq. (4.4) wrt. r
                            dt2dr1 * rho2[i] + t2i * drho2dr1 + 
                            dt3dr1 * rho3[i] + t3i * drho3dr1) - 2.0 * rho0[i] * drho0dr1 *  gamma[i]; 
        if (rho0[i] > 0.0) {
        dgammadr /= (rho0[i] * rho0[i]);
        }
        //--- 2nd deriv. of Eq. (4.4) wrt. r (index i)
        double argg1 = ddt1drdr1 * rho1[i] + 2.0 * dt1dr1 * drho1dr1 + t1i * ddrho1drdr1 +
                ddt2drdr1 * rho2[i] + 2.0 * dt2dr1 * drho2dr1 + t2i * ddrho2drdr1 +
                ddt3drdr1 * rho3[i] + 2.0 * dt3dr1 * drho3dr1 + t3i * ddrho3drdr1;
        double argg2 = drho0dr1 * drho0dr1 *  gamma[ i ] +
                rho0[i]  * ddrho0drdr1 * gamma[ i ] +
                2.0 * rho0[i] * drho0dr1 * dgammadr;
        double ddgammadrdr = argg1 - 2.0 * argg2;
        if (rho0[i] > 0.0) {
        ddgammadrdr /= (rho0[i] * rho0[i]);
        }        
        //--- total deriv
        double ddrhodrdr1    =  ddrho0drdr1 * G_array[i] + 
                         2.0 * drho0dr1 * dG_array[i] * dgammadr + 
                         rho0[i] * ddG_array[i] * dgammadr * dgammadr + 
                         rho0[i] * dG_array[i] * ddgammadrdr -
                         2.0 * drhodr1 * drho_bkgd_dr -
                         rho[ i ] * ddrho_bkgd_drdr;        
        ddrhodrdr1 /= rho_bkgd; //rho_bkgd defined?
        return ddrhodrdr1;
        


//       ddgamma1[i] = (dGdr - 2 * ddGdrdr * gamma[i]- 2 * dGdr * dgamma[i]) * denom +
//                    (G - 2 * dGdr * gamma[i]) * ddenom; //--- Eq. (4.36a): prefactor in the 1st term of the RHS 

        
        
//         ddrhodrdr1 = ddgamma1[i] * drho0dr1 + dgamma1[i] * ddrho0dr1 +
//                      ddgamma2[i] * (dt1dr1 * rho1[i] + t1i * drho1dr1 + dt2dr1 * rho2[i] + t2i * drho2dr1 + dt3dr1 * rho3[i] + t3i * drho3dr1) + 
//                      dgamma2[i] * (ddt1drdr1 * rho1[i] + dt1dr1 * drho1dr[i]+ 
//                                    dt1idr * drho1dr1 + t1i * ddrho1drdr1+
//                                    ddt2drdr1 * rho2[i] + dt2dr1 * drho2dr[i] + 
//                                    dt2idr * drho2dr1 + t2i * ddrho2drdr1 +
//                                    ddt3drdr1 * rho3[i] + dt3dr1 * drho3dr[i]+ 
//                                    dt3idr * drho3dr1 + t3i * ddrho3drdr1) -
//                      ddgamma3[i] * (shpi[0] * dt1dr1 + shpi[1] * dt2dr1 + shpi[2] * dt3dr1) -
//                      dgamma3[i] * (shpi[0] * ddt1drdr1 + shpi[1] * ddt2drdr1 + shpi[2] * ddt3drdr1); //shpi(r)? 
}





//-----------------------------------------------------------------------------
void
MEAM::Get_ddrhodrmdr( int i, int elti, //--- deriv. of Eq. 4.36(c) wrt. r
                        double* shpi, 
                        double t1i,  double t2i,  double t3i,
                        double dt1dr,  double dt2dr,  double dt3dr,
                        double* rho0, double* rho1, double* rho2, double* rho3,
                        double drho0dr,  double drho1dr,  double drho2dr,  double drho3dr, 
                        double* drho0drm,  double* drho1drm,  double* drho2drm,  double* drho3drm, 
                        double* ddrho0drmdr, double* ddrho1drmdr,  double* ddrho2drmdr,  double* ddrho3drmdr,
                        double* drhodrm,
                        double* ddrhodrmdr //--- modify
                       ){
          double LHS, dLHS;
          double dgamdr = (shpi[0] * dt1dr + shpi[1] * dt2dr + shpi[2] * dt3dr) / (Zarray[i] * Zarray[i]); //--- d\Gamma^{ref}/dr: deriv of Eq. (4.6)
          //
          double dgammadr =  (dt1dr * rho1[i] + t1i * drho1dr + //--- deriv. of Eq. (4.4) wrt. r
                              dt2dr * rho2[i] + t2i * drho2dr + 
                              dt3dr * rho3[i] + t3i * drho3dr) - 2.0 * rho0[i] * drho0dr *  gamma[i]; 
          if (rho0[i] > 0.0) {
            dgammadr /= (rho0[i] * rho0[i]);
          }
          //
          double drho_bkgd_dr = this->rho0_meam[elti] * Zarray[ i ] * dGbar_array[ i ] * dgamdr; //--- deriv of Eq. (4.5) wrt. r
          //
          for( int m = 0; m < 3; m++ ){
            LHS = drhodrm[m];
            dLHS = - LHS * ( drho_bkgd_dr * rho0[i] + rho_bkgd_array[i] * drho0dr ) +
                     ddG_array[i] * dgammadr * (t1i*drho1drm[m] + t2i*drho2drm[m] + t3i*drho3drm[m])+
                     dG_array[i] * ( dt1dr * drho1drm[m] + dt2dr * drho2drm[m] + dt3dr * drho3drm[m] +
                                     t1i * ddrho1drmdr[m] + t2i * ddrho1drmdr[m] + t3i * ddrho1drmdr[m] );
            dLHS /= (rho_bkgd_array[i] * rho0[i]);
            ddrhodrmdr[ m ] = dLHS;
          }
}
//-----------------------------------------------------------------------------
void
MEAM::Get_ddrhodrmdrn( int i, int elti, //--- deriv. of Eq. 4.36(c) wrt. rm
                        double* shpi, 
                        double t1i,  double t2i,  double t3i,
                        double* rho0, double* rho1, double* rho2, double* rho3,
                        double* drho0drm1,  double* drho1drm1,  double* drho2drm1,  double* drho3drm1, 
                                            double* ddrho1drmdrn1,  double* ddrho2drmdrn1,  double* ddrho3drmdrn1,
                        double* ddrhodrmdrn1 //--- modify
                       ){
          int ii = 0;
          double arg, dgammadrm[ 3 ] ; 
          //--- compute dgammadrm
          for( int m = 0; m < 3; m++ ){
            dgammadrm[ m ] =  (t1i * drho1drm1[m] + //--- deriv. of Eq. (4.4) wrt. rm
                               t2i * drho2drm1[m] + 
                               t3i * drho3drm1[m]); 
            dgammadrm[ m ] /= (rho0[i] * rho0[i]);
          }
         //--- compute ddrhodrmdrn1: deriv. of Eq. 4.36(c) wrt. rm
          for( int m = 0; m < 3; m++ ){
             arg = (t1i * drho1drm1[m] + t2i * drho2drm1[m] + t3i * drho3drm1[m] );
            for( int n = m; n < 3; n++ ){
               ddrhodrmdrn1[ii] = ddG_array[i] * dgammadrm[n] * arg +
                                  dG_array[i] * ( t1i * ddrho1drmdrn1[ii] + t2i * ddrho1drmdrn1[ii] + t3i * ddrho1drmdrn1[ii]);
               ddrhodrmdrn1[ii] /= (rho_bkgd_array[i] * rho0[i]);
               ii++;
            }
          }
}

//-----------------------------------------------------------------------------
double
MEAM::Get_ddrhodrds(  int i, int elti,
                      double* shpi,
                      double t1i, double t2i, double t3i,
                      double dt1dr1, double dt2dr1, double dt3dr1,
                      double dt1ds1, double dt2ds1, double dt3ds1,
                      double ddt1drds1, double ddt2drds1, double ddt3drds1, 
                      double* rho0, double* rho1, double* rho2, double* rho3, 
                      double drho0dr1, double drho1dr1, double drho2dr1, double drho3dr1, 
                      double drho0ds1, double drho1ds1, double drho2ds1, double drho3ds1, 
                      double ddrho0drds1, double ddrho1drds1, double ddrho2drds1, double ddrho3drds1,
                      double drhodr1, double drhods1
                     ){

        double dgamdr = (shpi[0] * dt1dr1 + shpi[1] * dt2dr1 + shpi[2] * dt3dr1) / (Zarray[i] * Zarray[i]); //--- d\Gamma^{ref}/dr: deriv of Eq. (4.6)
        double dgamds = (shpi[0] * dt1ds1 + shpi[1] * dt2ds1 + shpi[2] * dt3ds1) / (Zarray[i] * Zarray[i]); //--- d\Gamma^{ref}/ds: deriv of Eq. (4.6)
        double ddgamdrds = (shpi[0] * ddt1drds1 + shpi[1] * ddt2drds1 + shpi[2] * ddt3drds1) / (Zarray[i] * Zarray[i]); //--- d^2\Gamma^{ref}/drds: 2nd deriv of Eq. (4.6)
        double rho_bkgd = rho_bkgd_array[ i ];
        double drho_bkgd_dr = this->rho0_meam[elti] * Zarray[ i ] * dGbar_array[ i ] * dgamdr; //--- deriv of Eq. (4.5) wrt. r
        double drho_bkgd_ds = this->rho0_meam[elti] * Zarray[ i ] * dGbar_array[ i ] * dgamds; //--- deriv of Eq. (4.5) wrt. s
        double ddrho_bkgd_drds = this->rho0_meam[elti] * Zarray[ i ] * ( ddGbar_array[i] * dgamds * dgamdr + dGbar_array[ i ] * ddgamdrds ); //--- 2nd deriv of Eq. (4.5) wrt. r
        double dgammadr =  (dt1dr1 * rho1[i] + t1i * drho1dr1 + //--- deriv. of Eq. (4.4) wrt. r
                            dt2dr1 * rho2[i] + t2i * drho2dr1 + 
                            dt3dr1 * rho3[i] + t3i * drho3dr1) - 2.0 * rho0[i] * drho0dr1 *  gamma[i]; 
        double dgammads =  (dt1ds1 * rho1[i] + t1i * drho1ds1 + //--- deriv. of Eq. (4.4) wrt. s
                            dt2ds1 * rho2[i] + t2i * drho2ds1 + 
                            dt3ds1 * rho3[i] + t3i * drho3ds1) - 2.0 * rho0[i] * drho0ds1 *  gamma[i];                             
                            

        //--- d^2\Gamma/drds        
          double ddgammadrds   =  
                            (ddt1drds1 * rho1[i] + dt1dr1 * drho1ds1 + 
                             dt1ds1 * drho1dr1 + t1i * ddrho1drds1   +
                             ddt2drds1 * rho2[i] + dt2dr1 * drho2ds1 + 
                             dt2ds1 * drho2dr1 + t2i * ddrho2drds1   +
                             ddt3drds1 * rho3[i] + dt3dr1 * drho3ds1 +
                             dt3ds1 * drho3dr1 + t3i * ddrho3drds1)  -
                             2.0 * (drho0ds1 * drho0dr1 *  gamma[i] + rho0[i] * ddrho0drds1 *  gamma[i] + rho0[i] * drho0dr1 *  dgammads) -
                             dgammadr * (2.0 * rho0[i] * drho0ds1);
        if (!iszero(rho0[i])) {
           double rho0sq = (rho0[i] * rho0[i]);
           dgammadr /= rho0sq;
           dgammads /= rho0sq;
           ddgammadrds /= rho0sq;
        }
        else {
           dgammadr = dgammads = ddgammadrds = 0.0;
        }
        
        //--- total deriv
        double ddrhodrds1    =  ddrho0drds1 * G_array[i] + 
                         drho0dr1 * dG_array[i] * dgammads + 
                         drho0ds1 * dG_array[i] * dgammadr + 
                         rho0[i] * ddG_array[i] * dgammads * dgammadr + 
                         rho0[i] * dG_array[i] * ddgammadrds -
                         drhodr1 * drho_bkgd_ds -
                         drhods1 * drho_bkgd_dr -
                         rho[ i ] * ddrho_bkgd_drds;        
        ddrhodrds1 /= rho_bkgd;
        return ddrhodrds1;
}
//-----------------------------------------------------------------------------
void
MEAM::Get_ddrhodrmds( int i, int elti, //--- deriv. of Eq. 4.36(c) wrt. r
                        double* shpi, 
                        double t1i,  double t2i,  double t3i,
                        double dt1ds,  double dt2ds,  double dt3ds,
                        double* rho0, double* rho1, double* rho2, double* rho3,
                        double drho0ds,  double drho1ds,  double drho2ds,  double drho3ds, 
                        double* drho0drm,  double* drho1drm,  double* drho2drm,  double* drho3drm, 
                        double* ddrho0drmds, double* ddrho1drmds,  double* ddrho2drmds,  double* ddrho3drmds,
                        double* drhodrm,
                        double* ddrhodrmds //--- modify
                       ){
          double LHS, dLHS;
          double dgamds = (shpi[0] * dt1ds + shpi[1] * dt2ds + shpi[2] * dt3ds) / (Zarray[i] * Zarray[i]); //--- d\Gamma^{ref}/dr: deriv of Eq. (4.6)
          //
          double dgammads =  (dt1ds * rho1[i] + t1i * drho1ds + //--- deriv. of Eq. (4.4) wrt. r
                              dt2ds * rho2[i] + t2i * drho2ds + 
                              dt3ds * rho3[i] + t3i * drho3ds) - 2.0 * rho0[i] * drho0ds *  gamma[i]; 
          if (rho0[i] > 0.0) {
            dgammads /= (rho0[i] * rho0[i]);
          }
          else {
            dgammads = 0.0;
          }
          //
          double drho_bkgd_ds = this->rho0_meam[elti] * Zarray[ i ] * dGbar_array[ i ] * dgamds; //--- deriv of Eq. (4.5) wrt. r
          //
          for( int m = 0; m < 3; m++ ){
            LHS = drhodrm[m];
            dLHS = - LHS * ( drho_bkgd_ds * rho0[i] + rho_bkgd_array[i] * drho0ds ) +
                     ddG_array[i] * dgammads * (t1i*drho1drm[m] + t2i*drho2drm[m] + t3i*drho3drm[m])+
                     dG_array[i] * ( dt1ds * drho1drm[m] + dt2ds * drho2drm[m] + dt3ds * drho3drm[m] +
                                     t1i * ddrho1drmds[m] + t2i * ddrho1drmds[m] + t3i * ddrho1drmds[m] );
            dLHS /= (rho_bkgd_array[i] * rho0[i]);
            ddrhodrmds[ m ] = dLHS;
          }
}

double MEAM::GetModulus(int i, int j, double** x, int numneigh, int* firstneigh, int numneigh_full, int* firstneigh_full, int* type, int* fmap, double sij,
                        int alpha, int beta, int gamma, int lambda,  double r3,double ds, double dds, double recip,
                        double dUdrij, double dUdsij, double ddUddrij, double ddUdrijds, double ddUddsij,
                        double* dUdrijm, double* delij, double* ddUdrdrijm, double* ddUdrijmds, double* ddUdrmdrn){
     int nv2=0,m,n;
     for(m=0;m<alpha+1;m++){
       for(n=m;n<3;n++){
         if(m==alpha and n==gamma)
          break;
         nv2++;
       }
     }

   //
   double ddUdrijmds_tmp[ 3 ];
   if (iszero(ds)) {
     dUdsij = 0.0;
     ddUddsij = 0.0;
     for (m = 0; m < 3; m++){
        ddUdrijmds_tmp[m] = ddUdrijmds[m];
        ddUdrijmds[m] = 0.0;
     }
     ddUdrijds = 0.0;
    }
   
   double dsg_alpha_beta_drm[3];
   double recip2 = recip * recip;
   double dsg_alpha_beta_dr = ((-recip2*(dUdrij+dUdsij*ds)+recip*(ddUddrij+ddUdrijds*ds+dUdsij*dds))*delij[alpha]+ddUdrdrijm[alpha])*delij[beta];
   double dsg_alpha_beta_ds = (recip*(ddUdrijds+ddUddsij*ds)*delij[alpha]+ddUdrijmds[alpha])*delij[beta];
   dsg_alpha_beta_drm[gamma] = (recip*((ddUdrdrijm[gamma]+ddUdrijmds[gamma]*ds)*delij[alpha]+(dUdrij+dUdsij*ds)*(alpha == gamma ? 1 : 0))+ddUdrmdrn[nv2])*delij[beta];

   double mod2bdy = (recip*(dsg_alpha_beta_dr+dsg_alpha_beta_ds*ds)*delij[gamma]+dsg_alpha_beta_drm[gamma])*delij[lambda];
   double mod3bdy = 0.0;
  
  for (m = 0; m < 3; m++){
     ddUdrijmds[m] = ddUdrijmds_tmp[m];
  }
   
//    double arg1 = recip;
//      double arg2 = dUdrij + dUdsij * ds; // units of ds
//      double arg3 = dUdrijm[ alpha ];
//      double darg2 =(recip * ( ddUddrij + ddUdrijds * ds + ds * ( ddUdrijds + ddUddsij * ds ) + 1.0 * dUdsij * dds ) * delij[gamma]+
//             ( ddUdrdrijm[gamma]+ddUdrijmds[gamma]*ds))*delij[lambda];
//      double darg3 = (recip*(ddUdrdrijm[alpha]+ddUdrijmds[alpha]*ds)*delij[gamma]+ddUdrmdrn[nv2])*delij[lambda];  
//      return
//        (((-delij[gamma]*delij[lambda]/r3)*arg2+recip*darg2)*delij[alpha] + 
//         recip*arg2*(alpha == gamma ? 1 : 0)*delij[lambda]+ darg3)*delij[beta];

        //     Now compute forces on other atoms k due to change in sij     stiffness ?
        
        if (iszero(sij) || isone(sij) ) return mod2bdy; //: cont jn loop
        double dxik(0), dyik(0), dzik(0), dsij1, deljk[3], delki[3];
        double dxjk(0), dyjk(0), dzjk(0), dsij2;
        double rij2 = 1.0/recip2, rij=1.0/recip;
        double da,dcikj,dsikj,dsij,ddsij1drij,ddsij2drij,dsg_alpha_beta_drjk,dsg_alpha_beta_drik;
        int kn,k,eltk,elti,eltj;
        elti = fmap[type[i]];
        eltj = fmap[type[j]];
         
        for (kn = 0; kn < numneigh_full; kn++) {
          k = firstneigh_full[kn];
          eltk = fmap[type[k]];
          if (k != j && eltk >= 0) {
            double xik, xjk, cikj, sikj, dfc, ddfc, a;
            double dCikj1, dCikj2;
            double ddCikj1, ddCikj2;
            double delc, rik2, rjk2, rik, rjk;
            
            const double Cmax = this->Cmax_meam[elti][eltj][eltk];
            const double Cmin = this->Cmin_meam[elti][eltj][eltk];

            dsij1 = 0.0;
            dsij2 = 0.0;
            dsg_alpha_beta_drik=0.0;
            dsg_alpha_beta_drjk=0.0;
            if (!iszero(sij) && !isone(sij)) {
              const double rbound = rij2 * this->ebound_meam[elti][eltj];
              delc = Cmax - Cmin;
              dxjk = deljk[0]=x[k][0] - x[j][0];
              dyjk = deljk[1]=x[k][1] - x[j][1];
              dzjk = deljk[2]=x[k][2] - x[j][2];
              rjk2 = dxjk * dxjk + dyjk * dyjk + dzjk * dzjk;
              rjk = sqrt( rjk2 );
              if (rjk2 <= rbound) {
                dxik = delki[0]=x[k][0] - x[i][0];
                dyik = delki[1]=x[k][1] - x[i][1];
                dzik = delki[2]=x[k][2] - x[i][2];
                rik2 = dxik * dxik + dyik * dyik + dzik * dzik;
                rik = sqrt( rik2 );
                if (rik2 <= rbound) {
                  xik = rik2 / rij2;
                  xjk = rjk2 / rij2;
                  dxik=-2*rik2 /rij2/rij;
                  dxjk = -2*rjk2 / rij2/rij;
                  a = 1 - (xik - xjk) * (xik - xjk);
                  da = -2*(xik - xjk) * (dxik - dxjk);
                  if (!iszero(a)) {
                    cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
                    dcikj = (2.0 * (dxik + dxjk) + da) / a +(2.0 * (xik + xjk) + a - 2.0) *(-da)/ a / a;

                    if (cikj >= Cmin && cikj <= Cmax) {
                      assert(!iszero(delc));
                      cikj = (cikj - Cmin) / delc;
                      dcikj /= delc;

                      sikj = dfcut(cikj, dfc, ddfc);
                      dsikj = dfc * dcikj;
                      dCfunc2(rij2, rik2, rjk2, dCikj1, dCikj2); //--- 4.17b/rik, 4.17c/rjk
                      ddCfunc2(rij2, rik2, rjk2, ddCikj1, ddCikj2 );
                      dsij=ds;
                      assert(!iszero(sikj));
                      a = sij / delc * dfc / sikj;
                      da = (dsij*dfc/sikj+sij*ddfc*dcikj/sikj-dsikj*sij*dfc/sikj/sikj)/delc;

                      dsij1 = a * dCikj1; //--- 4.22b/rik: units of s/r^2
                      dsij2 = a * dCikj2; //--- 4.22c/rjk
                       
                      ddsij1drij = da * dCikj1+a * ddCikj1; //--- units of s/r^3
                      ddsij2drij = da * dCikj2+a * ddCikj2; //--- units of s/r^3//                    
                       assert(!isnan(da));
                       assert(!isnan(a));
                       assert(!isnan(dCikj2));
                       assert(!isnan(ddCikj2));

                      dsg_alpha_beta_drjk = recip * dUdsij * ddsij2drij * delij[alpha] * delij[beta];
                      dsg_alpha_beta_drik = recip * dUdsij * ddsij1drij * delij[alpha] * delij[beta];
                       assert(!isnan(recip));
                       assert(!isnan(dUdsij));
                       assert(!isnan(ddsij2drij));
                       assert(!isnan(delij[alpha]));
                       assert(!isnan(delij[beta]));
                    }
                  }
                }
              }
            }

              //
              //     Tabulate per-atom virial as symmetrized stress tensor
            if (!iszero(dsij1) || !iszero(dsij2) || !iszero(dsg_alpha_beta_drjk) || !iszero(dsg_alpha_beta_drik) ){ //modify!!!!!!!
               assert( !isnan(dsg_alpha_beta_drjk));
               assert( !isnan(dsg_alpha_beta_ds));
               assert( !isnan(dsij2));
               assert( !isnan(dsg_alpha_beta_drik));
               assert( !isnan(dsij1));
              mod3bdy += (dsg_alpha_beta_drjk + dsg_alpha_beta_ds * dsij2) * deljk[gamma] * deljk[lambda]+
                         (dsg_alpha_beta_drik + dsg_alpha_beta_ds * dsij1) * delki[gamma] * delki[lambda];
               
            }
          }
          //     end of k loop
        }
         return mod2bdy + mod3bdy;
     
 };


