#include "meam.h"
#include <cmath>
#include <algorithm>
#include "math_special.h"

using namespace LAMMPS_NS;


void
MEAM::meam_force(int i, int eflag_either, int eflag_global, int eflag_atom, int vflag_atom, double* eng_vdwl,
                 double* eatom, int /*ntype*/, int* type, int* fmap, double** scale, double** x, int numneigh, int* firstneigh,
                 int numneigh_full, int* firstneigh_full, int fnoffset, double** f, double** vatom)
{
  int j, jn, k, kn, kk, m, n, p, q, ii;
  int nv2, nv3, elti, eltj, eltk, ind;
  double xitmp, yitmp, zitmp, delij[3], delji[3], rij2, rij, rij3;
  double v[6], fi[3], fj[3], vm[21];
  double third, sixth;
  double pp, dUdrij, dUdsij, dUdrijm[3], force, forcem;
  double ddUddrij;
  double recip, phi, phip, phipp;
  double sij;
  double a1, a1i, a1j, a2, a2i, a2j;
  double a3i, a3j;
  double shpi[3], shpj[3];
  double ai, aj, ro0i, ro0j, invrei, invrej;
  double rhoa0j, drhoa0j, ddrhoa0j, rhoa0i, drhoa0i, ddrhoa0i;
  double rhoa1j, drhoa1j, ddrhoa1j, rhoa1i, drhoa1i, ddrhoa1i;
  double rhoa2j, drhoa2j, ddrhoa2j, rhoa2i, drhoa2i, ddrhoa2i;
  double rhoa3j, drhoa3j, ddrhoa3j, rhoa3i, drhoa3i, ddrhoa3i;
  double a3, a3a;
  double drho0dr1, drho0dr2, drho0ds1, drho0ds2;
  double drho0drm1[3], drho0drm2[3];
  double ddrho0drmdr1[3], ddrho0drmdr2[3];
  double ddrho0drdr1, ddrho0drdr2, ddrho1drdr1, ddrho1drdr2, ddrho2drdr1, ddrho2drdr2, ddrho3drdr1, ddrho3drdr2;
  double drho1dr1, drho1dr2, drho1ds1, drho1ds2;
  double drho1drm1[3], drho1drm2[3];
  double ddrho1drmdr1[3], ddrho1drmdr2[3];
  double drho2dr1, drho2dr2, drho2ds1, drho2ds2;
  double drho2drm1[3], drho2drm2[3];
  double ddrho2drmdr1[3], ddrho2drmdr2[3];
  double ddrho2drmdrn1[6], ddrho2drmdrn2[6];
  double drho3dr1, drho3dr2, drho3ds1, drho3ds2;
  double drho3drm1[3], drho3drm2[3];
  double ddrho3drmdr1[3], ddrho3drmdr2[3];
  double dt1dr1, dt1dr2, dt1ds1, dt1ds2;
  double dt2dr1, dt2dr2, dt2ds1, dt2ds2;
  double dt3dr1, dt3dr2, dt3ds1, dt3ds2;
  double drhodr1, drhodr2, drhods1, drhods2, drhodrm1[3], drhodrm2[3];
  double ddrhodrdr1, ddrhodrdr2, ddrhodrmdr1[3], ddrhodrmdr2[3];
  double arg;
  double arg1i1, arg1j1, arg1i2, arg1j2, arg1i3, arg1j3, arg3i3, arg3j3;
  double arg1i1_d, arg1j1_d, arg1i2_d, arg1j2_d, arg1i3_d, arg1j3_d, arg3i3_d, arg3j3_d;
  double dsij1, dsij2, force1, force2;
  double t1i, t2i, t3i, t1j, t2j, t3j;
  double scaleij;
  double ddrho3drmdrn1[6], ddrho3drmdrn2[6], ddrhodrmdrn1[6], ddrhodrmdrn2[6];
  double ddUdrdrijm[3], ddUdrijmdrijn[6];
  double stiff, stiff0, stiff1, stiff2;
  double n0, n1, n2;
  double ddt1drdr1,  ddt2drdr1,  ddt3drdr1;
  double ddt1drdr2,  ddt2drdr2,  ddt3drdr2;
  double ddrho1drmdrn1[6], ddrho1drmdrn2[6];
  
  third = 1.0 / 3.0;
  sixth = 1.0 / 6.0;

  //     Compute forces atom i

  elti = fmap[type[i]];
  if (elti < 0) return;

  xitmp = x[i][0];
  yitmp = x[i][1];
  zitmp = x[i][2];

  //     Treat each pair
  for (jn = 0; jn < numneigh; jn++) {
    j = firstneigh[jn];
    eltj = fmap[type[j]];
    scaleij = scale[type[i]][type[j]];

    if (!iszero(scrfcn[fnoffset + jn]) && eltj >= 0) {

      sij = scrfcn[fnoffset + jn] * fcpair[fnoffset + jn];
      delij[0] = x[j][0] - xitmp;
      delij[1] = x[j][1] - yitmp;
      delij[2] = x[j][2] - zitmp;
      delji[ 0 ] = -delij[ 0 ];
      delji[ 1 ] = -delij[ 1 ];
      delji[ 2 ] = -delij[ 2 ];
      rij2 = delij[0] * delij[0] + delij[1] * delij[1] + delij[2] * delij[2];
      if (rij2 < this->cutforcesq) {
        rij = sqrt(rij2);
        recip = 1.0 / rij;

        //     Compute phi and phip (potential function)
        ind = this->eltind[elti][eltj];
        pp = rij * this->rdrar;
        kk = (int)pp;
        kk = std::min(kk, this->nrar - 2);
        pp = pp - kk;
        pp = std::min(pp, 1.0);
        phi = ((this->phirar3[ind][kk] * pp + this->phirar2[ind][kk]) * pp + this->phirar1[ind][kk]) * pp + this->phirar[ind][kk]; //--- additional terms from the smoothing function
        phip = (this->phirar6[ind][kk] * pp + this->phirar5[ind][kk]) * pp + this->phirar4[ind][kk]; //--- (d/dr){\phi/S_{ij}}: polynomial smoothing function
        phipp = (this->phirar8[ind][kk]) * pp + this->phirar7[ind][kk]; //--- (d^2/dr^2){\phi/S_{ij}}

        if (eflag_either != 0) {
          double phi_sc = phi * scaleij; //--- scaled energy: scaleij = 1/zij0
          if (eflag_global != 0)
            *eng_vdwl = *eng_vdwl + phi_sc * sij;
          if (eflag_atom != 0) {
            eatom[i] = eatom[i] + 0.5 * phi_sc * sij;  //--- pair-wise term
            eatom[j] = eatom[j] + 0.5 * phi_sc * sij;
          }
        }

        //     write(1,*) "force_meamf: phi: ",phi
        //     write(1,*) "force_meamf: phip: ",phip

        //     Compute pair densities and derivatives
        invrei = 1.0 / this->re_meam[elti][elti];
        ai = rij * invrei - 1.0; //--- dimensionless distance
        ro0i = this->rho0_meam[elti];
        rhoa0i = ro0i * MathSpecial::fm_exp(-this->beta0_meam[elti] * ai); //--- Eq. (4.8)
        drhoa0i = -this->beta0_meam[elti] * invrei * rhoa0i; //--- drho/drij 
        ddrhoa0i = -this->beta0_meam[elti] * invrei * drhoa0i; //--- d^2rho/drij^2 
        rhoa1i = ro0i * MathSpecial::fm_exp(-this->beta1_meam[elti] * ai);
        drhoa1i = -this->beta1_meam[elti] * invrei * rhoa1i;
        ddrhoa1i = -this->beta1_meam[elti] * invrei * drhoa1i;
        rhoa2i = ro0i * MathSpecial::fm_exp(-this->beta2_meam[elti] * ai);
        drhoa2i = -this->beta2_meam[elti] * invrei * rhoa2i;
        ddrhoa2i = -this->beta2_meam[elti] * invrei * drhoa2i;
        rhoa3i = ro0i * MathSpecial::fm_exp(-this->beta3_meam[elti] * ai);
        drhoa3i = -this->beta3_meam[elti] * invrei * rhoa3i;
        ddrhoa3i = -this->beta3_meam[elti] * invrei * drhoa3i;

        if (elti != eltj) {
          invrej = 1.0 / this->re_meam[eltj][eltj];
          aj = rij * invrej - 1.0;
          ro0j = this->rho0_meam[eltj];
          rhoa0j = ro0j * MathSpecial::fm_exp(-this->beta0_meam[eltj] * aj);
          drhoa0j = -this->beta0_meam[eltj] * invrej * rhoa0j;
          ddrhoa0j = -this->beta0_meam[eltj] * invrej * drhoa0j;
          rhoa1j = ro0j * MathSpecial::fm_exp(-this->beta1_meam[eltj] * aj);
          drhoa1j = -this->beta1_meam[eltj] * invrej * rhoa1j;
          ddrhoa1j = -this->beta1_meam[eltj] * invrej * drhoa1j;
          rhoa2j = ro0j * MathSpecial::fm_exp(-this->beta2_meam[eltj] * aj);
          drhoa2j = -this->beta2_meam[eltj] * invrej * rhoa2j;
          ddrhoa2j = -this->beta2_meam[eltj] * invrej * drhoa2j;
          rhoa3j = ro0j * MathSpecial::fm_exp(-this->beta3_meam[eltj] * aj);
          drhoa3j = -this->beta3_meam[eltj] * invrej * rhoa3j;
          ddrhoa3j = -this->beta3_meam[eltj] * invrej * drhoa3j;
        } else {
          rhoa0j = rhoa0i;
          drhoa0j = drhoa0i;
          ddrhoa0j = ddrhoa0i;
          rhoa1j = rhoa1i;
          drhoa1j = drhoa1i;
          ddrhoa1j = ddrhoa1i;
          rhoa2j = rhoa2i;
          drhoa2j = drhoa2i;
          ddrhoa2j = ddrhoa2i;
          rhoa3j = rhoa3i;
          drhoa3j = drhoa3i;
          ddrhoa3j = ddrhoa3i;
        }

        const double t1mi = this->t1_meam[elti];
        const double t2mi = this->t2_meam[elti];
        const double t3mi = this->t3_meam[elti];
        const double t1mj = this->t1_meam[eltj];
        const double t2mj = this->t2_meam[eltj];
        const double t3mj = this->t3_meam[eltj];

        if (this->ialloy == 1) {
          rhoa1j  *= t1mj;
          rhoa2j  *= t2mj;
          rhoa3j  *= t3mj;
          rhoa1i  *= t1mi;
          rhoa2i  *= t2mi;
          rhoa3i  *= t3mi;
          drhoa1j *= t1mj;
          drhoa2j *= t2mj;
          drhoa3j *= t3mj;
          drhoa1i *= t1mi;
          drhoa2i *= t2mi;
          drhoa3i *= t3mi;
          ddrhoa1j *= t1mj;
          ddrhoa2j *= t2mj;
          ddrhoa3j *= t3mj;
          ddrhoa1i *= t1mi;
          ddrhoa2i *= t2mi;
          ddrhoa3i *= t3mi;
        }

        nv2 = 0;
        nv3 = 0;
        arg1i1 = 0.0;
        arg1j1 = 0.0;
        arg1i2 = 0.0;
        arg1j2 = 0.0;
        arg1i3 = 0.0;
        arg1j3 = 0.0;
        arg3i3 = 0.0;
        arg3j3 = 0.0;
        for (n = 0; n < 3; n++) {
          for (p = n; p < 3; p++) {
            for (q = p; q < 3; q++) {
              arg = delij[n] * delij[p] * delij[q] * this->v3D[nv3];
              arg1i3 +=  arho3[i][nv3] * arg; //--- arho3 is Y_{3i\sigma\beta\gamma} Eq.(4.27c)
              arg1j3 +=  - arho3[j][nv3] * arg; 
              arg1i3_d +=  darho3dr[i][nv3] * arg;
              arg1j3_d +=  - darho3dr[j][nv3] * arg;
              nv3 +=  1;
            }
            arg = delij[n] * delij[p] * this->v2D[nv2];
            arg1i2 +=  arho2[i][nv2] * arg;
            arg1i2_d +=  darho2dr[i][nv2] * arg;
            arg1j2 +=  arho2[j][nv2] * arg;
            arg1j2_d += darho2dr[j][nv2] * arg;
            nv2 +=  1;
          }
          arg1i1 += arho1[i][n] * delij[n]; //--- 4.30(a) in sandia report:  arho1[i][n] is Y_{1i\sigma}
          arg1j1 += - arho1[j][n] * delij[n];
          arg3i3 += arho3b[i][n] * delij[n];
          arg3j3 += - arho3b[j][n] * delij[n];
          arg1i1_d +=  darho1dr[i][n] * delij[n]; 
          arg3i3_d +=  darho3bdr[i][n] * delij[n];          
          arg1j1_d +=  - darho1dr[i][n] * delij[n]; 
          arg3j3_d +=  - darho3bdr[j][n] * delij[n];
        }

        //     rho0 terms
        drho0dr1 = drhoa0j * sij; //--- 4.26(a)
        drho0dr2 = drhoa0i * sij;
        ddrho0drdr1 = ddrhoa0j * sij; 
        ddrho0drdr2 = ddrhoa0i * sij;
        for (m = 0; m < 3; m++) {
          drho0drm1[m] = 0.0; 
          drho0drm2[m] = 0.0; 
          ddrho0drmdr1[m] = 0.0; 
          ddrho0drmdr2[m] = 0.0; 
          
        }
        //     rho1 terms
        a1 = 2 * sij / rij;
        drho1dr1 = a1 * (drhoa1j - rhoa1j / rij) * arg1i1; //--- 4.30(a)
        drho1dr2 = a1 * (drhoa1i - rhoa1i / rij) * arg1j1;
        ddrho1drdr1 = Get_ddrho1drdr( i, //--- deriv of 4.30(a) wrt rij
                                     rij,  sij, 
                                     rhoa1j,  drhoa1j,  ddrhoa1j,
                                     arho2b,
                                     arg1i1,
                                     arg1i1_d
                    );
        ddrho1drdr2 = Get_ddrho1drdr( j, 
                                     rij,  sij, 
                                     rhoa1i,  drhoa1i,  ddrhoa1i,
                                     arho2b,
                                     arg1j1,
                                     arg1j1_d
                    )  ;
        a1 = 2.0 * sij / rij;
        for (m = 0; m < 3; m++) {
          drho1drm1[m] = a1 * rhoa1j * arho1[i][m]; //--- 4.30(c)
          drho1drm2[m] = -a1 * rhoa1i * arho1[j][m]; //--- negative sign??? should be in arho1[j]??
        }
        //
        Get_ddrho1drmdr( i, //--- deriv of 4.30(c) wrt r
                        rij,  sij,  delij,
                        rhoa1j,  drhoa1j,
                        arho1,  darho1dr,
                        ddrho1drmdr1 //--- modify 
                    );
        Get_ddrho1drmdr( j,
                        rij,  sij,  delji,
                        rhoa1i,  drhoa1i,
                        arho1,  darho1dr,
                        ddrho1drmdr2 //--- modify 
                    );
         // 
        Get_ddrho1drmdrn( i, //--- deriv of 4.30(c) wrt rn
                        rij,  sij,
                        rhoa1j,
                        ddrho1drmdrn1 //--- modify 
                    );
        Get_ddrho1drmdrn( j,
                        rij,  sij,
                        rhoa1i,
                        ddrho1drmdrn2 //--- modify 
                    );
          
        //     rho2 terms
        a2 = 2 * sij / rij2;
        drho2dr1 = a2 * (drhoa2j - 2 * rhoa2j / rij) * arg1i2 - 2.0 / 3.0 * arho2b[i] * drhoa2j * sij; //--- 4.30(d): arho2b is W_{2i}      
        drho2dr2 = a2 * (drhoa2i - 2 * rhoa2i / rij) * arg1j2 - 2.0 / 3.0 * arho2b[j] * drhoa2i * sij;
        //--- 2nd derivative wrt rij (atom j)
        ddrho2drdr1 = Get_ddrho2drdr( i,  //--- deriv. of 4.30(d) wrt r
                     rij,  sij, 
                     rhoa2j,  drhoa2j,  ddrhoa2j,
                     arho2b,
                     arg1i2,
                     arg1i2_d
                    );
        ddrho2drdr2 = Get_ddrho2drdr( j, 
                     rij,  sij, 
                     rhoa2i,  drhoa2i,  ddrhoa2i,
                     arho2b,
                     arg1j2,
                     arg1j2_d
                    );
        //
        a2 = 4 * sij / rij2;
        for (m = 0; m < 3; m++) {
          drho2drm1[m] = 0.0;
          drho2drm2[m] = 0.0;
          for (n = 0; n < 3; n++) {
            drho2drm1[m] = drho2drm1[m] + arho2[i][this->vind2D[m][n]] * delij[n]; //--- 4.30(f): arho2 is Y_{2i\sigma\alpha}
            drho2drm2[m] = drho2drm2[m] - arho2[j][this->vind2D[m][n]] * delij[n];
          }
          //
          drho2drm1[m] = a2 * rhoa2j * drho2drm1[m];
          drho2drm2[m] = -a2 * rhoa2i * drho2drm2[m];
        }
        Get_ddrho2drmdr( i, //--- deriv of 4.30(f) wrt r
                         rij,  sij, delij,
                        rhoa2j,  drhoa2j,
                        arho2,
                        darho2dr,
                        ddrho2drmdr1 //--- modify 
                    );
         Get_ddrho2drmdr( j,
                         rij,  sij, delji,
                        rhoa2i,  drhoa2i,
                        arho2,
                        darho2dr,
                        ddrho2drmdr2 //--- modify
                    );
        //
        Get_ddrho2drmdrn( i, //--- deriv of 4.30(f) wrt rn
                         rij,  sij, delij,
                         rhoa2j,
                         arho2,
                         ddrho2drmdrn1 //--- modify 
                         );
        Get_ddrho2drmdrn( j, //--- deriv of 4.30(f) wrt rn
                         rij,  sij, delji,
                         rhoa2i,
                         arho2,
                         ddrho2drmdrn2 //--- modify 
                        )  ;      
        //
        //     rho3 terms
        rij3 = rij * rij2;
        a3 = 2 * sij / rij3;
        a3a = 6.0 / 5.0 * sij / rij;
        drho3dr1 = a3 * (drhoa3j - 3 * rhoa3j / rij) * arg1i3 - a3a * (drhoa3j - rhoa3j / rij) * arg3i3; //--- 4.30(g)
        drho3dr2 = a3 * (drhoa3i - 3 * rhoa3i / rij) * arg1j3 - a3a * (drhoa3i - rhoa3i / rij) * arg3j3;
        ddrho3drdr1 = Get_ddrho3drdr( rij, sij, //--- deriv. of 4.30(g) wrt r   where it was used???
                                    rhoa3j,  drhoa3j,  ddrhoa3j,
                                    arg1i3,  arg3i3,
                                    arg1i3_d, arg3i3_d
                                   );
        ddrho3drdr2 = Get_ddrho3drdr( rij, sij,
                                    rhoa3i,  drhoa3i,  ddrhoa3i,
                                    arg1j3,  arg3j3,
                                    arg1j3_d, arg3j3_d
                                   )   ;       
          
        a3 = 6 * sij / rij3;
        a3a = 6 * sij / (5 * rij);
        for (m = 0; m < 3; m++) {
          drho3drm1[m] = 0.0;
          drho3drm2[m] = 0.0;
          //
          nv2 = 0;
          for (n = 0; n < 3; n++) {
            for (p = n; p < 3; p++) {
              arg = delij[n] * delij[p] * this->v2D[nv2];
              drho3drm1[m] = drho3drm1[m] + arho3[i][this->vind3D[m][n][p]] * arg; //--- 4.30(i)
              drho3drm2[m] = drho3drm2[m] + arho3[j][this->vind3D[m][n][p]] * arg;
              nv2 = nv2 + 1;
            }
          }
          drho3drm1[m] = (a3 * drho3drm1[m] - a3a * arho3b[i][m]) * rhoa3j;
          drho3drm2[m] = (-a3 * drho3drm2[m] + a3a * arho3b[j][m]) * rhoa3i; //--- negative sign??? 
        }
        Get_ddrho3drmdr( i,
                         rij,  sij, delij,
                         rhoa3j,
                         drhoa3j,
                         darho3dr,
                         ddrho3drmdr1); //--- modify ddrho3drmdr1[m] where it was used????
        Get_ddrho3drmdr( j,
                         rij,  sij, delji,
                         rhoa3i,
                         drhoa3i,
                         darho3dr,
                         ddrho3drmdr2); //--- modify ddrho3drmdr2[m]
        //
        Get_ddrho3drmdrn( i, //--- deriv. of 4.30(i) wrt rn
                         rij,  sij, delij,
                         rhoa3j, 
                         ddrho3drmdrn1);
        Get_ddrho3drmdrn( j,
                         rij,  sij, delji,
                         rhoa3i, 
                         ddrho3drmdrn2);
        
        //     Compute derivatives of weighting functions t wrt rij
        t1i = t_ave[i][0];
        t2i = t_ave[i][1];
        t3i = t_ave[i][2];
        t1j = t_ave[j][0];
        t2j = t_ave[j][1];
        t3j = t_ave[j][2];

        if (this->ialloy == 1) {   //--- not included in the report ?????????

          a1i = fdiv_zero(drhoa0j * sij, tsq_ave[i][0]); 
          a1j = fdiv_zero(drhoa0i * sij, tsq_ave[j][0]);
          a2i = fdiv_zero(drhoa0j * sij, tsq_ave[i][1]);
          a2j = fdiv_zero(drhoa0i * sij, tsq_ave[j][1]);
          a3i = fdiv_zero(drhoa0j * sij, tsq_ave[i][2]);
          a3j = fdiv_zero(drhoa0i * sij, tsq_ave[j][2]);

          dt1dr1 = a1i * (t1mj - t1i * MathSpecial::square(t1mj));
          dt1dr2 = a1j * (t1mi - t1j * MathSpecial::square(t1mi));
          dt2dr1 = a2i * (t2mj - t2i * MathSpecial::square(t2mj));
          dt2dr2 = a2j * (t2mi - t2j * MathSpecial::square(t2mi));
          dt3dr1 = a3i * (t3mj - t3i * MathSpecial::square(t3mj));
          dt3dr2 = a3j * (t3mi - t3j * MathSpecial::square(t3mi));

        } else if (this->ialloy == 2) {

          dt1dr1 = 0.0;
          dt1dr2 = 0.0;
          dt2dr1 = 0.0;
          dt2dr2 = 0.0;
          dt3dr1 = 0.0;
          dt3dr2 = 0.0;

        } else {

          ai = 0.0;
          if (!iszero(rho0[i]))
            ai = drhoa0j * sij / rho0[i];
          aj = 0.0;
          if (!iszero(rho0[j]))
            aj = drhoa0i * sij / rho0[j];

          dt1dr1 = ai * (t1mj - t1i); //--- 4.32(a)
          dt1dr2 = aj * (t1mi - t1j);
          dt2dr1 = ai * (t2mj - t2i);
          dt2dr2 = aj * (t2mi - t2j);
          dt3dr1 = ai * (t3mj - t3i);
          dt3dr2 = aj * (t3mi - t3j);
          
          ddt1drdr1 = sij * ( - dt1dr1 * drhoa0j + ( t1mj - t1i ) * ddrhoa0j ) - dt1dr1 * drho0dr1; //--- deriv of 4.32(a) wrt. r
          ddt1drdr1 /= rho0[i];
          ddt2drdr1 = sij * ( - dt2dr1 * drhoa0j + ( t2mj - t2i ) * ddrhoa0j ) - dt2dr1 * drho0dr1;
          ddt2drdr1 /= rho0[i];
          ddt3drdr1 = sij * ( - dt3dr1 * drhoa0j + ( t3mj - t3i ) * ddrhoa0j ) - dt3dr1 * drho0dr1;
          ddt3drdr1 /= rho0[i];
            
          ddt1drdr2 = sij * ( - dt1dr2 * drhoa0i + ( t1mi - t1j ) * ddrhoa0i ) - dt1dr2 * drho0dr2; //--- index j
          ddt1drdr2 /= rho0[j];
          ddt2drdr2 = sij * ( - dt2dr2 * drhoa0i + ( t2mi - t2j ) * ddrhoa0i ) - dt2dr2 * drho0dr2;
          ddt2drdr2 /= rho0[j];
          ddt3drdr2 = sij * ( - dt3dr2 * drhoa0i + ( t3mi - t3j ) * ddrhoa0i ) - dt3dr2 * drho0dr2;
          ddt3drdr2 /= rho0[j] ;                     
        }
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        //------- Compute derivatives of total density wrt rij, sij and rij(3)  
        //------- total density is weighted average of zero and higher order densities
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        get_shpfcn(this->lattce_meam[elti][elti], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shpi);
        get_shpfcn(this->lattce_meam[eltj][eltj], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shpj);

        drhodr1 = dgamma1[i] * drho0dr1 + //--- index i: Eq. 4.36(a), dgamma1: defined in "meam_dens_final.cpp"
          dgamma2[i] * (dt1dr1 * rho1[i] + t1i * drho1dr1 + dt2dr1 * rho2[i] + t2i * drho2dr1 +
                        dt3dr1 * rho3[i] + t3i * drho3dr1) -
          dgamma3[i] * (shpi[0] * dt1dr1 + shpi[1] * dt2dr1 + shpi[2] * dt3dr1); 
        drhodr2 = dgamma1[j] * drho0dr2 +
          dgamma2[j] * (dt1dr2 * rho1[j] + t1j * drho1dr2 + dt2dr2 * rho2[j] + t2j * drho2dr2 +
                        dt3dr2 * rho3[j] + t3j * drho3dr2) -
          dgamma3[j] * (shpj[0] * dt1dr2 + shpj[1] * dt2dr2 + shpj[2] * dt3dr2);
        //---                
        ddrhodrdr1 = Get_ddrhodrdr(i, elti, //--- deriv. of Eq. 4.36(a) wrt. r
                                    shpi, 
                                    t1i,  t2i,  t3i,
                                    dt1dr1,  dt2dr1,  dt3dr1,
                                    ddt1drdr1,  ddt2drdr1,  ddt3drdr1,
                                    rho0, rho1, rho2, rho3, 
                                    drho0dr1,  drho1dr1,  drho2dr1,  drho3dr1, 
                                    ddrho0drdr1, ddrho1drdr1,  ddrho2drdr1,  ddrho3drdr1,
                                    drhodr1
                                  );
        //--- index j
        ddrhodrdr2 = Get_ddrhodrdr(j, eltj,
                                   shpj, 
                                   t1j,  t2j,  t3j,
                                   dt1dr2,  dt2dr2,  dt3dr2,
                                   ddt1drdr2,  ddt2drdr2,  ddt3drdr2,
                                   rho0, rho1, rho2, rho3, 
                                   drho0dr2,  drho1dr2,  drho2dr2,  drho3dr2, 
                                   ddrho0drdr2, ddrho1drdr2,  ddrho2drdr2,  ddrho3drdr2,
                                   drhodr2
                                  );
        
        //--- deriv. wrt. rij(3)
        for (m = 0; m < 3; m++) {
          drhodrm1[m] = 0.0;
          drhodrm2[m] = 0.0;
          drhodrm1[m] = dgamma2[i] * (t1i * drho1drm1[m] + t2i * drho2drm1[m] + t3i * drho3drm1[m]); //--- index i: Eq. 4.36(c)
          drhodrm2[m] = dgamma2[j] * (t1j * drho1drm2[m] + t2j * drho2drm2[m] + t3j * drho3drm2[m]);
        }
        Get_ddrhodrmdr( i, elti, //--- deriv. of Eq. 4.36(c) wrt. r
                        shpi, 
                        t1i,  t2i,  t3i,
                        dt1dr1,  dt2dr1,  dt3dr1,
                        rho0, rho1, rho2, rho3,
                        drho0dr1,  drho1dr1,  drho2dr1,  drho3dr1, 
                        drho0drm1,  drho1drm1,  drho2drm1,  drho3drm1, 
                        ddrho0drmdr1, ddrho1drmdr1,  ddrho2drmdr1,  ddrho3drmdr1,
                        drhodrm1,
                        ddrhodrmdr1 //--- modify
                       );
        Get_ddrhodrmdr( j, eltj, //--- deriv. of Eq. 4.36(c) wrt. r
                        shpj, 
                        t1j,  t2j,  t3j,
                        dt1dr2,  dt2dr2,  dt3dr2,
                        rho0, rho1, rho2, rho3,
                        drho0dr2,  drho1dr2,  drho2dr2,  drho3dr2, 
                        drho0drm2,  drho1drm2,  drho2drm2,  drho3drm2, 
                        ddrho0drmdr2, ddrho1drmdr2,  ddrho2drmdr2,  ddrho3drmdr2,
                        drhodrm2,
                        ddrhodrmdr2 //--- modify
                       );
         
        Get_ddrhodrmdrn(  i,  elti, //--- deriv. of Eq. 4.36(c) wrt. rm
                         shpi, 
                         t1i,   t2i,   t3i,
                         rho0,  rho1,  rho2,  rho3,
                         drho0drm1,   drho1drm1,  drho2drm1,  drho3drm1, 
                                      ddrho1drmdrn1,  ddrho2drmdrn1, ddrho3drmdrn1,
                         ddrhodrmdrn1 //--- modify
                       );
        Get_ddrhodrmdrn(  j,  eltj, //--- deriv. of Eq. 4.36(c) wrt. rm
                         shpj, 
                         t1j,   t2j,   t3j,
                         rho0,  rho1,  rho2,  rho3,
                         drho0drm2,   drho1drm2,  drho2drm2,  drho3drm2, 
                                      ddrho1drmdrn2,  ddrho2drmdrn2, ddrho3drmdrn2,
                         ddrhodrmdrn2 //--- modify
                       );
          
        //     Compute derivatives wrt sij, but only if necessary     wrt s??????????????/
        if (!iszero(dscrfcn[fnoffset + jn])) {
          drho0ds1 = rhoa0j;
          drho0ds2 = rhoa0i;
          a1 = 2.0 / rij;
          drho1ds1 = a1 * rhoa1j * arg1i1;
          drho1ds2 = a1 * rhoa1i * arg1j1;
          a2 = 2.0 / rij2;
          drho2ds1 = a2 * rhoa2j * arg1i2 - 2.0 / 3.0 * arho2b[i] * rhoa2j;
          drho2ds2 = a2 * rhoa2i * arg1j2 - 2.0 / 3.0 * arho2b[j] * rhoa2i;
          a3 = 2.0 / rij3;
          a3a = 6.0 / (5.0 * rij);
          drho3ds1 = a3 * rhoa3j * arg1i3 - a3a * rhoa3j * arg3i3;
          drho3ds2 = a3 * rhoa3i * arg1j3 - a3a * rhoa3i * arg3j3;

          if (this->ialloy == 1) {
            a1i = fdiv_zero(rhoa0j, tsq_ave[i][0]);
            a1j = fdiv_zero(rhoa0i, tsq_ave[j][0]);
            a2i = fdiv_zero(rhoa0j, tsq_ave[i][1]);
            a2j = fdiv_zero(rhoa0i, tsq_ave[j][1]);
            a3i = fdiv_zero(rhoa0j, tsq_ave[i][2]);
            a3j = fdiv_zero(rhoa0i, tsq_ave[j][2]);

            dt1ds1 = a1i * (t1mj - t1i * MathSpecial::square(t1mj));
            dt1ds2 = a1j * (t1mi - t1j * MathSpecial::square(t1mi));
            dt2ds1 = a2i * (t2mj - t2i * MathSpecial::square(t2mj));
            dt2ds2 = a2j * (t2mi - t2j * MathSpecial::square(t2mi));
            dt3ds1 = a3i * (t3mj - t3i * MathSpecial::square(t3mj));
            dt3ds2 = a3j * (t3mi - t3j * MathSpecial::square(t3mi));

          } else if (this->ialloy == 2) {

            dt1ds1 = 0.0;
            dt1ds2 = 0.0;
            dt2ds1 = 0.0;
            dt2ds2 = 0.0;
            dt3ds1 = 0.0;
            dt3ds2 = 0.0;

          } else {

            ai = 0.0;
            if (!iszero(rho0[i]))
              ai = rhoa0j / rho0[i];
            aj = 0.0;
            if (!iszero(rho0[j]))
              aj = rhoa0i / rho0[j];

            dt1ds1 = ai * (t1mj - t1i);
            dt1ds2 = aj * (t1mi - t1j);
            dt2ds1 = ai * (t2mj - t2i);
            dt2ds2 = aj * (t2mi - t2j);
            dt3ds1 = ai * (t3mj - t3i);
            dt3ds2 = aj * (t3mi - t3j);
          }

          drhods1 = dgamma1[i] * drho0ds1 +
            dgamma2[i] * (dt1ds1 * rho1[i] + t1i * drho1ds1 + dt2ds1 * rho2[i] + t2i * drho2ds1 +
                          dt3ds1 * rho3[i] + t3i * drho3ds1) -
            dgamma3[i] * (shpi[0] * dt1ds1 + shpi[1] * dt2ds1 + shpi[2] * dt3ds1);
          drhods2 = dgamma1[j] * drho0ds2 +
            dgamma2[j] * (dt1ds2 * rho1[j] + t1j * drho1ds2 + dt2ds2 * rho2[j] + t2j * drho2ds2 +
                          dt3ds2 * rho3[j] + t3j * drho3ds2) -
            dgamma3[j] * (shpj[0] * dt1ds2 + shpj[1] * dt2ds2 + shpj[2] * dt3ds2);
        }

        //     Compute derivatives of energy wrt rij, sij and rij[3]
        dUdrij = phip * sij + frhop[i] * drhodr1 + frhop[j] * drhodr2; //--- Eq. 4.41(a)
        ddUddrij = phipp * sij + ( frhopp[i] * drhodr1 * drhodr1 + frhop[i] * ddrhodrdr1 ) + //--- 1st deriv. of Eq. 4.41(a) wrt r
                                 ( frhopp[j] * drhodr2 * drhodr2 + frhop[j] * ddrhodrdr2 );
        dUdsij = 0.0;
        if (!iszero(dscrfcn[fnoffset + jn])) {
          dUdsij = phi + frhop[i] * drhods1 + frhop[j] * drhods2; //--- Eq. 4.41(b)
        }
        nv2 = 0;
        for (m = 0; m < 3; m++) {
          dUdrijm[m] = frhop[i] * drhodrm1[m] + frhop[j] * drhodrm2[m]; //--- Eq. 4.41(c)
          ddUdrdrijm[m] = frhopp[i] * drhodr1 * drhodrm1[m] + frhop[i] * ddrhodrmdr1[m] + 
                          frhopp[j] * drhodr2 * drhodrm2[m] + frhop[i] * ddrhodrmdr2[m]; //--- deriv of Eq. 4.41(c) wrt r
          for (n = 0; n < 3; n++) {
            ddUdrijmdrijn[nv2] += frhopp[i] * drhodrm1[m] * drhodrm1[n] + frhop[i] * ddrhodrmdrn1[nv2]+
                                   frhopp[j] * drhodrm2[m] * drhodrm2[n] + frhop[j] * ddrhodrmdrn2[nv2];
            nv2++;
          }  
        }
        if (!isone(scaleij)) {
          dUdrij *= scaleij;
          dUdsij *= scaleij;
          dUdrijm[0] *= scaleij;
          dUdrijm[1] *= scaleij;
          dUdrijm[2] *= scaleij;
        }

        //     Add the part of the force due to dUdrij and dUdsij

        force = dUdrij * recip + dUdsij * dscrfcn[fnoffset + jn]; //-- recip = 1/r_{ij}
        for (m = 0; m < 3; m++) {
          forcem = delij[m] * force + dUdrijm[m]; //--- Eq. (4.40)
          f[i][m] = f[i][m] + forcem;
          f[j][m] = f[j][m] - forcem;
        }

        //--- add stiffness (units of u/r^2)
        stiff = ddUddrij - dUdrij * recip;
        stiff0 = 0.0; 
        stiff1 = 0.0;
        nv2 = 0;
        for (m = 0; m < 3; m++) {
          stiff0 += - dUdrijm[m] * delij[m];
          stiff1 += ddUdrdrijm[m] * delij[m];
          for (n = 0; n < 3; n++) {
            stiff2 += ddUdrijmdrijn[nv2] * delij[m] * delij[n];
            nv2++;
          }
        }
        stiff += ( ( stiff0 + stiff2 ) * recip  + stiff1 ) * recip;
        
        
        //     Tabulate per-atom virial as symmetrized stress tensor

        if (vflag_atom != 0) {
          fi[0] = delij[0] * force + dUdrijm[0];
          fi[1] = delij[1] * force + dUdrijm[1];
          fi[2] = delij[2] * force + dUdrijm[2];
          v[0] = -0.5 * (delij[0] * fi[0]);
          v[1] = -0.5 * (delij[1] * fi[1]);
          v[2] = -0.5 * (delij[2] * fi[2]);
          v[3] = -0.25 * (delij[0] * fi[1] + delij[1] * fi[0]);
          v[4] = -0.25 * (delij[0] * fi[2] + delij[2] * fi[0]);
          v[5] = -0.25 * (delij[1] * fi[2] + delij[2] * fi[1]);
          nv2 = 0;
          for (m = 0; m < 6; m++) {
            vatom[i][m] = vatom[i][m] + v[m];
            vatom[j][m] = vatom[j][m] + v[m];
            nv2++;
          }
          n0 = delij[0] * recip;
          n1 = delij[1] * recip;
          n2 = delij[2] * recip;
          //--- per-atom modulus
          vm[ 0 ]  = stiff * n0 * n0 * n0 * n0;
          vm[ 1 ]  = stiff * n0 * n0 * n1 * n1;
          vm[ 2 ]  = stiff * n0 * n0 * n2 * n2;
          vm[ 3 ]  = stiff * n0 * n0 * n0 * n1;
          vm[ 4 ]  = stiff * n0 * n0 * n0 * n2;
          vm[ 5 ]  = stiff * n0 * n0 * n1 * n2;
          //
          vm[ 6 ]  = stiff * n1 * n1 * n1 * n1;
          vm[ 7 ]  = stiff * n1 * n1 * n2 * n2;
          vm[ 8 ]  = stiff * n1 * n1 * n0 * n1;
          vm[ 9 ]  = stiff * n1 * n1 * n0 * n2;
          vm[ 10 ] = stiff * n1 * n1 * n1 * n2;
          //
          vm[ 11 ] = stiff * n2 * n2 * n2 * n2;
          vm[ 12 ] = stiff * n2 * n2 * n0 * n1;
          vm[ 13 ] = stiff * n2 * n2 * n0 * n2;
          vm[ 14 ] = stiff * n2 * n2 * n1 * n2;
          //
          vm[ 15 ] = stiff * n0 * n1 * n0 * n1;
          vm[ 16 ] = stiff * n0 * n1 * n0 * n2;
          vm[ 17 ] = stiff * n0 * n1 * n1 * n2;
          //
          vm[ 18 ] = stiff * n0 * n2 * n0 * n2;
          vm[ 19 ] = stiff * n0 * n2 * n1 * n2;
          //
          vm[ 20 ] = stiff * n1 * n2 * n1 * n2;
          //
          for (m = 0; m < 6; m++) {
            for (n = m; n < 6; n++) {
              vatom[i][nv2] += 0.5 * vm[nv2]; //satom defined?? analog to vatom //--- *r^2 to get energy  //defined?????
              vatom[j][nv2] += 0.5 * vm[nv2];
              nv2++;
            }
          }
        }

        //     Now compute forces on other atoms k due to change in sij     stiffness ??????????????

        if (iszero(sij) || isone(sij)) continue; //: cont jn loop

        double dxik(0), dyik(0), dzik(0);
        double dxjk(0), dyjk(0), dzjk(0);

        for (kn = 0; kn < numneigh_full; kn++) {
          k = firstneigh_full[kn];
          eltk = fmap[type[k]];
          if (k != j && eltk >= 0) {
            double xik, xjk, cikj, sikj, dfc, a;
            double dCikj1, dCikj2;
            double delc, rik2, rjk2;

            sij = scrfcn[jn+fnoffset] * fcpair[jn+fnoffset];
            const double Cmax = this->Cmax_meam[elti][eltj][eltk];
            const double Cmin = this->Cmin_meam[elti][eltj][eltk];

            dsij1 = 0.0;
            dsij2 = 0.0;
            if (!iszero(sij) && !isone(sij)) {
              const double rbound = rij2 * this->ebound_meam[elti][eltj];
              delc = Cmax - Cmin;
              dxjk = x[k][0] - x[j][0];
              dyjk = x[k][1] - x[j][1];
              dzjk = x[k][2] - x[j][2];
              rjk2 = dxjk * dxjk + dyjk * dyjk + dzjk * dzjk;
              if (rjk2 <= rbound) {
                dxik = x[k][0] - x[i][0];
                dyik = x[k][1] - x[i][1];
                dzik = x[k][2] - x[i][2];
                rik2 = dxik * dxik + dyik * dyik + dzik * dzik;
                if (rik2 <= rbound) {
                  xik = rik2 / rij2;
                  xjk = rjk2 / rij2;
                  a = 1 - (xik - xjk) * (xik - xjk);
                  if (!iszero(a)) {
                    cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
                    if (cikj >= Cmin && cikj <= Cmax) {
                      cikj = (cikj - Cmin) / delc;
                      sikj = dfcut(cikj, dfc);
                      dCfunc2(rij2, rik2, rjk2, dCikj1, dCikj2);
                      a = sij / delc * dfc / sikj;
                      dsij1 = a * dCikj1;
                      dsij2 = a * dCikj2;
                    }
                  }
                }
              }
            }

            if (!iszero(dsij1) || !iszero(dsij2)) {
              force1 = dUdsij * dsij1;
              force2 = dUdsij * dsij2;

              f[i][0] += force1 * dxik;
              f[i][1] += force1 * dyik;
              f[i][2] += force1 * dzik;
              f[j][0] += force2 * dxjk;
              f[j][1] += force2 * dyjk;
              f[j][2] += force2 * dzjk;
              f[k][0] -= force1 * dxik + force2 * dxjk;
              f[k][1] -= force1 * dyik + force2 * dyjk;
              f[k][2] -= force1 * dzik + force2 * dzjk;

              //     Tabulate per-atom virial as symmetrized stress tensor

              if (vflag_atom != 0) {
                fi[0] = force1 * dxik;
                fi[1] = force1 * dyik;
                fi[2] = force1 * dzik;
                fj[0] = force2 * dxjk;
                fj[1] = force2 * dyjk;
                fj[2] = force2 * dzjk;
                v[0] = -third * (dxik * fi[0] + dxjk * fj[0]);
                v[1] = -third * (dyik * fi[1] + dyjk * fj[1]);
                v[2] = -third * (dzik * fi[2] + dzjk * fj[2]);
                v[3] = -sixth * (dxik * fi[1] + dxjk * fj[1] + dyik * fi[0] + dyjk * fj[0]);
                v[4] = -sixth * (dxik * fi[2] + dxjk * fj[2] + dzik * fi[0] + dzjk * fj[0]);
                v[5] = -sixth * (dyik * fi[2] + dyjk * fj[2] + dzik * fi[1] + dzjk * fj[1]);

                for (m = 0; m < 6; m++) {
                  vatom[i][m] = vatom[i][m] + v[m];
                  vatom[j][m] = vatom[j][m] + v[m];
                  vatom[k][m] = vatom[k][m] + v[m];
                }
              }
            }
          }
          //     end of k loop
        }
      }
    }
    //     end of j loop
  }
}
