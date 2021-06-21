#ifndef LMP_MEAM_H
#define LMP_MEAM_H

#include <cmath>
#include <cstring>
#include "math_const.h"

#define maxelt 5


namespace LAMMPS_NS {
class Memory;
  
// struct MyStruct {
//   double r3,ds, dds, recip;
//   double dUdrij, dUdsij, ddUddrij, ddUdrijds, ddUddsij;
//   double* dUdrijm, delij, ddUdrdrijm, ddUdrijmds, ddUdrmdrn;
// };

typedef enum { FCC, BCC, HCP, DIM, DIA, DIA3, B1, C11, L12, B2, CH4, LIN, ZIG, TRI } lattice_t;

class MEAM
{
public:
  MEAM(Memory* mem);
  ~MEAM();

private:
  Memory* memory;

  // cutforce = force cutoff
  // cutforcesq = force cutoff squared

  double cutforce, cutforcesq;

  // Ec_meam = cohesive energy
  // re_meam = nearest-neighbor distance
  // B_meam = bulk modulus
  // ielt_meam = atomic number of element
  // A_meam = adjustable parameter
  // alpha_meam = sqrt(9*Omega*B/Ec)
  // rho0_meam = density scaling parameter
  // delta_meam = heat of formation for alloys
  // beta[0-3]_meam = electron density constants
  // t[0-3]_meam = coefficients on densities in Gamma computation
  // rho_ref_meam = background density for reference structure
  // ibar_meam(i) = selection parameter for Gamma function for elt i,
  // lattce_meam(i,j) = lattce configuration for elt i or alloy (i,j)
  // neltypes = maximum number of element type defined
  // eltind = index number of pair (similar to Voigt notation; ij = ji)
  // phir = pair potential function array
  // phirar[1-6] = spline coeffs
  // attrac_meam = attraction parameter in Rose energy
  // repuls_meam = repulsion parameter in Rose energy
  // nn2_meam = 1 if second nearest neighbors are to be computed, else 0
  // zbl_meam = 1 if zbl potential for small r to be use, else 0
  // emb_lin_neg = 1 if linear embedding function for rhob to be used, else 0
  // bkgd_dyn = 1 if reference densities follows Dynamo, else 0
  // Cmin_meam, Cmax_meam = min and max values in screening cutoff
  // rc_meam = cutoff distance for meam
  // delr_meam = cutoff region for meam
  // ebound_meam = factor giving maximum boundary of sceen fcn ellipse
  // augt1 = flag for whether t1 coefficient should be augmented
  // ialloy = flag for newer alloy formulation (as in dynamo code)
  // mix_ref_t = flag to recover "old" way of computing t in reference config
  // erose_form = selection parameter for form of E_rose function
  // gsmooth_factor = factor determining length of G smoothing region
  // vind[23]D = Voight notation index maps for 2 and 3D
  // v2D,v3D = array of factors to apply for Voight notation

  // nr,dr = pair function discretization parameters
  // nrar,rdrar = spline coeff array parameters

  // theta = angle between three atoms in line, zigzag, and trimer reference structures
  // stheta_meam = sin(theta/2) in radian used in line, zigzag, and trimer reference structures
  // ctheta_meam = cos(theta/2) in radian used in line, zigzag, and trimer reference structures

  double Ec_meam[maxelt][maxelt], re_meam[maxelt][maxelt];
  double A_meam[maxelt], alpha_meam[maxelt][maxelt], rho0_meam[maxelt];
  double delta_meam[maxelt][maxelt];
  double beta0_meam[maxelt], beta1_meam[maxelt];
  double beta2_meam[maxelt], beta3_meam[maxelt];
  double t0_meam[maxelt], t1_meam[maxelt];
  double t2_meam[maxelt], t3_meam[maxelt];
  double rho_ref_meam[maxelt];
  int ibar_meam[maxelt], ielt_meam[maxelt];
  lattice_t lattce_meam[maxelt][maxelt];
  int nn2_meam[maxelt][maxelt];
  int zbl_meam[maxelt][maxelt];
  int eltind[maxelt][maxelt];
  int neltypes;

  double** phir;

  double **phirar, **phirar1, **phirar2, **phirar3, **phirar4, **phirar5, **phirar6, **phirar7, **phirar8, **phirar9;

  double attrac_meam[maxelt][maxelt], repuls_meam[maxelt][maxelt];

  double Cmin_meam[maxelt][maxelt][maxelt];
  double Cmax_meam[maxelt][maxelt][maxelt];
  double rc_meam, delr_meam, ebound_meam[maxelt][maxelt];
  int augt1, ialloy, mix_ref_t, erose_form;
  int emb_lin_neg, bkgd_dyn;
  double gsmooth_factor;

  int vind2D[3][3], vind3D[3][3][3];                  // x-y-z to Voigt-like index
  int v2D[6], v3D[10];                                // multiplicity of Voigt index (i.e. [1] -> xy+yx = 2

  int nr, nrar;
  double dr, rdrar;

public:
  int nmax;
  double *rho, *rho0, *rho1, *rho2, *rho3, *frhop, *frhopp;
  double *gamma, *dgamma1, *dgamma2, *dgamma3, *arho2b; //, *darho2b;
  double *G_array, *dG_array, *ddG_array, *dGbar_array, *ddGbar_array, *rho_bkgd_array;
  double *Zarray;
  double **arho1, **arho2, **arho3, **arho3b, **t_ave, **tsq_ave;
//  double **darho1dr, **darho2dr, **darho3dr, **darho3bdr;
    
  int maxneigh;
  double *scrfcn, *dscrfcn, *ddscrfcn, *fcpair;

  //angle for trimer, zigzag, line reference structures
  double stheta_meam[maxelt][maxelt];
  double ctheta_meam[maxelt][maxelt];

protected:
  // meam_funcs.cpp

  //-----------------------------------------------------------------------------
  // Cutoff function
  //
  static double fcut(const double xi) {
    double a;
    if (xi >= 1.0)
      return 1.0;
    else if (xi <= 0.0)
      return 0.0;
    else {
      // ( 1.d0 - (1.d0 - xi)**4 )**2, but with better codegen
      a = 1.0 - xi;
      a *= a; a *= a;
      a = 1.0 - a;
      return a * a;
    }
  }

  //-----------------------------------------------------------------------------
  // Cutoff function and derivative
  //
  static double dfcut(const double xi, double& dfc, double& ddfc) {
    double a, a3, a4, a1m4;
    if (xi >= 1.0) {
      dfc = 0.0;
      ddfc = 0.0;
      return 1.0;
    } else if (xi <= 0.0) {
      dfc = 0.0;
      ddfc = 0.0;
      return 0.0;
    } else {
      a = 1.0 - xi;
      a3 = a * a * a;
      a4 = a * a3;
      a1m4 = 1.0-a4;

      dfc = 8 * a1m4 * a3;
      ddfc = 8*(-3*a*a*a1m4+4*a3*a3);
      return a1m4*a1m4;
    }
  }

  //-----------------------------------------------------------------------------
  // Derivative of Cikj w.r.t. rij
  //     Inputs: rij,rij2,rik2,rjk2
  //
  static double dCfunc(const double rij2, const double rik2, const double rjk2) {
    double rij4, a, asq, b,denom;

    rij4 = rij2 * rij2;
    a = rik2 - rjk2;
    b = rik2 + rjk2;
    asq = a*a;
    denom = rij4 - asq;
    denom = denom * denom;
    return -4 * (-2 * rij2 * asq + rij4 * b + asq * b) / denom; //---(4.17a)/rij: wrong sign for rij2 * asq in the report
  }

  //-----------------------------------------------------------------------------
  // Derivative of Cikj w.r.t. rik and rjk
  //     Inputs: rij,rij2,rik2,rjk2
  //
  static void dCfunc2(const double rij2, const double rik2, const double rjk2,
               double& dCikj1, double& dCikj2) {
    double rij4, rik4, rjk4, a, denom;

    rij4 = rij2 * rij2;
    rik4 = rik2 * rik2;
    rjk4 = rjk2 * rjk2;
    a = rik2 - rjk2;
    denom = rij4 - a * a;
    denom = denom * denom;
    dCikj1 = 4 * rij2 * (rij4 + rik4 + 2 * rik2 * rjk2 - 3 * rjk4 - 2 * rij2 * a) / denom; //---(4.17b)/rik
    dCikj2 = 4 * rij2 * (rij4 - 3 * rik4 + 2 * rik2 * rjk2 + rjk4 + 2 * rij2 * a) / denom; //---(4.17c)/rjk
  }

  //-----------------------------------------------------------------------------
  // 2nd Derivative of Cikj w.r.t. rij
  //     Inputs: rij,rij2,rik2,rjk2
  //
  static double ddCfunc(const double rij, const double rij2, const double rik2, const double rjk2) {
    double rij4, a, asq, b,denom, ddenom, rij3, dcikj, ddcikj;

    rij4 = rij2 * rij2;
    rij3 = rij2 * rij;
    a = rik2 - rjk2;
    b = rik2 + rjk2;
    asq = a*a;
    denom = rij4 - asq;
    denom = denom * denom;
    ddenom = 2*(rij4-asq)*(4*rij3);
    
    dcikj = rij * dCfunc(rij2, rik2, rjk2);
    ddcikj = -4 * (-2 * rij2 * asq + rij4 * b + asq * b) - 4 * rij * (-2 * 2 * rij * asq + 4*rij3 * b ) - dcikj * ddenom;
    ddcikj /= denom;
    
    return ddcikj;    
    
  }
  //-----------------------------------------------------------------------------
  // 2nd Derivative of Cikj w.r.t. rik and rjk
  //     Inputs: rij,rij2,rik2,rjk2
  //
  static void ddCfunc2(const double rik, const double rjk, const double rij2, const double rik2, const double rjk2,
               double& ddCikj1, double& ddCikj2) {
    double rij4, rik4, rjk4, a, denom, ddenom_ik, ddenom_jk, rik3, rjk3, dCikj1, dCikj2;

    rij4 = rij2 * rij2;
    rik4 = rik2 * rik2;
    rjk4 = rjk2 * rjk2;
    rik3 = rik2 * rik;
    rjk3 = rjk2 * rjk;
    a = rik2 - rjk2;
    denom = rij4 - a * a;
    
    ddenom_ik = - 2*a * 2*rik;
    ddenom_ik *= (2 * denom);

    ddenom_jk = 2*a * 2*rjk;
    ddenom_jk *= (2 * denom);

    denom = denom * denom;

    dCfunc2(rij2, rik2, rjk2, dCikj1, dCikj2);

    dCikj1 *= rik; 
    ddCikj1 = 4 * rij2 * (rij4 + rik4 + 2 * rik2 * rjk2 - 3 * rjk4 - 2 * rij2 * a) +
              4 * rij2 * rik * ( 3 * rik3 + 4 * rik * rjk2 - 4 * rij2 * rik ) -
              ddenom_ik * dCikj1; // d(dCikj1)/drik
    ddCikj1 /= denom;
    
    dCikj2 *= rjk; 
    ddCikj2 = 4 * rij2 * (rij4 - 3 * rik4 + 2 * rik2 * rjk2 + rjk4 + 2 * rij2 * a) +
              4 * rij2 * rjk * ( 4 * rik2 * rjk + 4 * rjk3 - 4 * rij2 *rjk ) - 
              ddenom_jk * dCikj2; // d(dCikj2)/drjk
    ddCikj2 /= denom;
  }
  
  double G_gam(const double gamma, const int ibar, int &errorflag) const;
  double dG_gam(const double gamma, const int ibar, double &dG, double &ddG) const;
  double Get_ddrho1drdr(int i, 
                    double rij, double sij, 
                    double rhoa1j, double drhoa1j, double ddrhoa1j,
                    double arg1i1,
                    double arg1i1_d
                    );
  void Get_ddrho1drmdr(int i,
                       double rij, double sij, double* delij,
                       double rhoa1j, double drhoa1j,
                       double** arho1, double* darho1dri,
                       double* ddrho1drmdr1);
  void Get_ddrho1drmdrn(int i, //--- deriv of Eq. (4.30c) wrt rn
                       double rij, double sij,
                       double rhoa1j,
                       double* ddrho1drmdrn1);
  double Get_ddrho2drdr(int i, 
                    double rij, double sij, 
                    double rhoa2j, double drhoa2j, double ddrhoa2j,
                    double* arho2b,
                    double arg1i2,
                    double arg1i2_d);
  void Get_ddrho2drmdr(int i,
                       double rij, double sij, double* delij,
                       double rhoa2j, double drhoa2j,
                       double** arho2,
                       double* darho2dri,
                       double* ddrho2drmdr1);
  void Get_ddrho2drmdrn(int i,
                       double rij, double sij, double* delij,
                       double rhoa2j,
                       double** arho2,
                       double* ddrho2drmdrn1);
  double Get_ddrho3drdr( double rij, double sij, 
                    double rhoa3j, double drhoa3j, double ddrhoa3j,
                    double arg1i3, double arg3i3,
                    double arg1i3_d, double arg3i3_d
                    );
  void Get_ddrho3drmdr( int i,
                       double rij, double sij, double* delij,
                       double rhoa3j, 
                       double drhoa3j, 
                       double* darho3dri, double* darho3bdri,
                       double* ddrho3drmdr1 //--- modify 
                     );
  void Get_ddrho3drmdrn( int i,
                       double rij, double sij, double* delij,
                       double rhoa3j, 
                       double* ddrho3drmdrn1);
  double Get_ddrhodrdr(  int i, int elti,
                      double* shpi,
                      double t1i, double t2i, double t3i,
                      double dt1dr1, double dt2dr1, double dt3dr1,
                      double ddt1drdr1, double ddt2drdr1, double ddt3drdr1,
                      double* rho0, double* rho1, double* rho2, double* rho3, 
                      double drho0dr1, double drho1dr1, double drho2dr1, double drho3dr1, 
                      double ddrho0drdr1, double ddrho1drdr1, double ddrho2drdr1, double ddrho3drdr1,
                      double drhodr1 );
  void Get_ddrhodrmdr( int i, int elti, //--- deriv. of Eq. 4.36(c) wrt. r
                        double* shpi, 
                        double t1i,  double t2i,  double t3i,
                        double dt1dr1,  double dt2dr1,  double dt3dr1,
                        double* rho0, double* rho1, double* rho2, double* rho3,
                        double drho0dr1,  double drho1dr1,  double drho2dr1,  double drho3dr1, 
                        double* drho0drm1,  double* drho1drm1,  double* drho2drm1,  double* drho3drm1, 
                        double* ddrho0drmdr1, double* ddrho1drmdr1,  double* ddrho2drmdr1,  double* ddrho3drmdr1,
                        double* drhodrm1,
                        double* ddrhodrmdr1);
  void Get_ddrhodrmdrn( int i, int elti, //--- deriv. of Eq. 4.36(c) wrt. rm
                        double* shpi, 
                        double t1i,  double t2i,  double t3i,
                        double* rho0, double* rho1, double* rho2, double* rho3,
                        double* drho0drm1,  double* drho1drm1,  double* drho2drm1,  double* drho3drm1, 
                                            double* ddrho1drmdrn1,  double* ddrho2drmdrn1,  double* ddrho3drmdrn1,
                        double* ddrhodrmdrn1);
  double Get_ddrhodrds(  int i, int elti,
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
                     );
  void Get_ddrhodrmds( int i, int elti, //--- deriv. of Eq. 4.36(c) wrt. r
                        double* shpi, 
                        double t1i,  double t2i,  double t3i,
                        double dt1ds,  double dt2ds,  double dt3ds,
                        double* rho0, double* rho1, double* rho2, double* rho3,
                        double drho0ds,  double drho1ds,  double drho2ds,  double drho3ds, 
                        double* drho0drm,  double* drho1drm,  double* drho2drm,  double* drho3drm, 
                        double* ddrho0drmds, double* ddrho1drmds,  double* ddrho2drmds,  double* ddrho3drmds,
                        double* drhodrm,
                        double* ddrhodrmds //--- modify
                       );
  double GetModulus(int alpha, int beta, int gamma, int lambda,  double r3,double ds, double dds, double recip,
                        double dUdrij, double dUdsij, double ddUddrij, double ddUdrijds, double ddUddsij,
                        double* dUdrijm, double* delij, double* ddUdrdrijm, double* ddUdrijmds, double* ddUdrmdrn);
  
  static double zbl(const double r, const int z1, const int z2);
  double embedding(const double A, const double Ec, const double rhobar, double& dF, double& ddF ) const;
  static double erose(const double r, const double re, const double alpha, const double Ec, const double repuls, const double attrac, const int form);

  static void get_shpfcn(const lattice_t latt, const double sthe, const double cthe, double (&s)[3]);

  static int get_Zij2(const lattice_t latt, const double cmin, const double cmax,
                      const double sthe, double &a, double &S);
  static int get_Zij2_b2nn(const lattice_t latt, const double cmin, const double cmax, double &S);

protected:
  void meam_checkindex(int, int, int, int*, int*);
  void getscreen(int i, double* scrfcn, double* dscrfcn, double* ddscrfcn, double* fcpair, double** x, int numneigh,
                 int* firstneigh, int numneigh_full, int* firstneigh_full, int ntype, int* type, int* fmap);
  void calc_rho1(int i, int ntype, int* type, int* fmap, double** x, int numneigh, int* firstneigh,
                 double* scrfcn, double* fcpair);

  void alloyparams();
  void compute_pair_meam();
  double phi_meam(double, int, int);
  double phi_meam_series(const double scrn, const int Z1, const int Z2, const int a, const int b, const double r, const double arat);
  void compute_reference_density();
  void get_tavref(double*, double*, double*, double*, double*, double*, double, double, double, double,
                  double, double, double, int, int, lattice_t);
  void get_sijk(double, int, int, int, double*);
  void get_densref(double, int, int, double*, double*, double*, double*, double*, double*, double*, double*);
  void interpolate_meam(int);

public:
  //-----------------------------------------------------------------------------
  // convert lattice spec to lattice_t
  // only use single-element lattices if single=true
  // return false on failure
  // return true and set lat on success
  static bool str_to_lat(const char* str, bool single, lattice_t& lat)
  {
    if (strcmp(str,"fcc") == 0) lat = FCC;
    else if (strcmp(str,"bcc") == 0) lat = BCC;
    else if (strcmp(str,"hcp") == 0) lat = HCP;
    else if (strcmp(str,"dim") == 0) lat = DIM;
    else if (strcmp(str,"dia") == 0) lat = DIA;
    else if (strcmp(str,"dia3") == 0) lat = DIA3;
    else if (strcmp(str,"lin") == 0) lat = LIN;
    else if (strcmp(str,"zig") == 0) lat = ZIG;
    else if (strcmp(str,"tri") == 0) lat = TRI;
    else {
      if (single)
        return false;

      if (strcmp(str,"b1")  == 0) lat = B1;
      else if (strcmp(str,"c11") == 0) lat = C11;
      else if (strcmp(str,"l12") == 0) lat = L12;
      else if (strcmp(str,"b2")  == 0) lat = B2;
      else if (strcmp(str,"ch4")  == 0) lat = CH4;
      else if (strcmp(str,"lin")  == 0) lat =LIN;
      else if (strcmp(str,"zig")  == 0) lat = ZIG;
      else if (strcmp(str,"tri")  == 0) lat = TRI;
      else return false;
    }
    return true;
  }

  static int get_Zij(const lattice_t latt);
  void meam_setup_global(int nelt, lattice_t* lat, int* ielement, double* atwt, double* alpha,
                         double* b0, double* b1, double* b2, double* b3, double* alat, double* esub,
                         double* asub, double* t0, double* t1, double* t2, double* t3, double* rozero,
                         int* ibar);
  void meam_setup_param(int which, double value, int nindex, int* index /*index(3)*/, int* errorflag);
  void meam_setup_done(double* cutmax);
  void meam_dens_setup(int atom_nmax, int nall, int n_neigh);
  void meam_dens_init(int i, int ntype, int* type, int* fmap, double** x, int numneigh, int* firstneigh,
                      int numneigh_full, int* firstneigh_full, int fnoffset);
  void meam_dens_final(int nlocal, int eflag_either, int eflag_global, int eflag_atom, double* eng_vdwl,
                       double* eatom, int ntype, int* type, int* fmap, double** scale, int& errorflag);
  void meam_force(int i, int eflag_either, int eflag_global, int eflag_atom, int vflag_atom, double* eng_vdwl,
                  double* eatom, int ntype, int* type, int* fmap, double** scale, double** x, int numneigh, int* firstneigh,
                  int numneigh_full, int* firstneigh_full, int fnoffset, double** f, double** vatom);
};

// Functions we need for compat

static inline bool iszero(const double f) {
  return fabs(f) < 1e-20;
}

static inline bool isone(const double f) {
  return fabs(f-1.0) < 1e-20;
}

// Helper functions

static inline double fdiv_zero(const double n, const double d) {
  if (iszero(d))
    return 0.0;
  return n / d;
}

}
#endif
