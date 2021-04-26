#include "meam.h"
#include <iostream>

using namespace std
using namespace LAMMPS_NS;

void
MEAM::meam_dens_final(int nlocal, int eflag_either, int eflag_global, int eflag_atom, double* eng_vdwl,  //---- meam_dens_final called before force????
                      double* eatom, int /*ntype*/, int* type, int* fmap, double** scale, int& errorflag)
{
  int i, elti;
  int m;
  double rhob, G, dG, ddG, Gbar, dGbar, ddGbar, gam, shp[3], Z;
  double denom, rho_bkgd, Fl;
  double scaleii;

  //     Complete the calculation of density

  for (i = 0; i < nlocal; i++) {
    elti = fmap[type[i]];
    if (elti >= 0) {
      scaleii = scale[type[i]][type[i]];
      rho1[i] = 0.0;
      rho2[i] = -1.0 / 3.0 * arho2b[i] * arho2b[i];
      rho3[i] = 0.0;
      for (m = 0; m < 3; m++) {
        rho1[i] = rho1[i] + arho1[i][m] * arho1[i][m]; //--- Eq. (4.28 a)
        rho3[i] = rho3[i] - 3.0 / 5.0 * arho3b[i][m] * arho3b[i][m];
      }
      for (m = 0; m < 6; m++) {
        rho2[i] = rho2[i] + this->v2D[m] * arho2[i][m] * arho2[i][m];
      }
      for (m = 0; m < 10; m++) {
        rho3[i] = rho3[i] + this->v3D[m] * arho3[i][m] * arho3[i][m];
      }

      if (rho0[i] > 0.0) {
        if (this->ialloy == 1) {
          t_ave[i][0] = fdiv_zero(t_ave[i][0], tsq_ave[i][0]);
          t_ave[i][1] = fdiv_zero(t_ave[i][1], tsq_ave[i][1]);
          t_ave[i][2] = fdiv_zero(t_ave[i][2], tsq_ave[i][2]);
        } else if (this->ialloy == 2) {
          t_ave[i][0] = this->t1_meam[elti];
          t_ave[i][1] = this->t2_meam[elti];
          t_ave[i][2] = this->t3_meam[elti];
        } else {
          t_ave[i][0] = t_ave[i][0] / rho0[i];
          t_ave[i][1] = t_ave[i][1] / rho0[i];
          t_ave[i][2] = t_ave[i][2] / rho0[i];
        }
      }

      gamma[i] = t_ave[i][0] * rho1[i] + t_ave[i][1] * rho2[i] + t_ave[i][2] * rho3[i]; //--- \Gamma: Eq. (4.4) rho1 is (\bar{\rho}_i^{(k)})^2

      if (rho0[i] > 0.0) {
        gamma[i] = gamma[i] / (rho0[i] * rho0[i]); //--- \Gamma: Eq. (4.4) normalize by (\bar{\rho}_i^{(0)})^2
      }

      Z = get_Zij(this->lattce_meam[elti][elti]); //--- Z_{i0}
      if (not Z > 0)
        cout << "Z("<<i<<")="<<Z<<"\n";
      Zarray[ i ] = Z; //--- defined?
      
      G = G_gam(gamma[i], this->ibar_meam[elti], errorflag); //--- G(\Gamma)
      G_array[i] = G;
      
      if (errorflag != 0)
        return;

      get_shpfcn(this->lattce_meam[elti][elti], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shp);  //--- s_i^{(k)} in Eq. (4.6)

      if (this->ibar_meam[elti] <= 0) {
        Gbar = 1.0;
        dGbar = 0.0;
        ddGbar = 0.0;
        dGbar_array[ i ] = 0.0;
        ddGbar_array[ i ] = 0.0;
      } else {
        if (this->mix_ref_t == 1) {
          gam = (t_ave[i][0] * shp[0] + t_ave[i][1] * shp[1] + t_ave[i][2] * shp[2]) / (Z * Z); //--- \Gamma^{ref}: Eq. (4.6)
        } else {
          gam = (this->t1_meam[elti] * shp[0] + this->t2_meam[elti] * shp[1] + this->t3_meam[elti] * shp[2]) /
                (Z * Z);
        }
        Gbar = G_gam(gam, this->ibar_meam[elti], errorflag); //--- G(\Gamma^{ref})
      }
      rho[i] = rho0[i] * G; //--- Eq. (4.3): missing / rho_bkgd ??

      if (this->mix_ref_t == 1) {
        if (this->ibar_meam[elti] <= 0) {
          Gbar = 1.0;
          dGbar = 0.0;
          ddGbar = 0.0;
          dGbar_array[ i ] = 0.0;
          ddGbar_array[ i ] = 0.0;
        } else {
          gam = (t_ave[i][0] * shp[0] + t_ave[i][1] * shp[1] + t_ave[i][2] * shp[2]) / (Z * Z);
          Gbar = dG_gam(gam, this->ibar_meam[elti], dGbar, ddGbar );
          dGbar_array[ i ] = dGbar;
          ddGbar_array[ i ] = ddGbar;
        }
        rho_bkgd = this->rho0_meam[elti] * Z * Gbar; //--- Eq. (4.5)
      } else {
        if (this->bkgd_dyn == 1) {
          rho_bkgd = this->rho0_meam[elti] * Z;
        } else {
          rho_bkgd = this->rho_ref_meam[elti];
        }
      }
      rho_bkgd_array[i]=rho_bkgd; // rho_bkgd_array[i] defined??
      
      rhob = rho[i] / rho_bkgd;
      denom = 1.0 / rho_bkgd;
      
       
      G = dG_gam(gamma[i], this->ibar_meam[elti], dG, ddG); //--- dG/dgamma
      G_array[i] = G; //defined?
      dG_array[i] = dG; ///defined?
      ddG_array[i] = ddG; ///defined?
      
      
      dgamma1[i] = (G - 2 * dG * gamma[i]) * denom; //--- Eq. (4.36a): prefactor in the 1st term of the RHS 

      if (!iszero(rho0[i])) {
        dgamma2[i] = (dG / rho0[i]) * denom; //--- Eq. (4.36a): prefactor in the 2nd term of the RHS 
      } else {
        dgamma2[i] = 0.0;
      }

      //     dgamma3 is nonzero only if we are using the "mixed" rule for
      //     computing t in the reference system (which is not correct, but
      //     included for backward compatibility
      if (this->mix_ref_t == 1) {
        dgamma3[i] = rho0[i] * G * dGbar / (Gbar * Z * Z) * denom;
      } else {
        dgamma3[i] = 0.0;
      }

      Fl = embedding(this->A_meam[elti], this->Ec_meam[elti][elti], rhob, frhop[i], frhopp[i]);

      if (eflag_either != 0) {
        Fl *= scaleii;
        if (eflag_global != 0) {
          *eng_vdwl = *eng_vdwl + Fl;
        }
        if (eflag_atom != 0) {
          eatom[i] = eatom[i] + Fl;
        }
      }
    }
  }
}

