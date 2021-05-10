#include "meam.h"
#include <cmath>
#include "memory.h"
#include "math_special.h"
#include <iostream>

using namespace LAMMPS_NS;
using namespace std;

void
MEAM::meam_dens_setup(int atom_nmax, int nall, int n_neigh)
{
  int i, j;

  // grow local arrays if necessary
  if (atom_nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(rho0);
    memory->destroy(rho1);
    memory->destroy(rho2);
    memory->destroy(rho3);
    memory->destroy(frhop);
    memory->destroy(frhopp);
    memory->destroy(gamma);
    memory->destroy(dgamma1);
    memory->destroy(dgamma2);
    memory->destroy(dgamma3);
    memory->destroy(arho2b);
    memory->destroy(arho1);
    memory->destroy(arho2);
    memory->destroy(arho3);
    memory->destroy(arho3b);
//     memory->destroy(darho2b);
//     memory->destroy(darho1dr);
//     memory->destroy(darho2dr);
//     memory->destroy(darho3dr);
//     memory->destroy(darho3bdr);
    memory->destroy(t_ave);
    memory->destroy(tsq_ave);
    memory->destroy(Zarray);
    memory->destroy(G_array);
    memory->destroy(dG_array);
    memory->destroy(ddG_array);
    memory->destroy(dGbar_array);
    memory->destroy(ddGbar_array);
    memory->destroy(rho_bkgd_array);

    nmax = atom_nmax;

    memory->create(rho, nmax, "pair:rho");
    memory->create(rho0, nmax, "pair:rho0");
    memory->create(rho1, nmax, "pair:rho1");
    memory->create(rho2, nmax, "pair:rho2");
    memory->create(rho3, nmax, "pair:rho3");
    memory->create(frhop, nmax, "pair:frhop");
    memory->create(frhopp, nmax, "pair:frhopp");
    memory->create(gamma, nmax, "pair:gamma");
    memory->create(dgamma1, nmax, "pair:dgamma1");
    memory->create(dgamma2, nmax, "pair:dgamma2");
    memory->create(dgamma3, nmax, "pair:dgamma3");
    memory->create(arho2b, nmax, "pair:arho2b");
    memory->create(arho1, nmax, 3, "pair:arho1");
    memory->create(arho2, nmax, 6, "pair:arho2");
    memory->create(arho3, nmax, 10, "pair:arho3");
    memory->create(arho3b, nmax, 3, "pair:arho3b");
//     memory->create(darho2b, nmax, "pair:darho2b");
//     memory->create(darho1dr, nmax, 3, "pair:darho1dr");
//     memory->create(darho2dr, nmax, 6, "pair:darho2dr");
//     memory->create(darho3dr, nmax, 10, "pair:darho3dr");
//     memory->create(darho3bdr, nmax, 3, "pair:darho3bdr");
    memory->create(t_ave, nmax, 3, "pair:t_ave");
    memory->create(tsq_ave, nmax, 3, "pair:tsq_ave");
    memory->create(Zarray, nmax, "pair:Zarray");
    memory->create(G_array, nmax, "pair:G_array");
    memory->create(dG_array, nmax, "pair:dG_array");
    memory->create(ddG_array, nmax, "pair:dG_array");
    memory->create(dGbar_array, nmax, "pair:dGbar_array");
    memory->create(ddGbar_array, nmax, "pair:ddGbar_array");
    memory->create(rho_bkgd_array, nmax, "pair:rho_bkgd_array");

  }

  if (n_neigh > maxneigh) {
    memory->destroy(scrfcn);
    memory->destroy(dscrfcn);
    memory->destroy(ddscrfcn);
    memory->destroy(fcpair);
    maxneigh = n_neigh;
    memory->create(scrfcn, maxneigh, "pair:scrfcn");
    memory->create(dscrfcn, maxneigh, "pair:dscrfcn");
    memory->create(ddscrfcn, maxneigh, "pair:ddscrfcn");
    memory->create(fcpair, maxneigh, "pair:fcpair");
  }

  // zero out local arrays

  for (i = 0; i < nall; i++) {
    rho0[i] = 0.0;
    arho2b[i] = 0.0;
    arho1[i][0] = arho1[i][1] = arho1[i][2] = 0.0;
//    drho0dr[i] = 0.0;
//    darho2bdr[i] = 0.0;
//    darho1dr[i][0] = darho1dr[i][1] = darho1dr[i][2] = 0.0;
    for (j = 0; j < 6; j++){
      arho2[i][j] = 0.0;
//      darho2dr[i][j] = 0.0;
    }
    for (j = 0; j < 10; j++){
      arho3[i][j] = 0.0;
//      darho3dr[i][j] = 0.0;
    }
    arho3b[i][0] = arho3b[i][1] = arho3b[i][2] = 0.0;
//    darho3bdr[i][0] = darho3bdr[i][1] = darho3bdr[i][2] = 0.0;
    t_ave[i][0] = t_ave[i][1] = t_ave[i][2] = 0.0;
    tsq_ave[i][0] = tsq_ave[i][1] = tsq_ave[i][2] = 0.0;
  }
}

void
MEAM::meam_dens_init(int i, int ntype, int* type, int* fmap, double** x,
                     int numneigh, int* firstneigh,
                     int numneigh_full, int* firstneigh_full, int fnoffset)
{
  //     Compute screening function and derivatives
  getscreen(i, &scrfcn[fnoffset], &dscrfcn[fnoffset], &ddscrfcn[fnoffset], &fcpair[fnoffset], x, numneigh, firstneigh,
            numneigh_full, firstneigh_full, ntype, type, fmap);

  //     Calculate intermediate density terms to be communicated
  calc_rho1(i, ntype, type, fmap, x, numneigh, firstneigh, &scrfcn[fnoffset], &fcpair[fnoffset]);
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void
MEAM::getscreen(int i, double* scrfcn, double* dscrfcn, double* ddscrfcn, double* fcpair, double** x, int numneigh,
                int* firstneigh, int numneigh_full, int* firstneigh_full, int /*ntype*/, int* type, int* fmap)
{
  int jn, j, kn, k;
  int elti, eltj, eltk;
  double xitmp, yitmp, zitmp, delxij, delyij, delzij, rij2, rij;
  double xjtmp, yjtmp, zjtmp, delxik, delyik, delzik, rik2 /*,rik*/;
  double xktmp, yktmp, zktmp, delxjk, delyjk, delzjk, rjk2 /*,rjk*/;
  double xik, xjk, sij, fcij, sfcij, dfcij, ddfcij, sikj, dfikj, ddfikj, cikj;
  double arg1, arg1_d;
  double dsij, ddsij;
  double Cmin, Cmax, delc, /*ebound,*/ a, coef1, coef2;
  double dCikj, ddCikj;
  double rnorm, fc, dfc, ddfc, drinv;

  drinv = 1.0 / this->delr_meam;
  elti = fmap[type[i]];
  if (elti < 0) return;

  xitmp = x[i][0];
  yitmp = x[i][1];
  zitmp = x[i][2];

  for (jn = 0; jn < numneigh; jn++) {
    j = firstneigh[jn];

    eltj = fmap[type[j]];
    if (eltj < 0) continue;

    //     First compute screening function itself, sij
    xjtmp = x[j][0];
    yjtmp = x[j][1];
    zjtmp = x[j][2];
    delxij = xjtmp - xitmp;
    delyij = yjtmp - yitmp;
    delzij = zjtmp - zitmp;
    rij2 = delxij * delxij + delyij * delyij + delzij * delzij;

    if (rij2 > this->cutforcesq) {
      dscrfcn[jn] = 0.0;
      scrfcn[jn] = 0.0;
      fcpair[jn] = 0.0;
      ddscrfcn[jn] = 0.0;
      continue;
    }

    const double rbound = this->ebound_meam[elti][eltj] * rij2;
    rij = sqrt(rij2);
    rnorm = (this->cutforce - rij) * drinv;
    sij = 1.0;

    //     if rjk2 > ebound*rijsq, atom k is definitely outside the ellipse
    for (kn = 0; kn < numneigh_full; kn++) {
      k = firstneigh_full[kn];
      if (k == j) continue;
      eltk = fmap[type[k]];
      if (eltk < 0) continue;

      xktmp = x[k][0];
      yktmp = x[k][1];
      zktmp = x[k][2];

      delxjk = xktmp - xjtmp;
      delyjk = yktmp - yjtmp;
      delzjk = zktmp - zjtmp;
      rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
      if (rjk2 > rbound) continue;

      delxik = xktmp - xitmp;
      delyik = yktmp - yitmp;
      delzik = zktmp - zitmp;
      rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
      if (rik2 > rbound) continue;

      xik = rik2 / rij2;
      xjk = rjk2 / rij2;
      a = 1 - (xik - xjk) * (xik - xjk);
      //     if a < 0, then ellipse equation doesn't describe this case and
      //     atom k can't possibly screen i-j
      if (a <= 0.0) continue;

      cikj = (2.0 * (xik + xjk) + a - 2.0) / a; //--- Eq. (4.11d)
      Cmax = this->Cmax_meam[elti][eltj][eltk];
      Cmin = this->Cmin_meam[elti][eltj][eltk];
      if (cikj >= Cmax) continue;
      //     note that cikj may be slightly negative (within numerical
      //     tolerance) if atoms are colinear, so don't reject that case here
      //     (other negative cikj cases were handled by the test on "a" above)
      else if (cikj <= Cmin) {
        sij = 0.0;
        break;
      } else {
        delc = Cmax - Cmin;
        cikj = (cikj - Cmin) / delc; //--- func. arg. in Eq.(4.11c)
        sikj = fcut(cikj); //--- Eq.(4.11c)
      }
      sij *= sikj; //--- \bar{s_{ij}} in Eq.(4.11a)
    }

    fc = dfcut(rnorm, dfc, ddfc);
    fcij = fc;
    dfcij = dfc * drinv; 
    ddfcij = ddfc * drinv * drinv; 

    //     Now compute derivatives
    dscrfcn[jn] = 0.0;
    ddscrfcn[jn] = 0.0;
    arg1_d = 0.0;
    sfcij = sij * fcij; //--- 4.11a
    if (!iszero(sfcij) && !isone(sfcij)) {
      for (kn = 0; kn < numneigh_full; kn++) {
        k = firstneigh_full[kn];
        if (k == j) continue;
        eltk = fmap[type[k]];
        if (eltk < 0) continue;

        delxjk = x[k][0] - xjtmp;
        delyjk = x[k][1] - yjtmp;
        delzjk = x[k][2] - zjtmp;
        rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
        if (rjk2 > rbound) continue;

        delxik = x[k][0] - xitmp;
        delyik = x[k][1] - yitmp;
        delzik = x[k][2] - zitmp;
        rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
        if (rik2 > rbound) continue;

        xik = rik2 / rij2;
        xjk = rjk2 / rij2;
        a = 1 - (xik - xjk) * (xik - xjk);
        //     if a < 0, then ellipse equation doesn't describe this case and
        //     atom k can't possibly screen i-j
        if (a <= 0.0) continue;

        cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
        Cmax = this->Cmax_meam[elti][eltj][eltk];
        Cmin = this->Cmin_meam[elti][eltj][eltk];
        if (cikj >= Cmax) {
          continue;
          //     Note that cikj may be slightly negative (within numerical
          //     tolerance) if atoms are colinear, so don't reject that case
          //     here
          //     (other negative cikj cases were handled by the test on "a"
          //     above)
          //     Note that we never have 0<cikj<Cmin here, else sij=0
          //     (rejected above)
        } else {
          delc = Cmax - Cmin;
          cikj = (cikj - Cmin) / delc; //--- func. arg. in 4.20b
          sikj = dfcut(cikj, dfikj, ddfikj ); //--- dfikj is (4.20b), sikj is (4.11c)
          coef1 = dfikj / (delc * sikj);
          dCikj = dCfunc(rij2, rik2, rjk2); //--- (4.17)/rij
          ddCikj = ddCfunc(rij, rij2, rik2, rjk2);
          dscrfcn[jn] = dscrfcn[jn] + coef1 * dCikj; //--- (4.21)/rij: sum over k
          dCikj *= rij;
          arg1_d += (1.0/delc)*( -(dfikj*dfikj*dCikj*dCikj)/delc/sikj/sikj+  
                                (ddfikj*dCikj*dCikj/sikj) + 
                                (dfikj*ddCikj/sikj)  ) ;
        }
      }
      coef1 = sfcij;
      coef2 = sij * dfcij / rij; //--- scaled by rij
      arg1 = dscrfcn[jn] * rij;
      dsij = sij * arg1;
      ddsij = dsij * arg1 + sij * arg1_d;
      dscrfcn[jn] = dscrfcn[jn] * coef1 - coef2; //--- (4.22a)/rij: units of s/r^2
//      ddscrfcn[jn] = - drinv * dfcij * dsij + fcij * ddsij - drinv * ( dsij * dfcij- sij * ddfcij * drinv );
      ddscrfcn[jn] = - drinv * dfc * dsij + fcij * ddsij - drinv * ( dsij * dfc - sij * ddfc * drinv ); //--- units of s/r^2
    }

    scrfcn[jn] = sij;
    fcpair[jn] = fcij;
  }
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void
MEAM::calc_rho1(int i, int /*ntype*/, int* type, int* fmap, double** x, int numneigh, int* firstneigh,
                double* scrfcn, double* fcpair)
{
  int jn, j, m, n, p, elti, eltj;
  int nv2, nv3;
  double xtmp, ytmp, ztmp, delij[3], rij2, rij, sij;
  double ai, aj, rhoa0j, rhoa1j, rhoa2j, rhoa3j, A1j, A2j, A3j;
  double rscalei, rscalej;
//  double drhoa0j, drhoa1j, drhoa2j, drhoa3j, A1j_d, A2j_d, A3j_d;
  // double G,Gbar,gam,shp[3+1];
  double ro0i, ro0j;
  double rhoa0i, rhoa1i, rhoa2i, rhoa3i, A1i, A2i, A3i;
//  double drhoa0i, drhoa1i, drhoa2i, drhoa3i, A1i_d, A2i_d, A3i_d;

  elti = fmap[type[i]];
  xtmp = x[i][0];
  ytmp = x[i][1];
  ztmp = x[i][2];
  for (jn = 0; jn < numneigh; jn++) {
    if (!iszero(scrfcn[jn])) {
      j = firstneigh[jn];
      sij = scrfcn[jn] * fcpair[jn];
      delij[0] = x[j][0] - xtmp;
      delij[1] = x[j][1] - ytmp;
      delij[2] = x[j][2] - ztmp;
      rij2 = delij[0] * delij[0] + delij[1] * delij[1] + delij[2] * delij[2];
      if (rij2 < this->cutforcesq) {
        eltj = fmap[type[j]];
        rij = sqrt(rij2);
        ai = rij / this->re_meam[elti][elti] - 1.0;
        aj = rij / this->re_meam[eltj][eltj] - 1.0;
        rscalei = this->re_meam[elti][elti];
        rscalej = this->re_meam[eltj][eltj];
        
        ro0i = this->rho0_meam[elti];
        ro0j = this->rho0_meam[eltj];
        rhoa0j = ro0j * MathSpecial::fm_exp(-this->beta0_meam[eltj] * aj) * sij;
//        drhoa0j = -this->beta0_meam[eltj] * rhoa0j / rscalej;
        rhoa1j = ro0j * MathSpecial::fm_exp(-this->beta1_meam[eltj] * aj) * sij;
//        drhoa1j = -this->beta1_meam[eltj] * rhoa1j / rscalej;
        rhoa2j = ro0j * MathSpecial::fm_exp(-this->beta2_meam[eltj] * aj) * sij;
//        drhoa2j = -this->beta2_meam[eltj] * rhoa2j / rscalej;
        rhoa3j = ro0j * MathSpecial::fm_exp(-this->beta3_meam[eltj] * aj) * sij;
//        drhoa3j = -this->beta3_meam[eltj] * rhoa3j / rscalej;
        rhoa0i = ro0i * MathSpecial::fm_exp(-this->beta0_meam[elti] * ai) * sij;
//        drhoa0i = -this->beta0_meam[elti] * rhoa0i / rscalei;
        rhoa1i = ro0i * MathSpecial::fm_exp(-this->beta1_meam[elti] * ai) * sij;
//        drhoa1i = -this->beta1_meam[elti] * rhoa1i / rscalei;
        rhoa2i = ro0i * MathSpecial::fm_exp(-this->beta2_meam[elti] * ai) * sij;
//        drhoa2i = -this->beta2_meam[elti] * rhoa2i / rscalei;
        rhoa3i = ro0i * MathSpecial::fm_exp(-this->beta3_meam[elti] * ai) * sij;
//        drhoa3i = -this->beta3_meam[elti] * rhoa3i / rscalei;
        if (this->ialloy == 1) {
          rhoa1j = rhoa1j * this->t1_meam[eltj];
          rhoa2j = rhoa2j * this->t2_meam[eltj];
          rhoa3j = rhoa3j * this->t3_meam[eltj];
//          drhoa1j = drhoa1j * this->t1_meam[eltj];
//          drhoa2j = drhoa2j * this->t2_meam[eltj];
//          drhoa3j = drhoa3j * this->t3_meam[eltj];
          rhoa1i = rhoa1i * this->t1_meam[elti];
          rhoa2i = rhoa2i * this->t2_meam[elti];
          rhoa3i = rhoa3i * this->t3_meam[elti];
 //         drhoa1i = drhoa1i * this->t1_meam[elti];
 //         drhoa2i = drhoa2i * this->t2_meam[elti];
 //         drhoa3i = drhoa3i * this->t3_meam[elti];
        }
        rho0[i] = rho0[i] + rhoa0j;
        rho0[j] = rho0[j] + rhoa0i;
        // For ialloy = 2, use single-element value (not average)
        if (this->ialloy != 2) {
          t_ave[i][0] = t_ave[i][0] + this->t1_meam[eltj] * rhoa0j;
          t_ave[i][1] = t_ave[i][1] + this->t2_meam[eltj] * rhoa0j;
          t_ave[i][2] = t_ave[i][2] + this->t3_meam[eltj] * rhoa0j;
          t_ave[j][0] = t_ave[j][0] + this->t1_meam[elti] * rhoa0i;
          t_ave[j][1] = t_ave[j][1] + this->t2_meam[elti] * rhoa0i;
          t_ave[j][2] = t_ave[j][2] + this->t3_meam[elti] * rhoa0i;
        }
        if (this->ialloy == 1) {
          tsq_ave[i][0] = tsq_ave[i][0] + this->t1_meam[eltj] * this->t1_meam[eltj] * rhoa0j;
          tsq_ave[i][1] = tsq_ave[i][1] + this->t2_meam[eltj] * this->t2_meam[eltj] * rhoa0j;
          tsq_ave[i][2] = tsq_ave[i][2] + this->t3_meam[eltj] * this->t3_meam[eltj] * rhoa0j;
          tsq_ave[j][0] = tsq_ave[j][0] + this->t1_meam[elti] * this->t1_meam[elti] * rhoa0i;
          tsq_ave[j][1] = tsq_ave[j][1] + this->t2_meam[elti] * this->t2_meam[elti] * rhoa0i;
          tsq_ave[j][2] = tsq_ave[j][2] + this->t3_meam[elti] * this->t3_meam[elti] * rhoa0i;
        }
        arho2b[i] = arho2b[i] + rhoa2j;
        arho2b[j] = arho2b[j] + rhoa2i;
//        darho2b[i] = darho2b[i] + drhoa2j;
//        darho2b[j] = darho2b[j] + drhoa2i;

        A1j = rhoa1j / rij;
        A2j = rhoa2j / rij2;
        A3j = rhoa3j / (rij2 * rij);
        A1i = rhoa1i / rij;
        A2i = rhoa2i / rij2;
        A3i = rhoa3i / (rij2 * rij);
        nv2 = 0;
        nv3 = 0;
        for (m = 0; m < 3; m++) {
          arho1[i][m] = arho1[i][m] + A1j * delij[m]; //--- Eq. 4.27(a)
          arho1[j][m] = arho1[j][m] - A1i * delij[m];
//          darho1dr[i][m] = darho1dr[i][m] + A1j_d * delij[m]; //--- deriv. Eq. 4.27(a) wrt rij 
//          darho1dr[j][m] = darho1dr[j][m] - A1i_d * delij[m];

          arho3b[i][m] = arho3b[i][m] + rhoa3j * delij[m] / rij; //---  Eq. 4.27(e)
          arho3b[j][m] = arho3b[j][m] - rhoa3i * delij[m] / rij;
//           darho3bdr[i][m] = darho3bdr[i][m] + ( drhoa3j - rhoa3j / rij ) * delij[m] / rij; //--- deriv. Eq. 4.27(e) wrt rij
//           darho3bdr[j][m] = darho3bdr[j][m] - ( drhoa3i - rhoa3i / rij ) * delij[m] / rij;
         for (n = m; n < 3; n++) {
            arho2[i][nv2] = arho2[i][nv2] + A2j * delij[m] * delij[n]; //--- Eq. 4.27(b)
            arho2[j][nv2] = arho2[j][nv2] + A2i * delij[m] * delij[n];
//            darho2dr[i][nv2] = darho2dr[i][nv2] + A2j_d * delij[m] * delij[n]; //--- deriv. Eq. 4.27(b) wrt rij
//            darho2dr[j][nv2] = darho2dr[j][nv2] + A2i_d * delij[m] * delij[n];
            nv2 = nv2 + 1;
            for (p = n; p < 3; p++) {
              arho3[i][nv3] = arho3[i][nv3] + A3j * delij[m] * delij[n] * delij[p];
              arho3[j][nv3] = arho3[j][nv3] - A3i * delij[m] * delij[n] * delij[p];
//               darho3dr[i][nv3] = darho3dr[i][nv3] + A3j_d * delij[m] * delij[n] * delij[p]; //--- deriv. Eq. 4.27(c) wrt rij 
//               darho3dr[j][nv3] = darho3dr[j][nv3] - A3i_d * delij[m] * delij[n] * delij[p];
              nv3 = nv3 + 1;
            }
          }
        }
      }
    }
  }
}

