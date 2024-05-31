/*
 *  driver.cpp
 *
 *  Created by Patrick on 17.04.19.
 *  Modified by Yijun Wang and Jonas Luther on 10.05.24.
 *
 */
#include "driver.h"
#include <cmath>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <iostream>

using namespace std;

double Driver::slope (int limiterType, double s1, double s2)
{
    double r = s1/s2;
    switch (limiterType) {
        case 1:
            //////
            // minmod
            //////
            if      (s1*s2 <= 0.0)      return 0.0;
            else if (abs(s1) < abs(s2)) return s1;
            else                        return s2;
            break;
        case 2:
            //////
            // van Leer
            //////
            return (r+abs(r))/(1+abs(r)) * s2;
            break;
        case 3:
            //////
            // superbee
            //////
            return max(max(0.0,min(1.5*r,1.0)),min(r,1.5))*s2;
            break;
        case 4:
            //////
            // von Alba 2
            //////
            return 2.0*r/(r*r+1)*s2;
            break;
        case 5:
            //////
            // Koren
            //////
            return max(0.0,min(2.0*r,min((1.0+2.0*r)/3.0,2.)))*s2;
            break;
        default:
            //////
            // 1st order
            //////
            return 0.0;
    }
    return 0.0;
}

void  Driver::flux (double* flux,
                    int limiterType1, int limiterType2,
                    double dt_s, double dt_l, double dt_r,
                    double  hll, double  hl,  double  hr,  double  hrr,
                    double* ull, double* ul,  double* ur,  double* urr,
                    double  htl, double  hbl, double  htr, double  hbr,
                    double* utl, double* ubl, double* utr, double* ubr)
{
    flux_char(flux,
              limiterType1, limiterType2,
              dt_s,dt_l,dt_r,
              hll,hl,hr,hrr,
              ull,ul,ur,urr,
              htl,hbl,htr,hbr,
              utl,ubl,utr,ubr);
    return;
}

void  Driver::flux_char(double* flux,
                        int     limiterType1, int limiterType2,
                        double dt_s, double dt_l, double dt_r,
                        double  hll, double  hl,  double  hr,  double  hrr,
                        double* ull, double* ul,  double* ur,  double* urr,
                        double  htl, double  hbl, double  htr, double  hbr,
                        double* utl, double* ubl, double* utr, double* ubr)
{
    // transform into primitive variables
    double wll[4] = {ull[0], ull[1]/ull[0], ull[2]/ull[0], ull[4]};
    double wl [4] = {ul [0], ul [1]/ul [0], ul [2]/ul [0], ul [4]};
    double wr [4] = {ur [0], ur [1]/ur [0], ur [2]/ur [0], ur [4]};
    double wrr[4] = {urr[0], urr[1]/urr[0], urr[2]/urr[0], urr[4]};
    double wtl[4] = {utl[0], utl[1]/utl[0], utl[2]/utl[0], utl[4]};
    double wbl[4] = {ubl[0], ubl[1]/ubl[0], ubl[2]/ubl[0], ubl[4]};
    double wtr[4] = {utr[0], utr[1]/utr[0], utr[2]/utr[0], utr[4]};
    double wbr[4] = {ubr[0], ubr[1]/ubr[0], ubr[2]/ubr[0], ubr[4]};

    // linearize around state
    double r = (hr*wl[0] + hl*wr[0]) / (hl + hr);
    double u = (hr*wl[1] + hl*wr[1]) / (hl + hr);
    double v = (hr*wl[2] + hl*wr[2]) / (hl + hr);
    double p = (hr*wl[3] + hl*wr[3]) / (hl + hr);
    double c = flash_c(r, p);

    // characteristic variables
    double Vll[4] = {-r*c*wll[1] + wll[3], wll[0] - 1.0/(c*c)*wll[3], wll[2], r*c*wll[1] + wll[3]};
    double Vl [4] = {-r*c*wl [1] + wl [3], wl [0] - 1.0/(c*c)*wl [3], wl [2], r*c*wl [1] + wl [3]};
    double Vr [4] = {-r*c*wr [1] + wr [3], wr [0] - 1.0/(c*c)*wr [3], wr [2], r*c*wr [1] + wr [3]};
    double Vrr[4] = {-r*c*wrr[1] + wrr[3], wrr[0] - 1.0/(c*c)*wrr[3], wrr[2], r*c*wrr[1] + wrr[3]};
    double Vtl[4] = {-r*c*wtl[1] + wtl[3], wtl[0] - 1.0/(c*c)*wtl[3], wtl[2], r*c*wtl[1] + wtl[3]};
    double Vbl[4] = {-r*c*wbl[1] + wbl[3], wbl[0] - 1.0/(c*c)*wbl[3], wbl[2], r*c*wbl[1] + wbl[3]};
    double Vtr[4] = {-r*c*wtr[1] + wtr[3], wtr[0] - 1.0/(c*c)*wtr[3], wtr[2], r*c*wtr[1] + wtr[3]};
    double Vbr[4] = {-r*c*wbr[1] + wbr[3], wbr[0] - 1.0/(c*c)*wbr[3], wbr[2], r*c*wbr[1] + wbr[3]};

    // reconstruction
    double sl1[4], sr1[4], sl2[4], sr2[4];
    for(int k=0; k<4; k++) {
        // normal direction
        sl1[k] = slope(limiterType1, (Vl [k] - Vll[k]) / (0.5*(hll+hl)), (Vr [k] - Vl [k]) / (0.5*(hr+hl )));
        sr1[k] = slope(limiterType1, (Vr [k] - Vl [k]) / (0.5*(hr +hl)), (Vrr[k] - Vr [k]) / (0.5*(hr+hrr)));
        // transverse direction
        sl2[k] = slope(limiterType2, (Vtl[k] - Vl [k]) / (0.5*(htl+hl)), (Vl [k] - Vbl[k]) / (0.5*(hl+hbl)));
        sr2[k] = slope(limiterType2, (Vtr[k] - Vr [k]) / (0.5*(htr+hr)), (Vr [k] - Vbr[k]) / (0.5*(hr+hbr)));
    }

    double VF[4];
    // characteristic wave with speed u-c
    if      (u - c > 0) VF[0] = Vl[0] + sl1[0]*( 0.5*hl - 0.5*dt_l*(u-c)) - sl2[0]*0.5*dt_l*(v-c);
    else if (u - c < 0) VF[0] = Vr[0] + sr1[0]*(-0.5*hr - 0.5*dt_r*(u-c)) - sr2[0]*0.5*dt_r*(v-c);
    else                VF[0] = 0.5 * (Vl[0] + Vr[0]);
    // characteristic wave with speed u+c
    if      (u + c > 0) VF[3] = Vl[3] + sl1[3]*( 0.5*hl - 0.5*dt_l*(u+c)) - sl2[3]*0.5*dt_l*(v+c);
    else if (u + c < 0) VF[3] = Vr[3] + sr1[3]*(-0.5*hr - 0.5*dt_r*(u+c)) - sr2[3]*0.5*dt_r*(v+c);
    else                VF[3] = 0.5 * (Vl[1] + Vr[1]);
    // characteristic wave with speed u
    if (u>0) {
        for (int k=1; k<3; k++) {
            VF[k] = Vl[k] + sl1[k]*( 0.5*hl - 0.5*dt_l*u) - sl2[k]*0.5*dt_l*v;
        }
    } else if (u<0) {
        for (int k=1; k<3; k++) {
            VF[k] = Vr[k] + sr1[k]*(-0.5*hr - 0.5*dt_r*u) - sr2[k]*0.5*dt_r*v;
        }
    } else {
        for (int k=1; k<3; k++) {
            VF[k] = 0.5 * (Vl[k] + Vr[k]);
        }
    }

    double pf  = 0.5 * (VF[0] + VF[3]);
    double uf  = (VF[3] - pf) / (r*c);
    double rf  = VF[1] + pf / (c*c);
    double vf  = VF[2];
        
    flux[0] = dt_s* rf*uf;
    flux[1] = dt_s*(rf*uf*uf + pf);
    flux[2] = dt_s* rf*uf*vf;
    flux[3] = dt_s*(uf*(pf + 0.5*rf*(uf*uf + vf*vf) + flash_re(rf,pf)));

    return;
}

void Driver::derivatives(CellVertex &vertex, Grid &grid, vector<Cell> &cell, int i, int j)
{
    // cell index
    int ibl = grid._ic(i-1, j-1);
    int ibr = grid._ic(i  , j-1);
    int itl = grid._ic(i-1, j  );
    int itr = grid._ic(i  , j  );
    
    ///
    //  transform values
    ///
    double tbl     = cell[ibl].et();
    double tbr     = cell[ibr].et();
    double ttl     = cell[itl].et();
    double ttr     = cell[itr].et();
    double tbl_old = cell[ibl].et_old();
    double tbr_old = cell[ibr].et_old();
    double ttl_old = cell[itl].et_old();
    double ttr_old = cell[itr].et_old();
    
    double hl = cell[ibl].hx();
    double hr = cell[ibr].hx();
    double hb = cell[ibl].hy();
    double ht = cell[itl].hy();
    
    ///
    //  get neighbour u
    ///
    double ubl    [7], ubr    [7], utl    [7], utr    [7];
    double ubl_old[7], ubr_old[7], utl_old[7], utr_old[7];
    grid.find_u(ubl    , i-1, j-1, 0, cell[ibl].x());
    grid.find_u(ubr    , i  , j-1, 0, cell[ibr].x());
    grid.find_u(utl    , i-1, j  , 0, cell[itl].x());
    grid.find_u(utr    , i  , j  , 0, cell[itr].x());
    grid.find_u(ubl_old, i-1, j-1, 1, cell[ibl].x());
    grid.find_u(ubr_old, i  , j-1, 1, cell[ibr].x());
    grid.find_u(utl_old, i-1, j  , 1, cell[itl].x());
    grid.find_u(utr_old, i  , j  , 1, cell[itr].x());
    
    ///
    //  compute derivatives with time-corrected central difference (TCCD)
    ///
    // x-derivatives
    double uxl[3], uxr[3];
    TCCD(uxl, ubl, ubl_old, ubr, ubr_old, tbl, tbl_old, tbr, tbr_old, hl, hr);  // betwenn bottom left and right cell
    TCCD(uxr, utl, utl_old, utr, utr_old, ttl, ttl_old, ttr, ttr_old, hl, hr);  // between top left and right cell
    
    // y-derivatives
    double uyl[3], uyr[3];
    TCCD(uyl, ubl, ubl_old, utl, utl_old, tbl, tbl_old, ttl, ttl_old, hb, ht);  // between bottom and top left cell
    TCCD(uyr, ubr, ubr_old, utr, utr_old, tbr, tbr_old, ttr, ttr_old, hb, ht);  // between bottom and top right cell
    
    ///
    //  average to get values at cell vertices
    ///
    for (int k=0; k<3; k++) {
        vertex.ux()[k] = 0.5*(uxl[k] + uxr[k]);
        vertex.uy()[k] = 0.5*(uyl[k] + uyr[k]);
    }
}

void Driver::TCCD(double* ux,
                  double* ul, double* ul_old, double* ur, double* ur_old,
                  double  tl, double  tl_old, double  tr, double  tr_old,
                  double  hl, double  hr)  // time corrected central difference
{
    double dt1=0.0, dt2=0.0;
    double ur1=0.0, vr1=0.0, Tr1=0.0;   // velocity - u
    double ur2=0.0, vr2=0.0, Tr2=0.0;   // velocity - v
    double ur3=0.0, vr3=0.0, Tr3=0.0;   // temperature
    double a1=0.0, a2=0.0, a3=0.0;

    if (tl < tr - _eps) {    // right cell is more recent
        dt1 = tr - tl;
        dt2 = tr - tl_old;
        ur1 = ur[1]/ur[0];
        vr1 = ur[2]/ur[0];
        Tr1 = ur[5];
        ur2 = ul[1]/ul[0];
        vr2 = ul[2]/ul[0];
        Tr2 = ul[5];
        ur3 = ul_old[1]/ul_old[0];
        vr3 = ul_old[2]/ul_old[0];
        Tr3 = ul_old[5];
        a1  =  2/(hl+hr);
        a2  = -2/(hl+hr) * dt2/(dt2-dt1);
        a3  =  2/(hl+hr) * dt1/(dt2-dt1);
    }
    else if (tl - _eps > tr){   // left cell more recent
        dt1 = tl - tr;
        dt2 = tl - tr_old;
        ur1 = ul[1]/ul[0];
        vr1 = ul[2]/ul[0];
        Tr1 = ul[5];
        ur2 = ur[1]/ur[0];
        vr2 = ur[2]/ur[0];
        Tr2 = ur[5];
        ur3 = ur_old[1]/ur_old[0];
        vr3 = ur_old[2]/ur_old[0];
        Tr3 = ur_old[5];
        a1  = -2/(hl+hr);
        a2  =  2/(hl+hr) * dt2/(dt2-dt1);
        a3  = -2/(hl+hr) * dt1/(dt2-dt1);
    }
    else {
        a1  = -2/(hl+hr);
        a2  =  2/(hl+hr);
        ur1 = ul[1]/ul[0];
        vr1 = ul[2]/ul[0];
        Tr1 = ul[5];
        ur2 = ur[1]/ur[0];
        vr2 = ur[2]/ur[0];
        Tr2 = ur[5];
    }

    ux[0] = a1*ur1 + a2*ur2 + a3*ur3;   // dudx dudy
    ux[1] = a1*vr1 + a2*vr2 + a3*vr3;   // dvdx dvdy
    ux[2] = a1*Tr1 + a2*Tr2 + a3*Tr3;   // dTdx dTdy
}

void Driver::flux_viscous(double* flux,
                          int     dim,
                          double  dt,
                          double* ul, double* ur,
                          double  hl, double  hr,
                          CellVertex &v1, CellVertex &v2)
{
    // viscosity and conductivity
    double mu_e     = mu();
    double lambda_e = lambda();
    
    // velocities and derivatives on cell interface
    double u = 0.5 * (ul[1]/ul[0] + ur[1]/ur[0]);
    double v = 0.5 * (ul[2]/ul[0] + ur[2]/ur[0]);
    double dudx = 0.5 * (v1.ux()[0] + v2.ux()[0]);
    double dvdx = 0.5 * (v1.ux()[1] + v2.ux()[1]);
    double dTdx = 0.5 * (v1.ux()[2] + v2.ux()[2]);
    double dudy = 0.5 * (v1.uy()[0] + v2.uy()[0]);
    double dvdy = 0.5 * (v1.uy()[1] + v2.uy()[1]);
    double dTdy = 0.5 * (v1.uy()[2] + v2.uy()[2]);
    
    ///
    //  x-flux
    ///
    if (dim==1) {
        flux[1] -= dt * mu_e* (4.0/3.0*dudx - 2.0/3.0*dvdy);
        flux[2] -= dt * mu_e* (        dudy +         dvdx);
        flux[3] -= dt *(mu_e*((4.0/3.0*dudx - 2.0/3.0*dvdy)*u + (dudy + dvdx)*v) + lambda_e*dTdx);
    }
    ///
    //  y-flux
    ///
    if (dim==2) {
        flux[1] -= dt * mu_e* (        dudy +         dvdx);
        flux[2] -= dt * mu_e* (4.0/3.0*dvdy - 2.0/3.0*dudx);
        flux[3] -= dt *(mu_e*((4.0/3.0*dvdy - 2.0/3.0*dudx)*v + (dudy + dvdx)*u) + lambda_e*dTdy);
    }
    return;
}


void Driver::setupCase(Control          &c,
                       vector <Cell>    &cell,
                       vector <CellMap> &cellMap,
                       bool viscous,
                       int  testCase,
                       int  nx,
                       int  ny,
                       int  &dim)
{
    double lx       = c.getDouble("lx");         // domain length in x-direction
    double ly       = c.getDouble("ly");         // domain length in y-direction
    bool   Xrefine  = c.getBool  ("Xrefine");    // refine in x-direction?
    bool   Yrefine  = c.getBool  ("Yrefine");    // refine in y-direction?
    double fineness = c.getDouble("fineness");
    double Ma       = c.getDouble("Ma");
    double Re       = c.getDouble("Re");
    double theta    = c.getDouble("theta");
    double rho_ref  = c.getDouble("rho_ref");    
    p_ref           = c.getDouble("p_ref");
    u_ref           = c.getDouble("u_ref");

    dim = 2;
    ///
    //  initialize fields
    ///
    double x = 0.0;
    for (int i=0; i<nx; i++) {
        double hx;
        if (Xrefine == 1) {
            hx = lx/double(nx)*(1+fineness*cos(((double(i)+0.5)/double(nx))*2.0*3.1415));
        }
        else {
            hx = lx/double(nx);
        }
        x += 0.5*hx;
        double y  = 0.0;
        for (int j=0; j<ny; j++) {
            int ic = i*ny+j;
            double hy;
            if (Yrefine == 1) {
                if (testCase==3) // boundary layer
                    hy = ly * (fineness* ((double(j)+1)*(double(j)+2) - double(j)*(double(j)+1)) / ((double(ny)+1)*double(ny)) + (1-fineness)/double(ny));
                else
                    hy = ly/double(ny)*(1+fineness*cos(((double(j)+0.5)/double(ny))*2.0*3.1415));
            }
            else {
                hy = ly/double(ny);
            }
            y += 0.5*hy;
            cell[ic].i () = ic;
            cell[ic].x () = x;
            cell[ic].y () = y;
            cell[ic].hx() = hx;
            cell[ic].hy() = hy;
            // case specific
            switch (testCase) {
                case 1:{
                    ///
                    //  shock tube
                    ///
                    dim = 1;
                    cell[ic].u()[1] = 0.0;
                    cell[ic].u()[2] = 0.0;
                    if (x<0.5) {
                        cell[ic].u()[0] = rho_ref*10.0;
                        cell[ic].u()[4] = p_ref*10.0;
                    }
                    else {
                        cell[ic].u()[0] = rho_ref;
                        cell[ic].u()[4] = p_ref;
                    }
                    cell[ic].u()[3] = flash_re(cell[ic].u()[0], cell[ic].u()[4]);
                    cell[ic].u()[3] += 0.5 * (cell[ic].u()[1]*cell[ic].u()[1] + cell[ic].u()[2]*cell[ic].u()[2]) / cell[ic].u()[0];
                    break;
                }
                case 2:{
                    ///
                    //  mixing layer
                    ///
                    _mu_a = sqrt(2.0) / (4.0*Re);
                    cell[ic].u()[0] = 1.0;
                    cell[ic].u()[1] = 0.5*sqrt(2.0)*sin(2.0*3.1416*(x-y));
                    cell[ic].u()[2] = 0.5*sqrt(2.0)*sin(2.0*3.1416*(x-y));
                    cell[ic].u()[4] = 1/(2.0*_gamma_a*Ma*Ma);
                    cell[ic].u()[3] = flash_re(cell[ic].u()[0], cell[ic].u()[4]);
                    cell[ic].u()[3] += 0.5 * (cell[ic].u()[1]*cell[ic].u()[1] + cell[ic].u()[2]*cell[ic].u()[2]) / cell[ic].u()[0];
                    break;
                }
                case 3:{
                    ///
                    //  shock boundary layer
                    ///
                    double u1 = Ma*flash_c(rho_ref, p_ref);
                    _mu_a = rho_ref*u1*ly / Re;
                    
                    cell[ic].u()[0] = rho_ref;
                    cell[ic].u()[1] = rho_ref*u1;
                    cell[ic].u()[2] = 0;
                    cell[ic].u()[4] = p_ref;
                    cell[ic].u()[3] = flash_re(cell[ic].u()[0], cell[ic].u()[4]);
                    cell[ic].u()[3] += 0.5 * (cell[ic].u()[1]*cell[ic].u()[1] + cell[ic].u()[2]*cell[ic].u()[2]) / cell[ic].u()[0];
                    break;
                }
            }
            cellMap[ic].cell = &cell[ic];
            y += 0.5*hy;

            cell[ic].et_old() = 0;
        }
        x += 0.5*hx;
    }
    return;
}

void Driver::timestepSize(Control          &contr,
                          vector <Cell>    &cell,
                          vector <CellMap> &cellMap,
                          Grid &grid,
                          bool viscous,
                          int  nx,
                          int  ny,
                          int  &dim,
                          double &dt_max,
                          double &dt_min)
{
    double cfl        = contr.getDouble("cfl");        // local cfl number
    double alpha_diff = contr.getDouble("alpha_diff"); // local diffusion number
    bool   localTS    = contr.getBool  ("localTS");
    
    for (int ic=0; ic<nx*ny; ic++) {
        if(cell[ic].u()[0] < 0.0001) cell[ic].u()[0] = 0.0001;
        //
        // for compressible Euler
        //
        double cc = cell[ic].u()[6];    // speed of sound
        cell[ic].dt() = cfl*min(cell[ic].hx()/(abs(cell[ic].u()[1])/cell[ic].u()[0] + cc),
                                cell[ic].hy()/(abs(cell[ic].u()[2])/cell[ic].u()[0] + cc));
        //
        // viscous flow
        //
        if (viscous) {
            double Dmax = max(mu(), lambda());
            
            cell[ic].dt() = min(cell[ic].dt(),
                                alpha_diff/Dmax * min(cell[ic].hx()*cell[ic].hx(),
                                                      cell[ic].hy()*cell[ic].hy()));
        }
    }
    dt_max = max_element(cellMap.begin(), cellMap.end(), CellMap::compare_dt)->dt() + 1e-15;
    dt_min = min_element(cellMap.begin(), cellMap.end(), CellMap::compare_dt)->dt();
    for (int ic=0; ic<nx*ny; ic++) {
        if (cell[ic].dt() > dt_max-1e-15) cell[ic].dt() = dt_max - 1e-15;
    }
    for (int ic=0; ic<nx*ny; ic++) {
        double dt_old = cell[ic].dt();
        if (!localTS) {
            cell[ic].l () = 0;
            cell[ic].dt() = dt_min;
            dt_max = dt_min;
        }
        else {
            if (abs(cell[ic].dt()-dt_max)<1e-15){
                cell[ic].l() = 0;
            }
            else{
                cell[ic].l() = int(ceil(log(cell[ic].dt()/dt_max)/log(0.5)));
            }
            cell[ic].dt() = dt_max * pow(0.5, cell[ic].l());
        }
        cell[ic].dt_ratio() = cell[ic].dt()/dt_old;
    }
    return;
}

void Driver::run2D_acti(string dir)
{
    //////
    //  read parameters from control file
    //////
    Control c("control.txt", dir);                              // control object
    c.read();
    int    testCase     = c.getInt   ("testCase");              // test case
    int    limiterType1 = c.getInt   ("limiterType1");          // limiter type normal direction
    int    limiterType2 = c.getInt   ("limiterType2");          // limiter type cross  direction
    int    nx           = c.getInt   ("nx");                    // number of cells in x-direction
    int    ny           = c.getInt   ("ny");                    // number of cells in y-direction
    int    plotInterval = c.getInt   ("plotInterval");          // plot interval
    int    itmax        = c.getInt   ("itmax");                 // max time steps
    bool   viscous      = c.getBool  ("viscous");               // viscous flow?
    double T            = c.getDouble("T");                     // end    time
    double Ma           = c.getDouble("Ma");
    double Re           = c.getDouble("Re");
    double T_hat        = c.getDouble("T_hat");
    
    if (testCase==2) {
        T = sqrt(2.0) / (4*3.1416*3.1416) * Re * T_hat;
        cout << "Solution time for mixing layer case: " << T << endl;
    }
    
    ///
    //  initialization
    ///
    /// 1. solution output
    ofstream result, case_info;
    
    /// 2. Grids and fields
    vector <Cell>       cell   ( nx   * ny   );
    vector <CellMap>    cellMap( nx   * ny   );
    vector <CellVertex> vertex ((nx+1)*(ny+1));
    Grid grid(c, cell.data());

    /// 3. Setup initial conditions
    int  dim;
    setupCase(c,cell,cellMap,viscous,testCase,nx,ny,dim);
    
    for (int ic=0; ic<nx*ny; ic++) {
        flash(cell[ic].u());
        for (int k=0; k<7; k++) cell[ic].uold()[k] = cell[ic].u()[k];
    }

    // write a case info file
    case_info.open(dir + "case_info.csv");
    case_info << "case,nx,ny,dim,T_hat" << endl;
    case_info << testCase << "," << nx << "," << ny << "," << dim << "," << T_hat << endl;

    //////
    // counters - initialization
    //////
    long long int counter_calc_flux = 0;
    long long int counter_update_u  = 0;
    
    double dt_max, dt_min;
    double t      = 0.0;
    int    it     = 0;
    int    it_sub = 0;
    
    std::clock_t start;
    double time;
    cout << "simulation starts" << endl;
    start = std::clock();
    
    ///
    //  time step loop
    ///
    while (t<T) {
        it++;
        double dt_factor = 1.0;
        if (t >= T || it > itmax) break;
        ///
        //  time step size for each cell
        ///
        timestepSize(c,cell,cellMap,grid,viscous,nx,ny,dim,dt_max,dt_min);
        
        if (t+dt_max> T) {
            dt_factor = (T-t)/dt_max;
        }
        double dt = dt_max;
        t += dt*dt_factor;
        
        //  initialize every global time step
        for (int ic=0; ic<nx*ny; ic++) {
            for (int k=0; k<4; k++) {
                cell[ic].fl()[k] = 0.0;
                cell[ic].fr()[k] = 0.0;
                cell[ic].gl()[k] = 0.0;
                cell[ic].gr()[k] = 0.0;
            }
            cell[ic].st() = 0.0;
            cell[ic].et() = cell[ic].dt();
        }
        
        while (min_element(cellMap.begin(),cellMap.end(),CellMap::compare_st)->st() < dt*(1.0-pow(0.5,20))) { // HACK: 1 - (1/2)^20
            ///
            //  sort cells for ACTI
            ///
            // sort by end time, small->large
            sort(cellMap.begin(), cellMap.end(), CellMap::compare_et);
            
            // find all cells with et = min(et)
            vector<CellMap>::iterator cellMap_end_et_min = cellMap.end();
            for (auto im=cellMap.begin(); im!=cellMap.end(); im++) {
                if (cell[im->i()].et() > cell[cellMap[0].i()].et()+1e-15) {cellMap_end_et_min = im; break;}
            }
            
            // sort by time step level, high->low
            sort(cellMap.begin(), cellMap_end_et_min, CellMap::compare_l);
            it_sub++;
            
            ///
            //  compute derivatives (stored on cell vertices)
            ///
            for (auto im=cellMap.begin(); im!=cellMap_end_et_min; im++) {
                int i   = im->i()/ny;
                int j   = im->i()%ny;

                // index of vertices of cell (i, j)
                int ibl =  i   *(ny+1) + j;     // bottom-left vertex
                int ibr = (i+1)*(ny+1) + j;     // bottom-right vertex
                int itl =  i   *(ny+1) + j+1;   // top-left vertex
                int itr = (i+1)*(ny+1) + j+1;   // top-right vertex

                if (vertex[ibl].update_times < it_sub) {    // to avoid repeated calculation
                    derivatives(vertex[ibl], grid, cell, i  , j  );
                    vertex[ibl].update_times = it_sub;
                }

                if (vertex[ibr].update_times < it_sub) {
                    derivatives(vertex[ibr], grid, cell, i+1, j  );
                    vertex[ibr].update_times = it_sub;
                }

                if (vertex[itl].update_times < it_sub) {
                    derivatives(vertex[itl], grid, cell, i  , j+1);
                    vertex[itl].update_times = it_sub;
                }

                if (vertex[itr].update_times < it_sub) {
                    derivatives(vertex[itr], grid, cell, i+1, j+1);
                    vertex[itr].update_times = it_sub;
                }
            }         
                        
            //////
            // Update all cells with et = min(et)
            //////
            for (auto im=cellMap.begin(); im!=cellMap_end_et_min; im++) {
                int i  = im->i()/ny;
                int j  = im->i()%ny;
                int ic = grid._ic(i  , j  );   // center
                int il = grid._ic(i-1, j  );   // left
                int ir = grid._ic(i+1, j  );   // right
                int jl = grid._ic(i  , j-1);   // bottom
                int jr = grid._ic(i  , j+1);   // top
                
                // 4 vertices of cell (i, j)
                CellVertex vbl = vertex[ i   *(ny+1) + j  ];
                CellVertex vbr = vertex[(i+1)*(ny+1) + j  ];
                CellVertex vtl = vertex[ i   *(ny+1) + j+1];
                CellVertex vtr = vertex[(i+1)*(ny+1) + j+1];
                
                // update time
                cell[ic].et_old() = cell[ic].et();
                cell[ic].st() += cell[ic].dt();
                cell[ic].et() += cell[ic].dt();
                
                //////
                // left flux
                //////
                if (i==0 || cell[il].st()<(cell[ic].st()-_eps)) {
                    // extra points in the stencil
                    int ill = grid._ic(i-2, j  );
                    int itl = grid._ic(i-1, j+1);
                    int ibl = grid._ic(i-1, j-1);
                    
                    // stencil for flux computation
                    double uc[7], ull[7], ul[7], ur[7], ubl[7], ubr[7], utl[7], utr[7];
                    grid.find_u(uc , i  , j,   0, cell[ic].x());
                    grid.find_u(ull, i-2, j,   0, cell[ill].x());
                    grid.find_u(ul , i-1, j,   0, cell[il].x());
                    grid.find_u(ur , i+1, j,   0, cell[ir].x());
                    grid.find_u(ubl, i-1, j-1, 0, cell[ibl].x());
                    grid.find_u(ubr, i  , j-1, 0, cell[jl].x());
                    grid.find_u(utl, i-1, j+1, 0, cell[itl].x());
                    grid.find_u(utr, i  , j+1, 0, cell[jr].x());
                    
                    double dt_l = 2*(cell[ic].st() - cell[il].st()) - cell[ic].dt();
                    double dt_r = cell[ic].dt();
                    
                    flux(cell[ic].fl(),
                         limiterType1, limiterType2,
                         cell[ic].dt()*dt_factor,dt_l,dt_r,
                         cell[ill].hx(),cell[il].hx(),cell[ic].hx(),cell[ir].hx(),
                         ull,ul,uc,ur,
                         cell[itl].hy(),cell[ibl].hy(),cell[jr].hy(),cell[jl].hy(),
                         utl,ubl,utr,ubr);
                    
                    // add viscous flux
                    if (viscous) {
                        flux_viscous(cell[ic].fl(),
                                     1,
                                     cell[ic].dt()*dt_factor,
                                     ul, uc,
                                     cell[il].hx(), cell[ic].hx(),
                                     vtl, vbl);
                    }
                    
                    // ensure flux conservation
                    if (i>0) for (int k=0; k<4; k++) cell[il].fr()[k] += cell[ic].fl()[k];
                    counter_calc_flux++;
                }
                
                //////
                // right flux
                //////
                if (i==nx-1 || cell[ir].st()<(cell[ic].st()-_eps)) {
                    // extra points in the stencil
                    int irr = grid._ic(i+2, j  );
                    int itr = grid._ic(i+1, j+1);
                    int ibr = grid._ic(i+1, j-1);
                    
                    // stencil for flux computation
                    double uc[7], ul[7], ur[7], urr[7], ubl[7], ubr[7], utl[7], utr[7];
                    grid.find_u(uc , i  , j,   0, cell[ic].x());
                    grid.find_u(ul , i-1, j,   0, cell[il].x());
                    grid.find_u(ur , i+1, j,   0, cell[ir].x());
                    grid.find_u(urr, i+2, j,   0, cell[irr].x());
                    grid.find_u(ubl, i  , j-1, 0, cell[jl].x());
                    grid.find_u(ubr, i+1, j-1, 0, cell[ibr].x());
                    grid.find_u(utl, i  , j+1, 0, cell[jr].x());
                    grid.find_u(utr, i+1, j+1, 0, cell[itr].x());

                    double dt_l = cell[ic].dt();
                    double dt_r = 2*(cell[ic].st() - cell[ir].st()) - cell[ic].dt();
                    
                    flux(cell[ic].fr(),
                         limiterType1, limiterType2,
                         cell[ic].dt()*dt_factor,
                         dt_l,dt_r,
                         cell[il].hx(),cell[ic].hx(),cell[ir].hx(),cell[irr].hx(),
                         ul,uc,ur,urr,
                         cell[jr].hy(),cell[jl].hy(),cell[itr].hy(),cell[ibr].hy(),
                         utl,ubl,utr,ubr);
                    
                    // add viscous flux
                    if (viscous) {
                        flux_viscous(cell[ic].fr(),
                                     1,
                                     cell[ic].dt()*dt_factor,
                                     uc, ur,
                                     cell[ic].hx(), cell[ir].hx(),
                                     vtr, vbr);
                    }
                    
                    if (i<nx-1) for (int k=0; k<4; k++) cell[ir].fl()[k] += cell[ic].fr()[k];
                    counter_calc_flux++;
                }
                
                //////
                // bottom flux
                //////
                if (dim==2 && (j==0 || cell[jl].st()<(cell[ic].st()-_eps))) {
                    // extra points in the stencil
                    int jll = grid._ic(i  , j-2);
                    int ibl = grid._ic(i-1, j-1);
                    int ibr = grid._ic(i-1, j+1);
                    
                    // stencil for flux computation
                    double uc[7], ull[7], ul[7], ur[7], ubl[7], ubr[7], utl[7], utr[7];
                    double gl[4];
                    grid.find_u(uc , i,   j  , 0, cell[ic].x());
                    grid.find_u(ull, i,   j-2, 0, cell[jll].x());
                    grid.find_u(ul , i,   j-1, 0, cell[jl].x());
                    grid.find_u(ur , i,   j+1, 0, cell[jr].x());
                    grid.find_u(ubl, i-1, j-1, 0, cell[ibl].x());
                    grid.find_u(ubr, i-1, j  , 0, cell[il].x());
                    grid.find_u(utl, i+1, j-1, 0, cell[ibr].x());
                    grid.find_u(utr, i+1, j  , 0, cell[ir].x());
                    
                    // swap u and v
                    grid.swap(ull); grid.swap(ubl);
                    grid.swap(ul ); grid.swap(ubr);
                    grid.swap(uc ); grid.swap(utl);
                    grid.swap(ur ); grid.swap(utr);
                    
                    double dt_l = 2*(cell[ic].st() - cell[jl].st()) - cell[ic].dt();
                    double dt_r = cell[ic].dt();
                    
                    flux(gl,
                         limiterType1, limiterType2,
                         cell[ic].dt()*dt_factor,
                         dt_l,dt_r,
                         cell[jll].hy(),cell[jl].hy(),cell[ic].hy(),cell[jr].hy(),
                         ull,ul,uc,ur,
                         cell[ibr].hx(),cell[ibl].hx(),cell[ir].hy(),cell[il].hx(),
                         utl,ubl,utr,ubr);
                    
                    // assign and swap fluxes
                    grid.swap_f(cell[ic].gl(), gl);
                    
                    // add viscous flux
                    if (viscous) {
                        grid.swap(ul);
                        grid.swap(uc);
                        
                        flux_viscous(cell[ic].gl(),
                                     2,
                                     cell[ic].dt()*dt_factor,
                                     ul, uc,
                                     cell[jl].hy(), cell[ic].hy(),
                                     vbl, vbr);
                    }
                    
                    if (j>0) for (int k=0; k<4; k++) cell[jl].gr()[k] += cell[ic].gl()[k];
                    counter_calc_flux++;
                }
                
                //////
                // top flux
                //////
                if (dim==2 && (j==ny-1 || cell[jr].st()<(cell[ic].st()-_eps))) {
                    // extra points in the stencil
                    int jrr = grid._ic(i  , j+2);
                    int itl = grid._ic(i-1, j+1);
                    int itr = grid._ic(i+1, j+1);
                    
                    // stencil for flux computation
                    double uc[7], ul[7], ur[7], urr[7], ubl[7], ubr[7], utl[7], utr[7];
                    double gr[4];
                    grid.find_u(uc , i,   j  , 0, cell[ic].x());
                    grid.find_u(ul , i,   j-1, 0, cell[jl].x());
                    grid.find_u(ur , i,   j+1, 0, cell[jr].x());
                    grid.find_u(urr, i,   j+2, 0, cell[jrr].x());
                    grid.find_u(ubl, i-1, j  , 0, cell[il].x());
                    grid.find_u(ubr, i-1, j+1, 0, cell[itl].x());
                    grid.find_u(utl, i+1, j  , 0, cell[ir].x());
                    grid.find_u(utr, i+1, j+1, 0, cell[itr].x());
                    
                    // swap u and v
                    grid.swap(uc ); grid.swap(ubl);
                    grid.swap(ul ); grid.swap(ubr);
                    grid.swap(ur ); grid.swap(utl);
                    grid.swap(urr); grid.swap(utr);
                    
                    double dt_l = cell[ic].dt();
                    double dt_r = 2*(cell[ic].st() - cell[jr].st()) - cell[ic].dt();
                    
                    flux(gr,
                         limiterType1, limiterType2,
                         cell[ic].dt()*dt_factor,
                         dt_l,dt_r,
                         cell[jl].hy(),cell[ic].hy(),cell[jr].hy(),cell[jrr].hy(),
                         ul,uc,ur,urr,
                         cell[ir].hx(),cell[il].hx(),cell[itr].hy(),cell[itl].hx(),
                         utl,ubl,utr,ubr);
                    
                    // assign and swap fluxes
                    grid.swap_f(cell[ic].gr(), gr);
                    
                    // add viscous flux
                    if (viscous) {
                        grid.swap(uc);
                        grid.swap(ur);
                        
                        flux_viscous(cell[ic].gr(),
                                     2,
                                     cell[ic].dt()*dt_factor,
                                     uc, ur,
                                     cell[ic].hy(), cell[jr].hy(),
                                     vtl, vtr);
                    }
                    
                    if (j<ny-1) for (int k=0; k<4; k++) cell[jr].gl()[k] += cell[ic].gr()[k];
                    counter_calc_flux++;
                }
                
                //////
                // cell updates
                //////
                counter_update_u++;
                for (int k=0; k<4; k++) {
                    cell[ic].unew()[k] = cell[ic].u()[k]
                                        + (cell[ic].fl()[k] - cell[ic].fr()[k])/cell[ic].hx()
                                        + (cell[ic].gl()[k] - cell[ic].gr()[k])/cell[ic].hy();
                }
                flash(cell[ic].unew());

                for (int k=0; k<4; k++) {
                    cell[ic].fl()[k]  = 0.0;
                    cell[ic].fr()[k]  = 0.0;
                    cell[ic].gl()[k]  = 0.0;
                    cell[ic].gr()[k]  = 0.0;
                }
            } // end of for (auto im=cellMap.begin(); im!=cellMap_end_et_min; im++)
            
            ///
            //  updates
            ///
            for (auto im=cellMap.begin(); im!=cellMap_end_et_min; im++) {
                for(int k=0; k<7; k++) {
                    cell[im->i()].uold()[k] = cell[im->i()].u()[k];
                    cell[im->i()].u()[k]    = cell[im->i()].unew()[k];
                }
            }
        } // end for while
        
        ///
        //  Intermediate solution to file
        ///
        if (it%plotInterval == 0 || t==T || it==itmax) {
            string filename = dir + "Data_" + to_string(t) + ".csv";
            result.open(filename);

            // write headers
            result <<"i,j,x,y,rho,u,v,p,T,Ma,dt/dtmax,div_r,omega,level,rhou,rhoE,hx,hy" << endl;

            for (int i=0; i<nx; i++) {
                for (int j=0; j<ny; j++) {
                    int ic = i*ny+j;
                    int il = grid._ic(i-1, j  );   // left
                    int ir = grid._ic(i+1, j  );   // right
                    int jl = grid._ic(i  , j-1);   // bottom
                    int jr = grid._ic(i  , j+1);   // top

                    double r     = cell[ic].u()[0];
                    double omega = (cell[ir].u()[2]/cell[ir].u()[0] - cell[il].u()[2]/cell[il].u()[0]) / (cell[ic].hx()+0.5*(cell[ir].hx()+cell[il].hx())) -
                                   (cell[jr].u()[1]/cell[jr].u()[0] - cell[jl].u()[1]/cell[jl].u()[0]) / (cell[ic].hy()+0.5*(cell[jr].hy()+cell[jl].hy()));
                    double div_r = sqrt(pow((cell[ir].u()[0] - cell[il].u()[0]) / (cell[ic].hx()+0.5*(cell[ir].hx()+cell[il].hx())), 2.0)+
                                        pow((cell[jr].u()[0] - cell[jl].u()[0]) / (cell[ic].hy()+0.5*(cell[jr].hy()+cell[jl].hy())), 2.0));
                    double u_abs = sqrt(pow(cell[ic].u()[1]/r, 2.0) + pow(cell[ic].u()[2]/r, 2.0));

                    result
                    << i << ","                    // col 1  - i
                    << j << ","                    // col 2  - j
                    << cell[ic].x() << ","         // col 3  - x
                    << cell[ic].y() << ","         // col 4  - y
                    << cell[ic].u()[0]      << "," // col 5  - rho
                    << cell[ic].u()[1]/r    << "," // col 6  - u
                    << cell[ic].u()[2]/r    << "," // col 7  - v
                    << cell[ic].u()[4]      << "," // col 8  - p
                    << cell[ic].u()[5]      << "," // col 9  - T
                    << u_abs/cell[ic].u()[6]<< "," // col 10 - Ma
                    << cell[ic].dt_ratio()  << "," // col 11 - dt/dtmax
                    << div_r                << "," // col 12 - density gradient
                    << omega                << "," // col 13 - vorticity
                    << cell[ic].l()         << "," // col 14 - level
                    << cell[ic].u()[1]      << "," // col 15 - x momentum
                    << cell[ic].u()[3]      << "," // col 16 - energy
                    << cell[ic].hx()        << "," // col 17 - grid size x
                    << cell[ic].hy()               // col 18 - grid size y
                    << endl;
                }
            }
            result.close();
            cout << "it = " << it << "\t it_sub = " << it_sub << "\t time = " << t << endl;
        }
    } // while (t<T)
    
    time = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
    cout << endl;
    cout << "run time = " << time << endl;
    cout << "time     = " << t << endl;
    cout << "number of flux calculations      = " << counter_calc_flux << endl;
    cout << "number of updates                = " << counter_update_u  << endl;
    return;
}
