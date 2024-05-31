/*
 *  driver.h
 *
 *  Created by Patrick on 17.04.19.
 *  Modified by Yijun Wang and Jonas Luther on 19.01.24.
 *
 */
#ifndef driver_h
#define driver_h

#include <ctime>
#include <stdio.h>
#include <vector>
#include <string>
#include "math.h"
#include "control.h"
#include "cell.h"

class Driver {
public:
    void run2D_acti(string dir);
protected:
    const double _eps = 1e-20;

    // default air properties
    double _gamma_a  = 1.4;
    double _cp_a     = 1005;
    double _cv_a     = _cp_a/_gamma_a;
    double _R_a      = _cp_a - _cv_a;
    double _mu_a     = 1.872e-5;       // dynamic viscosity
    double _lambda_a = 0.02588;        // thermal conductivity

    double u_ref, p_ref;

    double  mu      () {return _mu_a;}
    double  lambda  () {return _lambda_a;}

    ///
    //  flash calculations
    ///
    double flash_re(double rho, double p) {return p / (_gamma_a-1);}
    double  flash_c(double rho, double p) {return sqrt(_gamma_a * p / rho);}
    
    void flash(double* u) {
        double rho = u[0];
        double e   =(u[3] - 0.5 * (u[1]*u[1] + u[2]*u[2]) / rho) / rho;
        
        if (e < _eps) e = _eps;

        u[4] = (_gamma_a - 1) * rho * e;    // p
        u[5] = e / _cv_a;                   // T
        u[6] = flash_c(rho, u[4]);          // c
    }
    
    ///
    ///
    ///
    double  slope       (int limiterType, double s1, double s2);
    void flux           (double* flux,
                         int limiterType1, int limiterType2,
                         double dt_s, double dt_l, double dt_r,
                         double  hll, double  hl,  double  hr,  double  hrr,
                         double* ull, double* ul,  double* ur,  double* urr,
                         double  htl, double  hbl, double  htr, double  hbr,
                         double* utl, double* ubl, double* utr, double* ubr);
    void flux_char      (double* flux,
                         int limiterType1, int limiterType2,
                         double dt_s, double dt_l, double dt_r,
                         double  hll, double  hl,  double  hr,  double  hrr,
                         double* ull, double* ul,  double* ur,  double* urr,
                         double  htl, double  hbl, double  htr, double  hbr,
                         double* utl, double* ubl, double* utr, double* ubr);
    void flux_viscous   (double* flux,
                         int     dim,
                         double  dt,
                         double* ul, double* ur,
                         double  hl, double  hr,
                         CellVertex &v1, CellVertex &v2);
    void derivatives    (CellVertex &vertex,
                         Grid       &grid,
                         vector<Cell> &cell,
                         int i, int j);
    void TCCD           (double* ux,
                         double* ul, double* ul_old, double* ur, double* ur_old,
                         double  tl, double  tl_old, double  tr, double  tr_old,
                         double  hl, double  hr);  // time corrected central difference
    void setupCase      (Control &c,
                         vector <Cell>    &cell,
                         vector <CellMap> &cellMap,
                         bool viscous,
                         int  testCase,
                         int  nx,
                         int  ny,
                         int  &dim);
    void timestepSize   (Control &c,
                         vector <Cell>    &cell,
                         vector <CellMap> &cellMap,
                         Grid &grid,
                         bool viscous,
                         int  nx,
                         int  ny,
                         int  &dim,
                         double &dt_max,
                         double &dt_min);
};

#endif /* driver_h */
