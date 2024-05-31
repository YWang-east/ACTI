/*
 *  cell.h
 *
 *  Created  by Patrick Jenny on 17.04.2019
 *  Modified by Yijun Wang and Jonas Luther on 19.01.24.
 *
 */
#ifndef _Cell_h
#define _Cell_h

#include "control.h"
#include <math.h>

///
//  sits on cell center
///
class Cell {
public:
    Cell () {}
    double* u   (){return _u; }
    double* uold(){return _uold;}
    double* unew(){return _unew;}
    double* fl  (){return _fl;}
    double* fr  (){return _fr;}
    double* gl  (){return _gl;}
    double* gr  (){return _gr;}
    double& x   (){return _x; }
    double& y   (){return _y; }
    double& hx  (){return _hx;}
    double& hy  (){return _hy;}
    double& dt  (){return _dt;}
    double& dt_ratio(){return _dt_ratio;}
    double& et  (){return _et;}
    double& st  (){return _st;}
    int&    l   (){return _l; }
    int&    i   (){return _i; }
    double& et_old   (){return _et_old;}
    double& lambda   (){return _lambda;}
protected:
    double _u   [7];
    double _uold[7];   
    double _unew[7];   //for cell updates
    double _fl[4];
    double _fr[4];
    double _gl[4];
    double _gr[4];
    double _x;
    double _y;
    double _hx;
    double _hy;
    double _dt;
    double _dt_ratio;
    double _lambda;
    double _et;
    double _st;
    double _et_old;
    int    _l;
    int    _i;
};

///
//  to order cells for ACTI
///
class CellMap {
public:
    CellMap () {}
    static bool compare_st (CellMap c1, CellMap c2) {return c1.st()<c2.st();} // start time
    static bool compare_et (CellMap c1, CellMap c2) {return c1.et()<c2.et();} // end time
    static bool compare_dt (CellMap c1, CellMap c2) {return c1.dt()<c2.dt();}
    static bool compare_l  (CellMap c1, CellMap c2) {return c1.l ()>c2.l ();}
    double  st(){return cell->st();}
    double  et(){return cell->et();}
    double  dt(){return cell->dt();}
    int     l (){return cell->l ();}
    int     i (){return cell->i ();}
    Cell*   cell;
};

///
//  sits on cell vertices
///
class CellVertex {
public:
    double* ux(){return _ux;}   // dudx, dvdx, dTdx, dadx
    double* uy(){return _uy;}   // dudy, dvdy, dTdy, dady
    int     update_times;
protected:
    double _ux[3];
    double _uy[3];
};

///
//  Grid wrapper - handle boundary conditions
///
class Grid {
protected:
    Cell*  cell;
    
    int    nx, ny;
    double lx, ly;
    
    // boundary condition types
    int    testCase;
    bool   caseDefault;
    int    left;
    int    right;
    int    bottom;
    int    top;
    
    // inflow
    double Re, Ma, theta;
    double rho_ref;
    double u_ref, v_ref, p_ref;

public:
    
    Grid(Control &c, Cell *cell) : cell(cell)
    {
        nx          = c.getInt   ("nx");
        ny          = c.getInt   ("ny");
        lx          = c.getDouble("lx");
        ly          = c.getDouble("ly");
        Re          = c.getDouble("Re");
        Ma          = c.getDouble("Ma");
        theta       = c.getDouble("theta");
        rho_ref     = c.getDouble("rho_ref");
        u_ref       = c.getDouble("u_ref");
        v_ref       = c.getDouble("v_ref");
        p_ref       = c.getDouble("p_ref");
        testCase    = c.getInt   ("testCase");
        left        = c.getInt   ("left");
        right       = c.getInt   ("right");
        bottom      = c.getInt   ("bottom");
        top         = c.getInt   ("top");
    }
    int BC_left     (){return left  ;}
    int BC_right    (){return right ;}
    int BC_bottom   (){return bottom;}
    int BC_top      (){return top   ;}
    int _i          (int i);    // bounded x-index
    int _j          (int j);    // bounded y-index
    int _ic         (int i, int j){return _i(i)*ny + _j(j);}
    
    void find_u     (double* u, int i, int j, bool old, double x);
    void swap       (double* u);
    void swap_f     (double* f1, double* f0);
   
};

#endif
