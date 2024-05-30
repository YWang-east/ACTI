#include "cell.h"

int Grid::_i(int i) {
    ///
    // left bound
    ///
    if (i < 0) {
        switch (left) {
            case 2:     // periodic
                i += nx;
                break;
            default:    
                i  = -i - 1;
        }
    }
    ///
    // right bound
    ///
    if (i > nx-1) {
        switch (right) {
            case 2:     // periodic
                i -= nx;
                break;
            default:    // Neumann
                i  = 2*nx - i - 1;
        }
    }
    return i;
}

int Grid::_j(int j) {
    ///
    // bottom bound
    ///
    if (j < 0) {
        switch (bottom) {
            case 2:     // periodic
                j += ny;
                break;
            default:    
                j  = -j - 1;
        }
    }
    ///
    // top bound
    ///
    if (j > ny-1) {
        switch (top) {
            case 2:     // periodic
                j -= ny;
                break;
            default:    
                j  = 2*ny - j - 1;
        }
    }
    return j;
}

void Grid::find_u(double *u, int i, int j, bool old, double x) {
    int ic = _ic(i, j);
    
    if (old) for (int k=0; k<7; k++) u[k] = cell[ic].uold()[k];
    else     for (int k=0; k<7; k++) u[k] = cell[ic].u()[k];
    
    ///
    //  left boundary
    ///
    if (i < 0) {
        if (left == 1) {   // inflow
            u[0] = rho_ref;
            u[1] = rho_ref* Ma*sqrt(1.4*p_ref/rho_ref);
            u[2] = 0.0;
            u[4] = p_ref;
            u[3] = 2.5*u[4] + 0.5/u[0]*(u[1]*u[1]+u[2]*u[2]);
            u[5] = u[4] / (0.4*1005/1.4*u[0]);
            u[6] = sqrt(1.4 * u[4]/u[0]);
        }
        else if (left == 3) {   // wall
            u[1] *= -1;
        }
    }
    ///
    //  right boundary
    ///
    if (i > nx-1) {
        if (right == 3) {   // wall
            u[1] *= -1;
        }
    }
    ///
    //  bottom boundary
    ///
    if (j < 0) {
        if    (bottom == 3) {   // wall
            u[2] *= -1;
        } else if (bottom == 4) { // no slip wall
            u[1] *= -1;
            u[2] *= -1;
        }
    }
    ///
    //  top boundary
    ///
    if (j > ny-1) {
        if (top == 1) {   // inflow
            if (testCase==3) {
                // oblique shock inlet
                double u1 = Ma*sqrt(1.4*p_ref/rho_ref);
                double r_ratio = 1.0/(5.0/6.0/(Ma*Ma*sin(theta*3.1416/180)*sin(theta*3.1416/180)) + 1.0/6.0);
                double p_ratio = 1 + 7.0/6.0*(Ma*Ma*sin(theta*3.1416/180)*sin(theta*3.1416/180)-1);
                double u2n     = u1*sin(theta*3.1416/180)/r_ratio;
                double u2t     = u1*cos(theta*3.1416/180);

                if (x>0.5*ly) {
                    u[0] = rho_ref*r_ratio;
                    u[1] = u[0]*(u2n*sin(theta*3.1416/180) + u2t*cos(theta*3.1416/180));
                    u[2] = u[0]*(u2n*cos(theta*3.1416/180) - u2t*sin(theta*3.1416/180));
                    u[4] = p_ref*p_ratio;
                } else {
                    u[0] = rho_ref;
                    u[1] = rho_ref*u1;
                    u[2] = 0.0;
                    u[4] = p_ref;
                }
                u[3] = 2.5*u[4] + 0.5/u[0]*(u[1]*u[1]+u[2]*u[2]);
                u[5] = u[4]/(0.4*1005/1.4*u[0]);
                u[6] = sqrt(1.4*u[4]/u[0]);
            }
        }
        else if (top == 3) {   // wall
            u[2] *= -1;
        }
    }
}

void Grid::swap(double *u) {
    double temp = u[2];
    u[2] = u[1];
    u[1] = temp;
}

void Grid::swap_f(double *f1, double *f0) {
    f1[0] = f0[0];
    f1[1] = f0[2];
    f1[2] = f0[1];
    f1[3] = f0[3];
}
