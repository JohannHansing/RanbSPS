#include "headers/CRod.h"

using namespace std;


CRod::CRod(){
    axis = -1;
    coord = Eigen::Vector3d::Zero();
}

CRod::CRod(int ax, double xi, double xj, bool ranU, boost::mt19937 *igen, bool Pointq, double dr_q){
    int i,j;
    i = ax +1;
    if (i==3) i=0;
    j=3-(i+ax);
    axis = ax;
    coord(axis) = 0.;
    coord(i) = xi;
    coord(j) = xj;
    _igen = igen;
    if (ranU){
        signs[0]=ran_sign();
        signs[1]=ran_sign();
        signs[2]=ran_sign();
        signs[3]=ran_sign();
        checksamesign();
    }
    if (Pointq) dz0_q = atob(0.,dr_q);
}

