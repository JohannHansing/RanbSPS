#include "headers/CRod.h"

using namespace std;


CRod::CRod(){
    axis = -1;
    coord = Eigen::Vector3d::Zero();
}

CRod::CRod(int ax, double xi, double xj, bool ranU, boost::mt19937 *igen, bool Pointq, double dr_q, bool mixU, double uratio_in, double Cratio_in){
    int i,j;
    i = ax +1;
    if (i==3) i=0;
    j=3-(i+ax);
    axis = ax;
    coord(axis) = 0.;
    coord(i) = xi;
    coord(j) = xj;
    _igen = igen;
    uratio=uratio_in;
    Cratio=Cratio_in;
    if (ranU){
        for (int l=0;l<4;l++){
            set_sign(signs[l]);
        }
        //signs[0]=ran_sign();
        //signs[1]=ran_sign();
        //signs[2]=ran_sign();
        //signs[3]=ran_sign();
        checksamesign();
    }
    else if (mixU){
        set_sign(signs[0]);//only this sign needs to be set. the rest is the same
        signs[1]=signs[0];
        signs[2]=signs[0];
        signs[3]=signs[0];
        checksamesign();
    }
    if (Pointq) dz0_q = atob(0.,dr_q);
}

