#include "headers/CRod.h"

using namespace std;


CRod::CRod(){
    axis = -1;
    coord = Eigen::Vector3d::Zero();
}

CRod::CRod(int ax, double xi, double xj, bool ranU, boost::mt19937 *igen, bool Pointq, double dr_q, bool mixU, double uratio, double Cratio){
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
    else if (mixU){
        // first determine if att or neut according to Cratio  : P(att) = C(att)/ [ C(rep) + C(att) ] = C(att)/C(rep) / [ 1 + C(att)/C(rep) ] = Cratio/(1+Cratio)
        // hence, if the below is true, we obtain the negative sign, since we are below P(att). Note that 0 < P(att) < 1
        if (atob(0.,1.) < Cratio/(1.+Cratio)) signs[0]=-1.*uratio;
        else signs[0] = 1.;// in this case the sign does not change, we obtain the repulsive u
        signs[1]=signs[0];
        signs[2]=signs[0];
        signs[3]=signs[0];
        checksamesign();
    }
    if (Pointq) dz0_q = atob(0.,dr_q);
}

