#include "headers/CRod.h"

using namespace std;


CRod::CRod(){
    axis = -1;
    coord = Eigen::Vector3d::Zero();
}

CRod::CRod(int ax, double xi, double xj, bool ranU, boost::mt19937 *igen, bool Pointq, double dr_q, bool mixU, double uratio_in, double Cratio_in, int N_patches){
    //first resize vectors appropriately
    this->N_patches = N_patches;
    signs.resize(N_patches);
    samesign.resize(N_patches-1);
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
        for (int l=0;l<N_patches;l++){
            set_sign(signs[l]);
        }
        checksamesign();
    }
    else if (mixU){
        set_sign(signs[0]);//only this sign needs to be set. the rest is the same
        for (int i=1; i<N_patches; i++){
            signs[i]=signs[0];
        }
        checksamesign();
    }
    if (Pointq) dz0_q = atob(0.,dr_q);
}

