#include "headers/CRod.h"

using namespace std;


CRod::CRod(){
    _sign = 0;
    axis = -1;
    coord[0] = 0;
    coord[1] =  0;
    coord[2] =  0;
}

CRod::CRod(int ax, double xi, double xj){
    _sign = 1;
    axis = ax;
    int ortho[2] = {1,2};
    if (axis == 1)    ortho[0]=2, ortho[1]=0;
    else if (axis == 2) ortho[0]=0, ortho[1]=1;
    coord[axis] = 0.;
    coord[ortho[0]] =  xi;
    coord[ortho[1]] =  xj;
}

