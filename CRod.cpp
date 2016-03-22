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
    int i,j;
    i = ax +1;
    if (i==3) i=0;
    j=3-(i+ax);
    _sign = 1;
    axis = ax;
    coord[axis] = 0.;
    coord[i] =  xi;
    coord[j] =  xj;
}

