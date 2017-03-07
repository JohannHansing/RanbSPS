/*
 * CPolymers.h
 *
 *  Created on: Aug 9, 2013
 *      Author: jh
 */

#ifndef CROD_H_
#define CROD_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <array>
#include <string>
#include <Eigen/Dense>
#include "boost/random.hpp"


using namespace std;


class CRod {//TODO include CPolymers into this here, by setting random sign!
private:
    //...

public:
    int axis; // rod parallel to axis 0 = x, axis 1 = y, etc.
    Eigen::Vector3d coord; // Coordinate of rod in 2D plane orthogonal to axis. the coord parallel to axis is always 0. (see initiaion)
    array <double,4> signs; //this array stores the signs of the charge along the polymer backbone (4 signs, for pbc similations)
    array <bool,2> samesign;
    boost::mt19937 *_igen;
    
    // ratio stuff
    double uratio;
    double Cratio;
    
    // offset for "lowest" point charge
    double dz0_q = 0;
    
    CRod();
    CRod(int ax, double xi, double xj, bool ranU, boost::mt19937 *igen, bool Pointq=false, double dr_q=0., bool mixU=false, double uratio_in=1., double Cratio_in=1);
    void shiftSigns(int exitmarker){
        // When the particle leaves the central cell and the system is shifted, the signs along the polymer need to be shifted, too.
        if (exitmarker==1){
            signs[0]=signs[1];
            signs[1]=signs[2];
            //signs[2]=ran_sign();
            set_sign(signs[2]);
        }
        else{
            signs[2]=signs[1];
            signs[1]=signs[0];
            //signs[0]=ran_sign();
            set_sign(signs[0]);
        }
        checksamesign();
    }
    
    
    void checksamesign(){
        // NOTE: Samesign check is not implemented for pbc!
        samesign[0] = ( signs[1]==signs[0] );
        samesign[1] = ( signs[1]==signs[2] );
    }
    
    int get_sign(){ return signs[1]; }
    
    void shiftqs(double deltaq, double L, int exitmarker){
        // delta needs to be defined from the start. For alexa488 in dextran(-) it is 2/1.5 = 1.33. L is the length of the rod, i.e. 3 * _boxsize
        if (exitmarker==1){
            dz0_q = fmod( L - dz0_q , deltaq ); //fmod, since % only works for integers
        }
        else{
            dz0_q = fmod( L - (deltaq - dz0_q) , deltaq );
        }
    }
    
    int ran_sign(){
        boost::variate_generator<boost::mt19937&, boost::uniform_int<>> zeroone(*_igen, boost::uniform_int<>(0, 1));
	    return (zeroone() * 2) - 1; //this calculation makes value either -1 or 1
    }
    
    double atob(double a, double b){
        boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<>> ran_gen(*_igen, boost::random::uniform_real_distribution<>(a, b));
        return ran_gen();
    }
    
    void set_sign(double & sign){
        // first determine if att or neut according to Cratio  : P(att) = C(att)/ [ C(rep) + C(att) ] = C(att)/C(rep) / [ 1 + C(att)/C(rep) ] = Cratio/(1+Cratio)
        // hence, if the below is true, we obtain the negative sign, since we are below P(att). Note that 0 < P(att) < 1
        if (atob(0.,1.) < Cratio/(1.+Cratio)) sign=-1.*uratio;
        else sign = 1.;// in this case the sign does not change, we obtain the repulsive u
    }
};



#endif /* CROD_H_ */
