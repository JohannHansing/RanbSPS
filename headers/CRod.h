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
#include <string>
#include "boost/random.hpp"

using namespace std;


class CRod {//TODO include CPolymers into this here, by setting random sign!
private:
    //...

public:
    int axis; // rod parallel to axis 0 = x, axis 1 = y, etc.
    double coord[3]; // Coordinate of rod in 2D plane orthogonal to axis. the coord parallel to axis is always 0. (see initiaion)
    array <int,3> signs; //this array stores the signs of the charge along the polymer backbone
    array <bool,2> samesign;
    boost::mt19937 *_igen;
    
    CRod();
    CRod(int ax, double xi, double xj, bool ranU, boost::mt19937 *igen);
    void shiftSigns(int exitmarker){
        // When the particle leaves the central cell and the system is shifted, the signs along the polymer need to be shifted, too.
        if (exitmarker==1){
            signs[0]=signs[1];
            signs[1]=signs[2];
            signs[2]=ran_sign();
        }
        if (exitmarker==-1){
            signs[2]=signs[1];
            signs[1]=signs[0];
            signs[0]=ran_sign();
        }
        checksamesign();
    }
    
    
    void checksamesign(){
        samesign[0] = (signs[1]==signs[0] );
        samesign[1] = (signs[1]==signs[2] );
    }
    
    int get_sign(){ return signs[1]; }
    
    int ran_sign(){
        boost::variate_generator<boost::mt19937&, boost::uniform_int<>> zeroone(*_igen, boost::uniform_int<>(0, 1));
	return (zeroone() * 2) - 1; //this calculation makes value either -1 or 1
    }
};



#endif /* CROD_H_ */
