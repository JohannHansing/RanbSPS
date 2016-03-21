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

using namespace std;


class CRod {//TODO include CPolymers into this here, by setting random sign!
private:
    //...

public:
    int _sign;
    int axis; // rod parallel to axis 0 = x, axis 1 = y, etc.
    double coord[3]; // Coordinate of rod in 2D plane orthogonal to axis. the coord parallel to axis is always 0. (see initiaion)
    
    CRod();
    CRod(int ax, double xi, double xj);
};



#endif /* CROD_H_ */
