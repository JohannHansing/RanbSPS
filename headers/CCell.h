/*
 * CCell.h
 *
 *  Created on: Feb 19, 2016
 *      Author: jh
 */

#ifndef CCELL_H_
#define CCELL_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <time.h>
#include "boost/random.hpp"
#include "CPolymers.h"

using namespace std;



class CCell {
private:
    //boost::mt19937 *_igen;   // Maybe take random number initiation from CConf at initiation
    vector<vector<CPolymers>> _polyvec; // vector to store polymer rods in cell, one vector stores polymers that are parallel to the same axis
        
    
public:
    CCell();
    CCell(boost::mt19937 &_igen, double n_rods, vector<vector<CPolymers>> inputpolyvec){
        _poolyvec = inputpolyvec
        
    }
    

public:
    CPolymers();
    void shiftPolySign(const int , const int);
    bool samesign(const int direction, const int plane_axis, const int edgenumber);
    int get_sign(const int plane_axis, const int edgenumber);
	void printSignFile(string a_w);

};



#endif /* CCELL_H_ */
