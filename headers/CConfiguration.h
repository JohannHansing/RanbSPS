#ifndef CCONFIGURATION_H_
#define CCONFIGURATION_H_


#include <array>
#include <string>
#include <vector>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include "CPolymers.h"

#define ifdebug(x) 

using namespace std;

class CConfiguration {
    /*Class where all the configuration variables such as potRange etc. and also most functions for the
     * simulation are stored
     */
private:
    //SCALING
    double _timestep;         //This is the RESCALED timestep! timestep = dt * kT / (frictionCoeffcient * particlesize)
    double _mu_sto;
    double _pradius;     //particle size is most relevant for scaling! (default = 1)
    double _boxsize[3];          // ALWAYS define boxsize through particlesize due to scaling!
    std::array<std::array<double,3>,3> _b_array;
    double _epsilon;

    //EXPONENTIAL Potential
    double _potRange;         // Avoid values like 10/3 = 3.33333... ALWAYS define the Range of the exp potential through boxsize due to scaling!
    double _potStrength;      // rescaled U_0 for exponential Potential
    double _rodDistance;
    CPolymers _poly;
    double _cutoffExpSq;
    double _r_cSq;    //cutoff for Lennard-Jones calculation (at minimum)


    //bool Parameters
    bool _potMod;             // true if the modified exponential potential version is used that is not 3*U_0 at the intersections.
    bool _LJPot;              // if true, then LJ interaction is calculated for steric hindrance
    bool _ranU;
    bool _hpi;
    bool _ranRod;

    //COUNTERS AND INIT VALUES
    double _boxCoord[3];
    unsigned int _wallcrossings[3]; //counts how many times the particle has crossed the [entry, opposite, side] wall, while
                                            //travelling through the lattice
    int _entryside[3];            //records through which side and in which direction the particle last entered a box. For the face of the cube in the x-z
                                            //plane at y=0, it is entryside[0] = 1, in the x-y plane at z=L it is entryside[1] = -1!
    double _resetpos;
    double _startpos[3];          //Stores where the particle starting position was. This is needed to calculate the mean square displacement
    double _prevpos[3];           //Stores previous particle position before particle is moved.

    int _min, _max;        // parameters for determining up to which order neighboring rods are considered for the potential

    //Particle parameters
    double _ppos[3];    //initialize particle position (DO IT LIKE resetpos FOR MOVEPARTICLEFOLLOW/RESET)
    double _upot;
    double _f_mob[3];   //store mobility and stochastic force
    double _f_sto[3];

    // Gamma distribution
    double _alpha = 10.;
    double _beta = 1.;

    boost::mt19937 *m_igen;                      //generate instance of random number generator "twister".
    double zerotoone(){
        boost::uniform_01<boost::mt19937&> dist(*m_igen);
        return dist();
    }
    
    double atob(double a, double b){
        return a + (b-a) * zerotoone()();
    }
    
    int ran_sign(){
    // the variate generator uses _igen (int rand number generator),
    // samples from uniform integer distribution 0, 1
        boost::variate_generator<boost::mt19937&, boost::uniform_int<>> zeroone(*m_igen, boost::uniform_int<>(0, 1));
	return (zeroone() * 2) - 1; //this calculation makes value either -1 or 1
    }
    
    double ran_gamma(){
        // http://www.boost.org/doc/libs/1_58_0/doc/html/boost/random/gamma_distribution.html
        // Mean of gamma dist is alpha* beta
        boost::variate_generator<boost::mt19937&    , boost::gamma_distribution<double> > ran_gen(
                *m_igen, boost::gamma_distribution<double>(_alpha, _beta));
        return ran_gen();
    }
    
    //TODO del
    void testgamma(){

        ofstream trajectoryfile;
        trajectoryfile.open(("tmp_gamma.txt"));
        for (int i=0;i<100000;i++){
            trajectoryfile << fixed << ran_gamma() << endl;
        }
        trajectoryfile.close();
    }

    void initRanb(){
        // For now, I just use a gamma distribution
        for (int i=0;i<3;i++){
            if (_pradius < 4.5){
                _boxsize[i] = 10;
            }
            else _boxsize[i] = 2*_pradius + 1;
            _b_array[i][1] = _boxsize[i];
            _b_array[i][0] = ran_gamma();
            _b_array[i][2] = ran_gamma();
        }
    }
    
    void updateRanb(int axis, int exitmarker){
        // ROTATE http://en.cppreference.com/w/cpp/algorithm/rotate
        //copy neighbor boxsize to _boxsize for particle box.
        ifdebug(cout << "Update Ranb\naxis " << axis << "  --  exitmarker " << exitmarker << endl;)
        double newb = ran_gamma();
        ifdebug(cout << "*" << newb << endl;)
        if (exitmarker==1){
            // rotation to the left
            std::rotate(_b_array[axis].begin(), _b_array[axis].begin() + 1, _b_array[axis].end());
            // assign new boxsize on right side
            _b_array[axis].back() = newb;
        }
        else {
            // rotation to the  right
            std::rotate(_b_array[axis].rbegin(), _b_array[axis].rbegin() + 1, _b_array[axis].rend());
            // assign new boxsize on right side
            //ifdebug(cout << "B4 _b_array[axis].front() = " << _b_array[axis].front() << endl;)
            _b_array[axis].front() = newb;
            //ifdebug(cout << "AFTER _b_array[axis].front() = " << _b_array[axis].front() << endl;)
        }
        //copy to _boxsize array
        _boxsize[axis] = _b_array[axis][1];
    }
    
    
    //TODOD
    // ############# ranRod Stuff ##################
    
    std::array<vector<CRod> , 3> _rodvec; // vector to store polymer rods in cell, one vector stores polymers that are parallel to the same axis

    
    void initRodsVec(){
        double xipos, xjpos;
        for (int axis=0;axis<3;axis++){//axis 0 is x axis.
            cellInterval_ai = - _b_array[i][0]
            cellInterval_aj = - _b_array[j][0]
            for (int abc=0;abc<3;abc++){
                for (int def=0;def<3;def++){
                    cellInterval_bi += _b_array[i][abc]
                    cellInterval_bj += _b_array[j][def]
                    xipos = atob(cellInterval_ai,cellInterval_bi);
                    xjpos = atob(cellInterval_aj,cellInterval_bj);
                    if ((abc ==1) && (def == 1)){
                        // In the central cell, the polymer goes to the origin, so that the particle has space to fit
                        xipos = 0;
                        xjpos = 0;
                    }
                    if ((abc ==1) && (def == 2)){
                        // The particle should have at least one escape path, so that it does net get stuck right from the start
                        xipos = 0;
                        xjpos = _b_array[j][1] + _b_array[j][2];
                    }
                    _rodarr[axis][abc][def] = CRod(axis, xipos, xjpos );
                }
            }
        }
    }
    
    //TODO
    void updateRodsVec(int crossaxis,int exitmarker){//exitmarker is -1 for negative direction, or 1 for positive
        //delete all polymers orthogonal to crossaxis, that are outside the box now
        //update other polymer positions
        int ortho[2] = {1,2};
        if (crossaxis == 1){
            ortho[0]=2; 
            ortho[1]=0;
        }
        else if (crossaxis == 2){
            ortho[0]=0; 
            ortho[1]=1;
        }
        // shift positions of rods
        int plane;
        CRod tmpRod = CRod();
        for (int oa=0;oa<2;oa++){
            plane = ortho[oa];
            //cout << "plane " << plane << endl;
            int nrods = _rodarr[plane].size();
            for (int i=0;i<nrods;i++){//need to count beckwards, due to erase function!
                //shift rod positions parallel to crossaxis. plane is direction that the shifted rods are parallel to.
                if (exitmarker == -1 && _rodarr[plane][i].coord[crossaxis] - _boxsize/2.  ) > 1.5*_boxsize){
                    // erase rods that have left the simulation box.
                    //cout << _rodvec[plane].size() << " XXX ";
                    _rodvec[plane].erase(_rodvec[plane].begin() + i);
                    //cout << _rodvec[plane].size() << endl;
                }
            }
            double reln = _n_rods / n_tries;// _n_rods / n_tries is probability of placing one of n_tries new
            assert((reln < 1.) && "Error: reln in updateRodsVec must be smaller than 1!");
            for (int j=0;j<3*n_tries;j++){// factor 3, since I reassign 3 cells per plane
                if (zerotoone() < reln ){  
                    tmpRod = CRod(plane,0.,0.);//Reset
                    //in direction parallel to crossaxis, choose new position in side cell 
                    tmpRod.coord[crossaxis] = (zerotoone()  + exitmarker) * _boxsize;
                    int ortho2 = 3 - (plane + crossaxis);
                    // in direction orthogonal to both plane and crossaxis
                    tmpRod.coord[ortho2] = (zerotoone() * 3 -1) * _boxsize;
                    _rodvec[plane].push_back(tmpRod);
                    //cout << _rodvec[plane].size() << endl;
                }
            }
        }
        avrods += _rodvec[0].size() + _rodvec[1].size() + _rodvec[2].size();
        avcount += 1;
    }
    void printAvRods(){
        cout << "nrods in yz plane mean: " << avrods/(3*avcount) << endl;
    }



private:
    void setRanNumberGen(double seed);
    void countWallCrossing(int crossaxis, int exitmarker);
    void calculateExpHPI(const double r, double& U, double& Fr);
    void calculateExpPotential(const double r, double& U, double& Fr);
    void calculateExpPotentialMOD(const double r, double& U, double& Fr, int index);
    void modifyPot(double& U, double& Fr, double dist);
    void calcLJPot(const double r, double &U, double &dU);
    void initPosHisto();





public:
    CConfiguration();
    CConfiguration(
        string distribution,double timestep,  double potRange,  double potStrength, const bool potMod,
        double psize, const bool posHisto, const bool steric, const bool ranU, bool hpi);
    void updateStartpos();
    void makeStep();
    void checkBoxCrossing();
    void calcStochasticForces();
    void calcMobilityForces();
    void calc_1RODONLY_MobilityForces();
    void calc_ZRODONLY_MobilityForces();
    void calc_YZRODONLY_MobilityForces();
    void saveXYZTraj(string name, const int& move, string a_w);

    double getUpot(){ return _upot; }
    double getDisplacement();
    unsigned int getwallcrossings(int i){ return _wallcrossings[i]; }
    bool testOverlap();
    void moveBack();
    //void addHistoValue();
    //void printHistoMatrix(string folder);
    //void positionHistogram(double block, double possq, double pposprev, int l, int posHisto[]);
    std::vector<double> getppos(){ // returns pointer to current particle position array
    	std::vector<double> pos (3);
    	for (int i = 0; i < 3; i++){
            pos[i] = _ppos[i] + _boxCoord[i];
    	}
    	return pos;
    }





};



#endif /* CCONFIGURATION_H_ */
