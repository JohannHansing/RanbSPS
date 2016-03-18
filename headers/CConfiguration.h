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


#define ifdebug(x) x

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
	double _hpi_u;
	double _hpi_k;


    //EXPONENTIAL Potential
    double _potRange;         // Avoid values like 10/3 = 3.33333... ALWAYS define the Range of the exp potential through boxsize due to scaling!
    double _potStrength;      // rescaled U_0 for exponential Potential
    double _rodDistance;
    CPolymers _poly;


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
    vector<vector<vector<int> > > _posHistoM;

    int _min, _max;        // parameters for determining up to which order neighboring rods are considered for the potential

    //Particle parameters
    double _ppos[3];    //initialize particle position (DO IT LIKE resetpos FOR MOVEPARTICLEFOLLOW/RESET)
    double _upot;
    double _f_mob[3];   //store mobility and stochastic force
    double _f_sto[3];


    boost::mt19937 *m_igen;                      //generate instance of random number generator "twister".
    double zerotoone(){
        boost::uniform_01<boost::mt19937&> dist(*m_igen);
        return dist();
    }



    int ran_sign(){
    // the variate generator uses _igen (int rand number generator),
    // samples from uniform integer distribution 0, 1
        boost::variate_generator<boost::mt19937&, boost::uniform_int<>> zeroone(*m_igen, boost::uniform_int<>(0, 1));
	return (zeroone() * 2) - 1; //this calculation makes value either -1 or 1
    }
    
    double ran_gamma(double alpha = 5., double beta = 1.){
        // http://www.boost.org/doc/libs/1_58_0/doc/html/boost/random/gamma_distribution.html
        // Mean of gamma dist is alpha* beta
        boost::variate_generator<boost::mt19937&    , boost::gamma_distribution<double> > ran_gen(
                *m_igen, boost::gamma_distribution<double>(alpha, beta));
        return ran_gen();
    }

    void initRanb(){
        // For now, I just use a gamma distribution
        for (int i=0;i<3;i++){
            _boxsize[i] = 10;
            _b_array[i][1] = _boxsize[i];
            _b_array[i][0] = ran_gamma();
            _b_array[i][2] = ran_gamma();
        }
    }
    
    void updateRanb(int axis, int exitmarker){
        // ROTATE http://en.cppreference.com/w/cpp/algorithm/rotate
        //copy neighbor boxsize to _boxsize for particle box.
        if (exitmarker=1){
            // rotation to the left
            std::rotate(_b_array[axis].begin(), _b_array[axis].begin() + 1, _b_array[axis].end());
            // assign new boxsize on right side
            _b_array[axis].back() = ran_gamma();
        }
        else {
            // rotation to the  right
            std::rotate(_b_array[axis].rbegin(), _b_array[axis].rbegin() + 1, _b_array[axis].rend());
            // assign new boxsize on right side
            ifdebug(cout << "B4 _b_array[axis].front() = " << _b_array[axis].front() << endl;)
            _b_array[axis].front() = ran_gamma();
            ifdebug(cout << "AFTER _b_array[axis].front() = " << _b_array[axis].front() << endl;)
        }
        //copy to _boxsize array
        _boxsize[axis] = _b_array[axis][1];
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
            double timestep,  double potRange,  double potStrength, double rodDistance, const bool potMod, double psize,
            const bool posHisto, const bool steric, const bool ranU,  bool hpi, double hpi_u, double hpi_k, bool ranRods, double n_rods);
    void resetParameters(double timestep, double potRange, double potStrength);
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
