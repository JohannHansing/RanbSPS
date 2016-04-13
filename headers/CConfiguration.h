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
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include "CPolymers.h"
#include "CRod.h"

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
    std::array<std::array<double,3>,3> _b_array_prev;
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
        ifdebug(if (a >= b) cout << "a " << a << "\nb " << b << endl;)
        boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<>> ran_gen(*m_igen, boost::random::uniform_real_distribution<>(a, b));
        return ran_gen();
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

    void testgamma(){//test function to check functionality of boost gamma generator --> It works!
        ofstream trajectoryfile;
        trajectoryfile.open(("tmp_gamma.txt"));
        for (int i=0;i<100000;i++){
            trajectoryfile << fixed << ran_gamma() << endl;
        }
        trajectoryfile.close();
    }

    double new_b(){
        // Function to return new random boxsize b - this assures that the new b is larger than a minimum value to avoid problems with the particle leaving the simulation box.
        double newb = ran_gamma();
        while ( newb < 1.5){
            newb = ran_gamma();
        }
        return newb;
    }

    void initRanb(){
        // For now, I just use a gamma distribution
        for (int i=0;i<3;i++){
            if (_pradius < 4.5)  _boxsize[i] = 10;
            else _boxsize[i] = 2*_pradius + 1;
            _b_array[i][1] = _boxsize[i];
            _b_array[i][0] = new_b();
            _b_array[i][2] = new_b();
            ifdebug(cout << _b_array[i][0] << "  " << _b_array[i][1] << "  " << _b_array[i][2] << endl;)
        }
    }

    void updateRanb(int axis, int exitmarker){
        // ROTATE http://en.cppreference.com/w/cpp/algorithm/rotate
        //copy neighbor boxsize to _boxsize for particle box.
        ifdebug(cout << "Update Ranb\naxis " << axis << "  --  exitmarker " << exitmarker << endl;)
        double newb = new_b();
        ifdebug(cout << "*" << newb << endl;)
        _b_array_prev = _b_array;
        if (exitmarker==1){
            // rotation to the left
            std::rotate(_b_array[axis].begin(), _b_array[axis].begin() + 1, _b_array[axis].end());
            // assign new boxsize on right side
            _b_array[axis].at(2) = newb;
        }
        else {
            // rotation to the  right
            std::rotate(_b_array[axis].rbegin(), _b_array[axis].rbegin() + 1, _b_array[axis].rend());
            // assign new boxsize on right side
            //ifdebug(cout << "B4 _b_array[axis].front() = " << _b_array[axis].front() << endl;)
            _b_array[axis].at(0) = newb;
            //ifdebug(cout << "AFTER _b_array[axis].front() = " << _b_array[axis].front() << endl;)
        }
        //copy to _boxsize array
        _boxsize[axis] = _b_array[axis][1];
    }


    //TODOD
    // ############# ranRod Stuff ##################
    std::array<std::array<std::array<CRod, 3>, 3>, 3> _rodarr; // vector to store polymer rods in cell, one vector stores polymers that are parallel to the same axis


    void initRodsArr(){
        int i,j;
        double xipos, xjpos, cellInterval_ai, cellInterval_aj;
        for (int axis=0;axis<3;axis++){//axis 0 is x axis.
            i = axis +1;
            if (i==3) i=0;
            j=3-(i+axis);
            cellInterval_ai = - _b_array[i][0];
            for (int abc=0;abc<3;abc++){//for i = axis + 1
                cellInterval_aj = - _b_array[j][0];
                for (int def=0;def<3;def++){// for j = axis + 2
                    xipos = atob(cellInterval_ai,cellInterval_ai+_b_array[i][abc]);
                    xjpos = atob(cellInterval_aj,cellInterval_aj+_b_array[j][def]);
                    if ((abc ==1) && (def == 1)){
                        // In the central cell, the polymer goes to the origin, so that the particle has space to fit
                        xipos = 0;
                        xjpos = 0;
                    }
                    //TODO
                    // if ((abc ==1) && (def == 2)){
//                         // The particle should have at least one escape path, so that it does net get stuck right from the start
//                         xipos = 0;
//                         xjpos = _b_array[j][1] + _b_array[j][2];
//                     }
                    _rodarr[axis][abc][def] = CRod(axis, xipos, xjpos );
                    cellInterval_aj += _b_array[j][def];
                }
                cellInterval_ai += _b_array[i][abc];
            }
            ifdebug(
                cout << axis << endl;
                prinRodPos(axis);
            )
        }
    }

    void initRodsRel(){
        int i,j;
        double xipos, xjpos, cellInterval_ai, cellInterval_aj;
        for (int axis=0;axis<3;axis++){//axis 0 is x axis.
            i = axis +1;
            if (i==3) i=0;
            j=3-(i+axis);
            for (int abc=0;abc<3;abc++){//for i = axis + 1
                for (int def=0;def<3;def++){// for j = axis + 2
                    xipos = atob(0,_b_array[i][abc]);
                    xjpos = atob(0,_b_array[j][def]);
                    if ((abc ==1) && (def == 1)){
                        // In the central cell, the polymer goes to the origin, so that the particle has space to fit
                        xipos = 0;
                        xjpos = 0;
                    }
                    //TODO
                    // if ((abc ==1) && (def == 2)){
//                         // The particle should have at least one escape path, so that it does net get stuck right from the start
//                         xipos = 0;
//                         xjpos = _b_array[j][1] + _b_array[j][2];
//                     }
                    _rodarr[axis][abc][def] = CRod(axis, xipos, xjpos );
                }
            }
            ifdebug(
                cout << axis << endl;
                prinRodPos(axis);
            )
        }
    }

    //TODO
    void updateRodsArr(int crossaxis,int exitmarker){//exitmarker is -1 for negative direction, or 1 for positive
        //delete all polymers orthogonal to crossaxis, that are outside the box now
        //update other polymer positions
        bool overlaps;
        double cellInterval_ai, cellInterval_aj;
        int i,j;
        i=crossaxis+1;
        if (i==3) i =0;
        j=3-(i+crossaxis);
        // rotate around rods in cells abc and def and reassign
        if (exitmarker == 1){
            // shift positions of rods USING PREVIOUS b ARRAY _b_array_prev!
            for (int abc=0; abc<3;abc++){
                for (int def=0; def<3;def++){
                    _rodarr[i][abc][def].coord[crossaxis] -= _b_array_prev[crossaxis][1];
                    _rodarr[j][abc][def].coord[crossaxis] -= _b_array_prev[crossaxis][1];
                }
            }
            rotate_left(_rodarr[j]);
            cellInterval_ai = - _b_array[i][0];
            cellInterval_aj = - _b_array[j][0];
            for (int abc=0;abc<3;abc++){
                rotate_left(_rodarr[i][abc]);
                // new rod positions
                //Example: -_b_array[j][0] , 0
                   //      0 , _b_array[j][1]
                   //      _b_array[j][1], _b_array[j][1]+ _b_array[j][2]
                overlaps=true;
                while (overlaps){
                    _rodarr[i][abc][2].coord[crossaxis] = atob(_b_array[crossaxis][1] , _b_array[crossaxis][1]+_b_array[crossaxis][2]);
                    _rodarr[i][abc][2].coord[j] = atob(cellInterval_aj , cellInterval_aj+_b_array[j][abc]);
                    overlaps= testTracerOverlap(crossaxis, j, _rodarr[i][abc][2].coord[crossaxis], _rodarr[i][abc][2].coord[j]);
                    //cout << "Repeat?";
                }
                overlaps=true;
                while (overlaps){
                    _rodarr[j][2][abc].coord[crossaxis] = atob(_b_array[crossaxis][1] , _b_array[crossaxis][1]+_b_array[crossaxis][2]);
                    _rodarr[j][2][abc].coord[i] = atob(cellInterval_ai , cellInterval_ai+_b_array[i][abc]);
                    overlaps= testTracerOverlap(crossaxis, i, _rodarr[j][2][abc].coord[crossaxis], _rodarr[j][2][abc].coord[i]);
                }
                cellInterval_aj+=_b_array[j][abc];
                cellInterval_ai+=_b_array[i][abc];
            }
        }
        else{
            // shift positions of rods
            for (int abc=0; abc<3;abc++){
                for (int def=0; def<3;def++){
                    _rodarr[i][abc][def].coord[crossaxis] += _b_array_prev[crossaxis][0];
                    _rodarr[j][abc][def].coord[crossaxis] += _b_array_prev[crossaxis][0];
                }
            }
            rotate_right(_rodarr[j]);
            cellInterval_ai = - _b_array[i][0];
            cellInterval_aj = - _b_array[j][0];
            for (int abc=0;abc<3;abc++){
                rotate_right(_rodarr[i][abc]);
                // new rod positions
                overlaps=true;
                while (overlaps){
                    _rodarr[i][abc][0].coord[crossaxis] = atob(-_b_array[crossaxis][0] , 0);
                    _rodarr[i][abc][0].coord[j] = atob(cellInterval_aj , cellInterval_aj+_b_array[j][abc]);
                    overlaps= testTracerOverlap(crossaxis, j, _rodarr[i][abc][0].coord[crossaxis], _rodarr[i][abc][0].coord[j]);
                }
                overlaps=true;
                while (overlaps){
                    _rodarr[j][0][abc].coord[crossaxis] = atob(-_b_array[crossaxis][0] , 0);
                    _rodarr[j][0][abc].coord[i] = atob(cellInterval_ai , cellInterval_ai+_b_array[i][abc]);
                    overlaps= testTracerOverlap(crossaxis, i, _rodarr[j][0][abc].coord[crossaxis], _rodarr[j][0][abc].coord[i]);
                }
                cellInterval_aj+=_b_array[j][abc];
                cellInterval_ai+=_b_array[i][abc];
            }
        }
        ifdebug(
            if (testOverlap()){
                cout << "\nERROR still overlap after newrod init!" << endl;
                //abort();
            }
        )
    }



    void updateRodsRel(int crossaxis,int exitmarker){//exitmarker is -1 for negative direction, or 1 for positive
        //delete all polymers orthogonal to crossaxis, that are outside the box now
        //update other polymer positions
        //TODO del
        // cout << "update rods\ncrossaxis " << crossaxis << " -- exitm " << exitmarker << endl;
        bool overlaps;
        double cellInterval_ai, cellInterval_aj;
        int i,j;
        i=crossaxis+1;
        if (i==3) i =0;
        j=3-(i+crossaxis);
        // rotate around rods in cells abc and def and reassign
        if (exitmarker == 1){
            rotate_left(_rodarr[j]);
            for (int abc=0;abc<3;abc++){
                rotate_left(_rodarr[i][abc]);
                // new rod positions
                //Example: -_b_array[j][0] , 0
                   //      0 , _b_array[j][1]
                   //      _b_array[j][1], _b_array[j][1]+ _b_array[j][2]
                overlaps=true;
                while (overlaps){
                    _rodarr[i][abc][2].coord[crossaxis] = atob(0,_b_array[crossaxis][2]);
                    _rodarr[i][abc][2].coord[j] = atob(0, _b_array[j][abc]);
                    //TODO overlaps= testTracerOverlap(crossaxis, j, _rodarr[i][abc][2].coord[crossaxis], _rodarr[i][abc][2].coord[j]);
                    overlaps=false;
                    //cout << "Repeat?";
                }
                overlaps=true;
                while (overlaps){
                    _rodarr[j][2][abc].coord[crossaxis] = atob(0, _b_array[crossaxis][2]);
                    _rodarr[j][2][abc].coord[i] = atob(0, _b_array[i][abc]);
                    //TODO overlaps= testTracerOverlap(crossaxis, i, _rodarr[j][2][abc].coord[crossaxis], _rodarr[j][2][abc].coord[i]);
                    overlaps=false;
                }
            }
        }
        else{
            rotate_right(_rodarr[j]);
            for (int abc=0;abc<3;abc++){
                rotate_right(_rodarr[i][abc]);
                // new rod positions
                overlaps=true;
                while (overlaps){
                    _rodarr[i][abc][0].coord[crossaxis] = atob(0.,_b_array[crossaxis][0]);
                    _rodarr[i][abc][0].coord[j] = atob(0.,_b_array[j][abc]);
                    //TODO overlaps= testTracerOverlap(crossaxis, j, _rodarr[i][abc][0].coord[crossaxis], _rodarr[i][abc][0].coord[j]);
                    overlaps=false;
                }
                overlaps=true;
                while (overlaps){
                    _rodarr[j][0][abc].coord[crossaxis] = atob(0.,_b_array[crossaxis][0]);
                    _rodarr[j][0][abc].coord[i] = atob(0.,_b_array[i][abc]);
                    //TODO overlaps= testTracerOverlap(crossaxis, i, -_b_array[crossaxis][0]+_rodarr[j][0][abc].coord[crossaxis], _rodarr[j][0][abc].coord[i]);
                    overlaps=false;
                }
            }
        }
    }

    bool testTracerOverlap(int i, int j, double ri, double rj){
        if ((pow( _ppos[i] - ri , 2 ) + pow( _ppos[j] - rj , 2 )) < _r_cSq){
            ifdebug(cout << "Overlap distance " << pow( _ppos[i] - ri , 2 ) + pow( _ppos[j] - rj , 2 ) << endl;)
            return true;
        }
        return false;
    }

    void prinRodPos(int axis){
        for (int irod=0;irod<_rodarr[axis].size();irod++){
            for (int jrod=0;jrod<_rodarr[axis].size();jrod++){
                double rx =_rodarr[axis][irod][jrod].coord[0];
                double ry =_rodarr[axis][irod][jrod].coord[1];
                double rz =_rodarr[axis][irod][jrod].coord[2];
                cout << ",[" << rx << "," << ry << "," << rz << "]";
            }
        }
        cout << "]," << endl;
    }

    bool tracerInBox(){
        for (int ax = 0;ax<3;ax++){
            if ((_ppos[ax] < 0.) || (_ppos[ax] > _boxsize[ax])){
                cout << "\nax " << ax << "\n_boxsize[ax] " << _boxsize[ax] << "\n_ppos[ax] " << _ppos[ax] << endl;
                return false;
            }
        }
        return true;
    }

    // bool rodinCell(){
//         for (int axis = 0;axis<3;axis++){
//             int i,j;
//             i=axis+1;
//             if (i==3) i =0;
//             j=3-(i+axis);
//             for (int irod=0;irod<_rodarr[axis].size();irod++){
//                 for (int jrod=0;jrod<_rodarr[axis].size();jrod++){
//                     double rax =_rodarr[axis][irod][jrod].coord[axis];
//                     double ri =_rodarr[axis][irod][jrod].coord[i];
//                     double rj =_rodarr[axis][irod][jrod].coord[j];
//                     if  (rax != 0.){
//                         cout << "\nax " << ax <<"\nError rax not zero, but rax = " << rax << endl;
//                         return false;
//                     }
//                     if  TODO ((ri < _b_array[i][irod]) || (_ppos[ax] > _boxsize[ax])){
//                         cout << "\nax " << ax << "\n_boxsize[ax] " << _boxsize[ax] << "\n_ppos[ax] " << _ppos[ax] << endl;
//                         return false;
//                     }
//                 }
//             }
//         }
//         return true;
//     }



    template<typename T, size_t N>
    void rotate_left(std::array<T,N> & arr){
        // rotation to the left
        std::rotate(arr.begin(), arr.begin() + 1, arr.end());
    }
    template<typename T, size_t N>
    void rotate_right(std::array<T,N> & arr){
        // rotation to the  right
        std::rotate(arr.rbegin(), arr.rbegin() + 1, arr.rend());
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
        double psize, const bool posHisto, const bool steric, const bool ranU, bool hpi, bool ranRod);
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
    std::vector<double> getppos_rel(){ // returns pointer to current particle position array
    	std::vector<double> pos (3);
    	for (int i = 0; i < 3; i++){
            pos[i] = _ppos[i];
    	}
    	return pos;
    }





};



#endif /* CCONFIGURATION_H_ */
