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
#include <Eigen/Dense>
#include "CRod.h"


#define ifdebug(x) 

using namespace std;

class CConfiguration {
    /*Class where all the configuration variables such as potRange etc. and also most functions for the
     * simulation are stored
     */
private:
    //MISC
    FILE* m_traj_file;
    
    //SCALING
    double _timestep;         //This is the RESCALED timestep! timestep = dt * kT / (frictionCoeffcient * particlesize)
    double _mu_sto;
    double _pradius;     //particle size is most relevant for scaling! (default = 1)
    double _boxsize[3];          // ALWAYS define boxsize through particlesize due to scaling!
    double _bdef = 10;     //default boxsize
    std::array<std::array<double,3>,3> _b_array;
    std::array<std::array<double,3>,3> _b_array_prev;
    double _epsilon;

    //EXPONENTIAL Potential
    double _potRange;         // Avoid values like 10/3 = 3.33333... ALWAYS define the Range of the exp potential through boxsize due to scaling!
    double _potStrength;      // rescaled U_0 for exponential Potential
    double _rodDistance;
    double _cutoffExpSq;
    double _cylLJSq;    //square of steric interaction parameter between cylinder and tracer particle _pradius + _polyrad
    
    


    //bool Parameters
    bool _potMod;             // true if the modified exponential potential version is used that is not 3*U_0 at the intersections.
    bool _LJPot;              // if true, then LJ interaction is calculated for steric hindrance
    bool _ranU;
    bool _hpi;
    bool _ranRod;
    bool _fixb;
    bool _rand;

    //COUNTERS AND INIT VALUES
    double _boxCoord[3];
    unsigned int _wallcrossings[3]; //counts how many times the particle has crossed the [entry, opposite, side] wall, while
                                            //travelling through the lattice
    int _entryside[3];            //records through which side and in which direction the particle last entered a box. For the face of the cube in the x-z
                                            //plane at y=0, it is entryside[0] = 1, in the x-y plane at z=L it is entryside[1] = -1!
    double _resetpos;
    double _startpos[3];          //Stores where the particle starting position was. This is needed to calculate the mean square displacement
    Eigen::Vector3d _prevpos;           //Stores previous particle position before particle is moved.
    std::array<std::array<double, 16>, 3> _distarr;

    int _min, _max;        // parameters for determining up to which order neighboring rods are considered for the potential

    //Particle parameters
    Eigen::Vector3d _ppos;    //initialize particle position (DO IT LIKE resetpos FOR MOVEPARTICLEFOLLOW/RESET)
    double _upot;
    Eigen::Vector3d _f_mob;   //store mobility and stochastic force
    Eigen::Vector3d _f_sto;
    // rod parameters
    double _polyrad;
    double _polydiamSq;
    
    
    
    //---------------------------- POTENTIALS ------------------------------
    double _uLJ;
        
    
    void calcLJPot(const double rSq, double& U, double& Fr, double stericSq){
        //Function to calculate the Lennard-Jones Potential
        double  por6 = pow((stericSq / rSq),3); //por6 stands for "p over r to the power of 6" . The 2 comes from the fact, that I need the particle radius, not the particle size
        ifdebug(if (4. * ( por6*por6 - por6 + 0.25 ) > 50) cout << "Very large LJ!!"<< endl;)
        U += 4. * ( por6*por6 - por6 + 0.25 );
        Fr +=  24. / ( rSq ) * ( 2. * por6*por6 - por6 );
    }
 


    // gamma distribution
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
    
    double getfixb(){
        // helper function to fix boxsize
        return _bdef;
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
        if (_fixb) {return _bdef;}
        
        double newb = ran_gamma();
        while ( newb < 0.5){
            newb = ran_gamma();
        }
        return newb;
    }

    void initRanb(){
        // For now, I just use a gamma distribution
        for (int i=0;i<3;i++){
            _boxsize[i] = _bdef;
            if ((_pradius > 4.5) && !_fixb )  _boxsize[i] = 2*_pradius + 1;
            _b_array[i][1] = _boxsize[i];
            _b_array[i][0] = new_b();
            _b_array[i][2] = new_b();
            ifdebug(cout << "INIT Ranb :" << _b_array[i][0] << "  " << _b_array[i][1] << "  " << _b_array[i][2] << endl;)
        }
    }

    void updateRanb(int axis, int exitmarker){
        _b_array_prev = _b_array;
        if (_fixb) {return;}
        // ROTATE http://en.cppreference.com/w/cpp/algorithm/rotate
        //copy neighbor boxsize to _boxsize for particle box.
        ifdebug(cout << "Update Ranb\naxis " << axis << "  --  exitmarker " << exitmarker << endl;)
        double newb = new_b();
        ifdebug(cout << "*" << newb << endl;)
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
        ifdebug(cout << "New b_array " << _b_array[axis][0] << "  " << _b_array[axis][1] << "  " << _b_array[axis][2] << endl;)
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
                    _rodarr[axis][abc][def] = CRod(axis, xipos, xjpos, _ranU, m_igen );
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
    
    //************** Rand ***************    
    double ran_norm(){
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > ran_gen(
                *m_igen, boost::normal_distribution<double>(0, _dvar));
        return ran_gen(); 
    }
    
    double _dvar = 1.;
   
    
    std::array<std::array<std::array<CRod, 4>, 4>, 3> _drods; // array to store polymer rods in simulation box, the outermost array stores polymers that are parallel to the same axis
    void initRand(){
        // Initialize system with random rod displacement d
        bool overlaps;
        double xipos, xjpos;
        double tmp[4] = {-_bdef,0.,_bdef,2*_bdef};
        for (int axis=0;axis<3;axis++){//axis 0 is x axis.
            int i,j;
            i=axis+1;
            if (i==3) i=0;
            j=3-(i+axis);
            for (int abcd=0;abcd<4;abcd++){
                for (int efgh=0;efgh<4;efgh++){
                    overlaps=true;
                    while (overlaps){
                        xipos = tmp[abcd] + ran_norm();
                        xjpos = tmp[efgh] + ran_norm();
                        overlaps= testTracerOverlap(i, j, xipos, xjpos);
                        //cout << "Repeat?";
                    }
                    CRod newRod = CRod(axis, xipos, xjpos, _ranU, m_igen );
                    _drods[axis][abcd][efgh] = newRod;
                }
            }
        }
        ifdebug(prinRodPos(0);)
    }
    
    void updateRand(int crossaxis,int exitmarker){
        //exitmarker is -1 for negative direction, or 1 for positive
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
            for (int abcd=0; abcd<4;abcd++){
                for (int efgh=0; efgh<4;efgh++){
                    _drods[i][abcd][efgh].coord[crossaxis] -= _bdef;
                    _drods[j][abcd][efgh].coord[crossaxis] -= _bdef;
                }
            }
            rotate_left(_drods[j]);
            cellInterval_ai = - _bdef;
            cellInterval_aj = - _bdef;
            for (int abcd=0;abcd<4;abcd++){
                rotate_left(_drods[i][abcd]);
                // new rod positions
                //Example: -_b_array[j][0] , 0
                   //      0 , _b_array[j][1]
                   //      _b_array[j][1], _b_array[j][1]+ _b_array[j][2]
                overlaps=true;
                while (overlaps){
                    _drods[i][abcd][3].coord[crossaxis] = 2*_bdef + ran_norm();
                    _drods[i][abcd][3].coord[j] = cellInterval_aj + ran_norm();
                    overlaps= testTracerOverlap(crossaxis, j, _drods[i][abcd][3].coord[crossaxis], _drods[i][abcd][3].coord[j])
                        //TODO overlap
                        || testRodOverlap(i, crossaxis, j, _drods[i][abcd][3].coord[crossaxis], _drods[i][abcd][3].coord[j]);
                    //cout << "Repeat?";
                }
                overlaps=true;
                while (overlaps){
                    _drods[j][3][abcd].coord[crossaxis] = 2*_bdef + ran_norm();
                    _drods[j][3][abcd].coord[i] = cellInterval_ai + ran_norm();
                    overlaps= testTracerOverlap(crossaxis, i, _drods[j][3][abcd].coord[crossaxis], _drods[j][3][abcd].coord[i])
                    //TODO overlap
                        || testRodOverlap(j, crossaxis, i, _drods[j][3][abcd].coord[crossaxis], _drods[j][3][abcd].coord[i]);
                }
                cellInterval_aj+=_bdef;
                cellInterval_ai+=_bdef;
            }
        }
        else{
            // shift positions of rods
            for (int abcd=0; abcd<4;abcd++){
                for (int efgh=0; efgh<4;efgh++){
                    _drods[i][abcd][efgh].coord[crossaxis] += _bdef;
                    _drods[j][abcd][efgh].coord[crossaxis] += _bdef;
                }
            }
            rotate_right(_drods[j]);
            cellInterval_ai = - _bdef;
            cellInterval_aj = - _bdef;
            for (int abcd=0;abcd<4;abcd++){
                rotate_right(_drods[i][abcd]);
                // new rod positions
                overlaps=true;
                while (overlaps){
                    _drods[i][abcd][0].coord[crossaxis] = -_bdef + ran_norm();
                    _drods[i][abcd][0].coord[j] = cellInterval_aj + ran_norm();
                    overlaps= testTracerOverlap(crossaxis, j, _drods[i][abcd][0].coord[crossaxis], _drods[i][abcd][0].coord[j])
                        //TODO overlap
                        || testRodOverlap(i, crossaxis, j, _drods[i][abcd][0].coord[crossaxis], _drods[i][abcd][0].coord[j]);
                }
                overlaps=true;
                while (overlaps){
                    _drods[j][0][abcd].coord[crossaxis] = -_bdef + ran_norm();
                    _drods[j][0][abcd].coord[i] = cellInterval_aj + ran_norm();
                    overlaps= testTracerOverlap(crossaxis, i, _drods[j][0][abcd].coord[crossaxis], _drods[j][0][abcd].coord[i])
                        //TODO overlap
                        || testRodOverlap(j, crossaxis, i, _drods[j][0][abcd].coord[crossaxis], _drods[j][0][abcd].coord[i]);
                }
                cellInterval_aj+=_bdef;
                cellInterval_ai+=_bdef;
            }
        }
        ifdebug(
            if (testOverlap()){
                cout << "\nERROR still overlap after newrod init!" << endl;
                //abort();
            }
        )
    }


    bool testTracerOverlap(int i, int j, double ri, double rj){
        //make this a little bigger, so that it's out of reach of LJ
        return ((pow( _ppos(i) - ri , 2 ) + pow( _ppos(j) - rj , 2 )) < 1.13*_cylLJSq);
    }

    // bool testTracerOverlap(int i, int j, double ri, double rj){
//         return ((pow( _ppos(i) - ri , 2 ) + pow( _ppos(j) - rj , 2 )) < 1.13*_cylLJSq);//make this a little bigger, so that it's out of reach of LJ
//     }
    
    //TODO overlap
    bool testRodOverlap(int rodaxis, int i, int j, double ri, double rj){
        double distSq;
        if (_polyrad==0) return false;
        for (auto & vrods : _drods[rodaxis]){
            for (auto & rod :  vrods){
                //rods axis need to have distance of a least a (i.e. polydiam)
                distSq = pow(rod.coord[i] - ri,2) + pow(rod.coord[j] - rj,2);
                if (  (distSq < _polydiamSq) && (distSq > 0.00001)  ){//The second clause is to avoid testing overlap with itself
                    return true;
                }
            }
        }
        return false;
    }

    void prinRodPos(int axis){
        for (int irod=0;irod<_drods[axis].size();irod++){
            for (int jrod=0;jrod<_drods[axis].size();jrod++){
                double rx =_drods[axis][irod][jrod].coord[0];
                double ry =_drods[axis][irod][jrod].coord[1];
                double rz =_drods[axis][irod][jrod].coord[2];
                cout << ",[" << rx << "," << ry << "," << rz << "]";
            }
        }
        cout << "]," << endl;
    }

    bool tracerInBox(){
        for (int ax = 0;ax<3;ax++){
            if ((_ppos(ax) < -0.1*_boxsize[ax]) || (_ppos(ax) > 1.1*_boxsize[ax])){
                cout << "\nax " << ax << "\n_boxsize[ax] " << _boxsize[ax] << "\n_ppos(ax) " << _ppos(ax) << endl;
                return false;
            }
        }
        return true;
    }




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
        string distribution,double timestep,  double potRange,  double potStrength, const bool rand,
        double psize, const bool posHisto, const bool steric, const bool ranU, bool ranRod, double dvar, double polydiam, string tmp5);
    void updateStartpos();
    void makeStep();
    void checkBoxCrossing();
    void calcStochasticForces();
    void calcMobilityForces();
    void calc_1RODONLY_MobilityForces();
    void calc_ZRODONLY_MobilityForces();
    void calc_YZRODONLY_MobilityForces();
    void saveXYZTraj(string name, int move, string flag);

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
            pos[i] = _ppos(i) + _boxCoord[i];
    	}
    	return pos;
    }
    std::vector<double> getppos_rel(){ // returns pointer to current particle position array
    	std::vector<double> pos (3);
    	for (int i = 0; i < 3; i++){
            pos[i] = _ppos(i);
    	}
    	return pos;
    }
    void writeDistances(ostream& distancesfile, unsigned int stepcount) {
    // So far this only writes the tracer particle position
	distancesfile << fixed << stepcount * _timestep << "\t";
        for (auto & arr : _distarr){
            for (auto & dist : arr){
                distancesfile << dist << " ";
            }
        }
        distancesfile << endl;
    }


// _______________________________ OLD REL STUFF _____________________________________________________


    // void updateRodsRel(int crossaxis,int exitmarker){//exitmarker is -1 for negative direction, or 1 for positive
    //     //delete all polymers orthogonal to crossaxis, that are outside the box now
    //     //update other polymer positions
    //     //TODO del
    //     // cout << "update rods\ncrossaxis " << crossaxis << " -- exitm " << exitmarker << endl;
    //     bool overlaps;
    //     double cellInterval_ai, cellInterval_aj;
    //     int i,j;
    //     i=crossaxis+1;
    //     if (i==3) i =0;
    //     j=3-(i+crossaxis);
    //     // rotate around rods in cells abc and def and reassign
    //     if (exitmarker == 1){
    //         rotate_left(_rodarr[j]);
    //         for (int abc=0;abc<3;abc++){
    //             rotate_left(_rodarr[i][abc]);
    //             // new rod positions
    //             //Example: -_b_array[j][0] , 0
    //                //      0 , _b_array[j][1]
    //                //      _b_array[j][1], _b_array[j][1]+ _b_array[j][2]
    //             overlaps=true;
    //             while (overlaps){
    //                 _rodarr[i][abc][2].coord[crossaxis] = atob(0,_b_array[crossaxis][2]);
    //                 _rodarr[i][abc][2].coord[j] = atob(0, _b_array[j][abc]);
    //                 //TODO overlaps= testTracerOverlap(crossaxis, j, _rodarr[i][abc][2].coord[crossaxis], _rodarr[i][abc][2].coord[j]);
    //                 overlaps=false;
    //                 //cout << "Repeat?";
    //             }
    //             overlaps=true;
    //             while (overlaps){
    //                 _rodarr[j][2][abc].coord[crossaxis] = atob(0, _b_array[crossaxis][2]);
    //                 _rodarr[j][2][abc].coord[i] = atob(0, _b_array[i][abc]);
    //                 //TODO overlaps= testTracerOverlap(crossaxis, i, _rodarr[j][2][abc].coord[crossaxis], _rodarr[j][2][abc].coord[i]);
    //                 overlaps=false;
    //             }
    //         }
    //     }
    //     else{
    //         rotate_right(_rodarr[j]);
    //         for (int abc=0;abc<3;abc++){
    //             rotate_right(_rodarr[i][abc]);
    //             // new rod positions
    //             overlaps=true;
    //             while (overlaps){
    //                 _rodarr[i][abc][0].coord[crossaxis] = atob(0.,_b_array[crossaxis][0]);
    //                 _rodarr[i][abc][0].coord[j] = atob(0.,_b_array[j][abc]);
    //                 //TODO overlaps= testTracerOverlap(crossaxis, j, _rodarr[i][abc][0].coord[crossaxis], _rodarr[i][abc][0].coord[j]);
    //                 overlaps=false;
    //             }
    //             overlaps=true;
    //             while (overlaps){
    //                 _rodarr[j][0][abc].coord[crossaxis] = atob(0.,_b_array[crossaxis][0]);
    //                 _rodarr[j][0][abc].coord[i] = atob(0.,_b_array[i][abc]);
    //                 //TODO overlaps= testTracerOverlap(crossaxis, i, -_b_array[crossaxis][0]+_rodarr[j][0][abc].coord[crossaxis], _rodarr[j][0][abc].coord[i]);
    //                 overlaps=false;
    //             }
    //         }
    //     }
    // }

//     void initRodsRel(){
//         int i,j;
//         double xipos, xjpos, cellInterval_ai, cellInterval_aj;
//         for (int axis=0;axis<3;axis++){//axis 0 is x axis.
//             i = axis +1;
//             if (i==3) i=0;
//             j=3-(i+axis);
//             for (int abc=0;abc<3;abc++){//for i = axis + 1
//                 for (int def=0;def<3;def++){// for j = axis + 2
//                     xipos = atob(0,_b_array[i][abc]);
//                     xjpos = atob(0,_b_array[j][def]);
//                     if ((abc ==1) && (def == 1)){
//                         // In the central cell, the polymer goes to the origin, so that the particle has space to fit
//                         xipos = 0;
//                         xjpos = 0;
//                     }
//                     //TODO
//                     // if ((abc ==1) && (def == 2)){
// //                         // The particle should have at least one escape path, so that it does net get stuck right from the start
// //                         xipos = 0;
// //                         xjpos = _b_array[j][1] + _b_array[j][2];
// //                     }
//                     _rodarr[axis][abc][def] = CRod(axis, xipos, xjpos );
//                 }
//             }
//             ifdebug(
//                 cout << axis << endl;
//                 prinRodPos(axis);
//             )
//         }
//     }


};



#endif /* CCONFIGURATION_H_ */
