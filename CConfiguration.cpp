#include "headers/CConfiguration.h"


using namespace std;



CConfiguration::CConfiguration(){
}

CConfiguration::CConfiguration(
        string distribution, double timestep,  double potRange,  double potStrength, const bool potMod,
        double psize, const bool posHisto, const bool steric, const bool ranU, bool hpi, bool ranRod){
    setRanNumberGen(0);
    _potRange = potRange;
    _potStrength = potStrength;
    _pradius = psize/2;
    _r_cSq = pow(1.122462 * _pradius,2);
    _cutoffExpSq = 5*_potRange;
    _timestep = timestep;
    _ranRod = ranRod;
    _potMod = potMod;
    _LJPot = (steric == false) && (psize != 0);
    _ranU = ranU;
    _poly = CPolymers();
    _upot = 0;
    _mu_sto = sqrt( 2 * _timestep );                 //timestep for stochastic force
    for (int i = 0; i < 3; i++){
        _ppos[i] = 5;
        _startpos[i] = 5;
        _entryside[i] = 0;
        _wallcrossings[i] = 0;
        _boxCoord[i] = 0;
        _prevpos[i] = 5;
    }
   if (distribution == "gamma2"){
        _alpha = 5.; _beta = 2.;
       cout << "_alpha = " << _alpha << endl;
   }
   if (distribution == "gamma4"){
        _alpha = 2.5; _beta = 4.;
       cout << "_alpha = " << _alpha << endl;
   }

    if (_ranU) {
       cout << "TODOOOOO" << endl;
    }
    //TODO del
    //testgamma();
    initRanb();
    
    if (_ranRod){
        initRodsArr();
    }
}

void CConfiguration::updateStartpos(){
    //This function is used if the particle should keep moving after each run, and not start at _resetpos again, like in the first run
    //This (hopefully) will give better averages without having to spend a lot of steps in the beginning of each run to get away from _resetpos
    for (int i = 0; i < 3; i++){
    _startpos[i] = _ppos[i] + _boxCoord[i];
    }
}



void CConfiguration::makeStep(){
    //move the particle according to the forces and record trajectory like watched by outsider
    for (int i = 0; i < 3; i++){
        _prevpos[i] = _ppos[i];
        _ppos[i] += _timestep * _f_mob[i] + _mu_sto * _f_sto[i];
        if (std::isnan(_ppos[i])){
            cout << "axis " << i << "\n_prevpos " << _prevpos[i] << "\n_ppos " << _ppos[i] << endl;
            cout << "_f_mob[i] " << _f_mob[i] << "\n_f_sto[i] " << _f_sto[i] << endl;
            abort();
        }
    }
}

void CConfiguration::checkBoxCrossing(){
    //should the particle cross the confinement of the cube, let it appear on the opposite side of the box
    int exitmarker = 0;
    for (int i = 0; i < 3; i++){
        exitmarker =0;
        if (_ppos[i] < 0){
            MISTAKE! _ppos[i] += _boxsize[i];
            _boxCoord[i] -= _boxsize[i];
            exitmarker = -1;
        }
        else if (_ppos[i] > _boxsize[i]){
            _ppos[i] -= _boxsize[i];
            _boxCoord[i] += _boxsize[i];
            exitmarker = 1;
            if (_ppos[i] > 10){
                cout << "Bad _ppos after boxcrossing. Aborting!" << endl;
                cout << "_ppos b4 " << _ppos[i] + _boxsize[i] << endl;
                cout << "_ppos after " << _ppos[i] << endl;
                abort();
            }
        }
        if (exitmarker!=0){
            updateRanb(i,exitmarker);
            if (_ranRod){
                updateRodsArr(i, exitmarker);
                ifdebug(
                    cout << "[["<<i<<"," << exitmarker <<"] ";
                    prinRodPos(0); // cout print rod pos!
                )
            }
            if (_ranU) _poly.shiftPolySign(i, exitmarker);
        }
    }
    //TODO del
    if (!tracerInBox()){
        cout << "\nTRACER NOT IN BOX!!!" << endl;
        abort();
    }
}






void CConfiguration::calcStochasticForces(){

    // the variate generator uses m_igen (int rand number generator),
    // samples from normal distribution with variance 1 (later sqrt(2) is multiplied)
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > ran_gen(
            *m_igen, boost::normal_distribution<double>(0, 1));

    for (int i = 0; i < 3; i++) {
        _f_sto[i] = ran_gen();
    }
}


void CConfiguration::calcMobilityForces(){
    //calculate mobility forces from potential Epot - Unified version that includes 2ndOrder if k is larger than or equal 0.2 b , except if ranPot is activated.
    double r_i = 0, r_k = 0;
    array<double,4> r_is, r_ks;
    vector<double> r_absSq;
    double utmp = 0, frtmp = 0;     //temporary "hilfsvariables"
    double Epot = 0;
    double z1, z2;
    if (_ranU){
        cout << "TODOOooOoOOoO calcMob for RanU" << endl;
        //z1 = 1/4 * _boxsize;
        //z2 = _boxsize - z1;   //z is in cylindrical coordinates. This indicates above/below which value the exp potential is modifed for random signs.
    }
    //reset mobility forces to zero
    for (int l = 0; l < 3; l++) {
        _f_mob[l] = 0;
    }

    for (int i = 0; i < 3; i++){
        int k = i + 1;   //k always one direction "further", i.e. if i = 0 = x-direction, then k = 1 = y-direction
        if ( k == 3 ) k = 0;
        int plane = 3 - (i+k); //this is the current plane of the cylindrical coordinates
        int n = 0;     // reset counter for index of next rod in plane  n = 0, 1, 2, 3 -> only needed for ranPot
        
        if (_ranRod){
            for (int abc=0;abc<3;abc++){
                for (int def=0;def<3;def++){
                    r_i = _ppos[i] - _rodarr[plane][abc][def].coord[i];
                    r_k = _ppos[k] - _rodarr[plane][abc][def].coord[k];
                    r_absSq.push_back( r_i * r_i + r_k * r_k);
                }
            }
        }
        else{
            //This creates the distance vectors from the rods to the tracer
            r_ks[0] = _ppos[k] + _b_array[k][0]; 
            r_is[0] = _ppos[i] + _b_array[i][0]; 
            for (int rodi = 0; rodi < 3; rodi++){
                r_ks[rodi+1] = r_ks[rodi] - _b_array[k][rodi]; 
                r_is[rodi+1] = r_is[rodi] - _b_array[i][rodi]; 
            }
            for (int nk = 0; nk < r_is.size(); nk++){
                for (int ni = 0; ni < r_ks.size(); ni++){
                    r_i = r_is[ni];
                    r_k = r_ks[nk];

                    r_absSq.push_back( r_i * r_i + r_k * r_k); //distance to the rods
                    //if (r_abs < 0.9*_pradius) cout << "Small rho distace: " << r_abs << endl;
                }
            }
        }
        
        for (int j=0;j<r_absSq.size();j++){
            calculateExpPotential(r_absSq.at(j), utmp, frtmp);


            // if (_ranU){
//                     utmp = utmp * _poly.get_sign(plane, n);
//                     frtmp = frtmp * _poly.get_sign(plane, n);
//                     if (_ppos[plane] > z2){
//                         if (! _poly.samesign(1, plane, n)){
//                             _f_mob[plane] += utmp * 4 / _boxsize;              //this takes care of the derivative of the potential modification and resulting force
//                             modifyPot(utmp, frtmp, _boxsize - _ppos[plane]);
//                         }
//                     }
//                     else if (_ppos[plane] < z1){
//                         if (! _poly.samesign(-1, plane, n)){
//                             _f_mob[plane] -= utmp * 4 / _boxsize;              //this takes care of the derivative of the potential modification and resulting force
//                             modifyPot(utmp, frtmp, _ppos[plane]);
//                         }
//                     }
//                     n++;  //index of next rod in curent plane
//                 }
            
            if (_LJPot && ( r_absSq.at(j) < _r_cSq )) calcLJPot(r_absSq.at(j), utmp, frtmp);
            
            //TODO del
            if (utmp > 10){
                cout << "utmp " << utmp << "\nrSq " << r_absSq.at(j) << "\nindex j "  << j << endl;
            }
            


            Epot += utmp;
            _f_mob[i] += frtmp * r_i;
            _f_mob[k] += frtmp * r_k;
        }
    }
    _upot = Epot;
}


void CConfiguration::saveXYZTraj(string name, const int& move, string a_w){
    // FILE *f = fopen(name.c_str(), a_w.c_str());
//
//     fprintf(f, "%d\n%s (%8.3f %8.3f %8.3f) t=%d \n", 1, "sim_name", _boxsize, _boxsize, _boxsize, move);
//
//     // polymer particles
//     fprintf(f, "%3s%9.3f%9.3f%9.3f \n","H", _ppos[0]+ _boxsize *_boxnumberXYZ[0], _ppos[1]+ _boxsize *_boxnumberXYZ[1], _ppos[2]+ _boxsize *_boxnumberXYZ[2]);   // absolute position
//     // fprintf(f, "%3s%9.3f%9.3f%9.3f \n","H", _ppos[0], _ppos[1], _ppos[2]);  // relative position in box
//     fclose(f);
}


void CConfiguration::setRanNumberGen(double seed){
    if (seed == 0) {
        m_igen = new boost::mt19937(static_cast<unsigned int>(time(NULL)));
        cout << "random seed is time!" << endl;
    } else {
        m_igen = new boost::mt19937(static_cast<unsigned int>(seed));
        cout << "random seed is " << seed << endl;
    }
}



void CConfiguration::moveBack(){
    //moves particle back to previous position
    for (int i = 0; i < 3; i++) {_ppos[i] = _prevpos[i];}
}
//****************************POTENTIALS**********************************************************//



void CConfiguration::calculateExpPotential(const double rSq, double& U, double& Fr){
    //function to calculate an exponential Potential U = U_0 * exp(-1 * r * k)
    // k is the interaction range. U_0 is the strength of the potential
    //which is attractive if direction = -1, and repulsive if direction = 1
    //The potential is weighted with kT!
    if (rSq < _cutoffExpSq && _potStrength != 0){
        const double r = sqrt(rSq);
        U = _potStrength * exp(-1 * r / _potRange);
        Fr = U / (_potRange * r);  //This is the force divided by the distance to the rod!
    }
    else{
        U=0;
        Fr=0;
    }
}


void CConfiguration::calculateExpHPI(const double r, double& U, double& Fr){
    cout << "NOTHING HERE" << endl;
    abort();
//	double u = _hpi_u * exp( - r / _hpi_k);
//	U += u;
//	Fr += u / (_hpi_k * r);
}


void CConfiguration::calculateExpPotentialMOD(const double r, double& U, double& Fr, int plane){
    cout << "NOT DEFINED!!!!!!!" << endl;
}

void CConfiguration::modifyPot(double& U, double& Fr, double dist){
    //function to modify the potential according to the distance along the polymer axis to the next neighbor,
    //in case the next neighboring polymer part is of opposite sign
    cout << "NOT DEFINED" << endl;
//     U = U * 4 * dist/_boxsize;
//     Fr = Fr * 4 * dist/_boxsize;
}

//****************************STERIC HINDRANCE****************************************************//

bool CConfiguration::testOverlap(){//TODO
    //Function to check, whether the diffusing particle of size psize is overlapping with any one of the rods (edges of the box)
    //most if borrowed from moveParticleAndWatch()
    bool overlaps = false;    
    double r_i = 0, r_k = 0;
    double r_abs = 0;
    
    if (_ranRod){
        int i,j;
        for (int axis=0;axis<3;axis++){
            for (int abc=0;abc<_rodarr[axis].size();abc++){
                for (int def=0;def<_rodarr[axis].size();def++){
                    i = axis +1;
                    if (i==3) i=0;
                    j=3-(i+axis);
                    if (testTracerOverlap(i, j, _rodarr[axis][abc][def].coord[i], _rodarr[axis][abc][def].coord[j])){
                        ifdebug(cout << "Overlap for\naxis" << axis << "\nabc " << abc << "\ndef " << def << endl;)
                        return true;
                    }
                }
            }
        }
    }
    else{
        for (int i = 0; i < 2; i++){
            for (int k = i+1; k < 3; k++){
                for (int ni = 0; ni < 2; ni++){
                    for (int nk = 0; nk < 2; nk++){
                        r_i = _ppos[i] - ni*_boxsize[i];
                        r_k = _ppos[k] - nk*_boxsize[k];       
                        r_abs = sqrt(r_i * r_i + r_k * r_k); //distance to the rods
                        if (r_abs < _pradius) overlaps = true;
                    }
                }
            }
        }
    }
    return overlaps;
}


void CConfiguration::calcLJPot(const double rSq, double& U, double& Fr){
    //Function to calculate the Lennard-Jones Potential
    double  por6 = pow((_pradius*_pradius / rSq ), 3);      //por6 stands for "p over r to the power of 6" . The 2 comes from the fact, that I need the particle radius, not the particle size
    ifdebug(if (4 * ( por6*por6 - por6 + 0.25 ) > 50) cout << "Very large LJ!!"<< endl;)
    U += 4 * ( por6*por6 - por6 + 0.25 );
    Fr +=  24 / ( rSq ) * ( 2 * por6*por6 - por6 );
}




//****************************OLD CODE****************************************************

//****************************POS HISTOGRAM****************************************************//

// void CConfiguration::initPosHisto(){
//     _posHistoM.resize(100);
//     for (int i = 0; i < 100; i++){
//         _posHistoM[i].resize(100);
//         for (int j = 0; j < 100; j++){
//             _posHistoM[i][j].resize(100, 0);  //initialize all the 100*100*100 matrix elements to zero!
//         }
//     }
// }
//
// void CConfiguration::addHistoValue(){
//     //adds a value to the position histogram
//     int x = _ppos[0] / _boxsize * 100;     //CAREFUL: THIS CAN'T BE DONE AT A POINT WHERE X MIGHT BE ZERO!!!
//     int y = _ppos[1] / _boxsize * 100;
//     int z = _ppos[2] / _boxsize * 100;
//     if ((x < 0) || (y < 0) || (z < 0) || (x > 99) || (y > 99) || (z > 99)){
//         cout << "The Position Histogram function 'conf.addHisto()' is in a bad place, since there is a negative position _ppos[]" << endl;
//         cout << "The current position is: " << _ppos[0] << " " << _ppos[1] << " " << _ppos[2] << endl;
//     }
//     _posHistoM[x][y][z] += 1;
// }
//
// void CConfiguration::printHistoMatrix(string folder){
//     //function to print the positionHistogram to a file called posHistoMatrix.txt
//     //The elements of the matrix are M[x][y][z]. First the z is counted from 1 to 100 in one row, then the y follows, then, after the 'X', the 100 x elements.
//
//     ofstream matrixfile;
//     matrixfile.open((folder + "/InstantValues/posHistoMatrix.txt").c_str());
//     int maxval = 0;
//
//     for (int i = 0; i < 100; i++){
//         for (int j = 0; j < 100; j++){
//             for (int k = 0; k < 100; k++){
//                 matrixfile << _posHistoM[i][j][k] << " ";
//                 if (maxval < _posHistoM[i][j][k] ) maxval = _posHistoM[i][j][k];
//             }
//             matrixfile << endl;
//         }
//         matrixfile << "X" << endl;
//     }
//     matrixfile << "MaxVal " << maxval;   //THis does not affect the grid files, since data is only copied to them when a "X" line comes!
//
//
//     matrixfile.close();
// }



void CConfiguration::calc_1RODONLY_MobilityForces(){   //see OldCode
}



void CConfiguration::calc_ZRODONLY_MobilityForces(){   //see OldCode
}



void CConfiguration::calc_YZRODONLY_MobilityForces(){   //see OldCode
}

double CConfiguration::getDisplacement(){
    double d = 0;
    for (int i = 0; i < 3; i++){
        d += pow((_timestep * _f_mob[i] + _mu_sto * _f_sto[i]), 2);
    }
    return sqrt(d);
}

/*
void CConfiguration::calcMobilityForces(){
    //calculate mobility forces from potential Epot
    double r_abs = 0;
    double r_i = 0, r_k = 0;
    double L1[] = {0, _boxsize, 0, _boxsize};  //these two arrays are only needed for the iteration.
    double L2[] = {0, 0, _boxsize, _boxsize};
    double utmp = 0, frtmp = 0;    //temporary "hilfsvariables"
    double Epot = 0;
    double z1, z2;
    if (_ranU){
        z1 = 1/4 * _boxsize;
        z2 = _boxsize - z1;   //z is in cylindrical coordinates. This indicates above/below which value the exp potential is modifed for random signs.
    }
    //reset mobility forces to zero
    for (int l = 0; l < 3; l++) {
        _f_mob[l] = 0;
    }
    const double r_c = _6root2 * _pradius;    //cutoff for Lennard-Jones calculation (at minimum)

    for (int i = 0; i < 3; i++){
            int k = i + 1;   //k always one direction "further", i.e. if i = 0 = x-direction, then k = 1 = y-direction
            if ( k == 3 ) k = 0;
            int plane = 3 - (i+k); //this is the current plane of the cylindrical coordinates
            for (int n = 0; n < 4; n++){
                r_i = _ppos[i] - L1[n];
                r_k = _ppos[k] - L2[n];
                //this is needed if we dont want the rods to cross each other to create a strong potential well
                if (plane == 0){
                    r_i -= _rodDistance;
                }
                else if (plane == 1){
                    r_k -= _rodDistance;
                    r_i -= _rodDistance;
                }
                r_abs = sqrt(r_i * r_i + r_k * r_k); //distance to the rods


                if (_potMod) calculateExpPotentialMOD(r_abs, utmp, frtmp, plane);
                else calculateExpPotential(r_abs, utmp, frtmp);


                if (_ranU){
                    utmp = utmp * _poly.get_sign(plane, n);
                    frtmp = frtmp * _poly.get_sign(plane, n);
                    if (_ppos[plane] > z2){
                        if (! _poly.samesign(1, plane, n)){
                            modifyPot(utmp, frtmp, _boxsize - _ppos[plane]);
                            _f_mob[plane] += utmp * 4 / _boxsize;              //this takes care of the derivative of the potential modification and resulting force
                        }
                    }
                    else if (_ppos[plane] < z1){
                        if (! _poly.samesign(-1, plane, n)){
                            modifyPot(utmp, frtmp, _ppos[plane]);
                            _f_mob[plane] -= utmp * 4 / _boxsize;              //this takes care of the derivative of the potential modification and resulting force
                        }
                    }
                }


                if (_LJPot && ( r_abs < r_c )) calcLJPot(r_abs, utmp, frtmp);


                Epot += utmp;
                _f_mob[i] += frtmp * r_i;
                _f_mob[k] += frtmp * r_k;
            }

    }
    _upot = Epot;
}


void CConfiguration::calcMobilityForces(){
    //calculate mobility forces from potential Epot
    double r_abs = 0;
    double r_i = 0, r_k = 0;
    double utmp = 0, frtmp = 0;     //temporary "hilfsvariables"
    double Epot = 0;
    //reset mobility forces to zero
    for (int l = 0; l < 3; l++) {
        _f_mob[l] = 0;
    }
    const double r_c = _6root2 * _pradius;    //cutoff for Lennard-Jones calculation (at minimum)

    for (int i = 0; i < 3; i++){
            int k = i + 1;   //k always one direction "further", i.e. if i = 0 = x-direction, then k = 1 = y-direction
            if ( k == 3 ) k = 0;
            int plane = 3 - (i+k); //this is the current plane of the cylindrical coordinates
            for (int ni = -1; ni < 3; ni++){
                for (int nk = -1; nk < 3; nk++){
                r_i = _ppos[i] - ni*_boxsize;
                r_k = _ppos[k] - nk*_boxsize;

                r_abs = sqrt(r_i * r_i + r_k * r_k); //distance to the rods


                if (_potMod) calculateExpPotentialMOD(r_abs, utmp, frtmp, plane);
                else calculateExpPotential(r_abs, utmp, frtmp);


                if (_LJPot && ( r_abs < r_c )) calcLJPot(r_abs, utmp, frtmp);


                Epot += utmp;
                _f_mob[i] += frtmp * r_i;
                _f_mob[k] += frtmp * r_k;
                }
            }

    }
    _upot = Epot;
}




void CConfiguration::countWallCrossing(int crossaxis, int exitmarker){

    if (_entryside[crossaxis] == exitmarker){                // case: entryside same as exitside
        _wallcrossings[0] += 1;
    }
    else if (_entryside[crossaxis] == (-1.0 * exitmarker)){  // case: entryside opposite exitside
        _wallcrossings[1] += 1;
    }
    else {                                                  // case: exiting "sideways", i.e. through one of the four other sides
        _wallcrossings[2] += 1;
    }
    for (int i = 0; i < 3; i++){                            // mark new entryside
        _entryside[i] = 0;
    }
    _entryside[crossaxis] = -1 * exitmarker;
}


*/

