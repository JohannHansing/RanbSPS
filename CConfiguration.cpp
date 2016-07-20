#include "headers/CConfiguration.h"


using namespace std;
using namespace Eigen;


CConfiguration::CConfiguration(){
}

CConfiguration::CConfiguration(
        string distribution, double timestep,  double potRange,  double potStrength, const bool rand,
        double psize, const bool posHisto, const bool steric, const bool ranU, bool ranRod, double dvar, double polydiam, string peptide){
    setRanNumberGen(0);
    _potRange = potRange;
    _potStrength = potStrength;
    _pradius = psize/2.;
    //TODO overlap
    _polyrad = polydiam/2.;//TODO test
    _polydiamSq = polydiam*polydiam;
    _cylLJSq = pow(_pradius + _polyrad,2);
    _cutoffExpSq = pow(6*_potRange + _pradius + _polyrad,2);
    _timestep = timestep;
    _ranRod = ranRod;
    _LJPot = (steric == false) && (psize != 0);
    _ranU = ranU;
    _rand = rand;
    _dvar = dvar;
    _upot = 0;
    _mu_sto = sqrt( 2 * _timestep );                 //timestep for stochastic force
    for (int i = 0; i < 3; i++){
        _ppos(i) = _bdef/2.;
        _startpos[i] = _ppos(i);
        _entryside[i] = 0;
        _wallcrossings[i] = 0;
        _boxCoord[i] = 0;
        _prevpos(i) = _ppos(i);
    }
    if (distribution=="gamma"){
       cout << "_alpha = " << _alpha << endl;
    }
    else if (distribution == "gamma2"){
        _alpha = 5.; _beta = 2.;
       cout << "_alpha = " << _alpha << endl;
    }
    else if (distribution == "gamma4"){
        _alpha = 2.5; _beta = 4.;
       cout << "_alpha = " << _alpha << endl;
    }
    else if (distribution == "fixb" ){
       _fixb = true;
       cout << "-----> b fixed !" << endl;
    }
    else {
        cout << "Invalid distribution!\nAborting!!!!!!!!!!" << endl;
        abort();
    }

    //TODO del
    //testgamma();
    initRanb();

    if (_ranRod){
        //TODO rel
        initRodsArr();
        //initRodsRel();
    }
    if (_rand){
        cout << "init Rand..." << endl;
        initRand();
    }
}

void CConfiguration::updateStartpos(){
    //This function is used if the particle should keep moving after each run, and not start at _resetpos again, like in the first run
    //This (hopefully) will give better averages without having to spend a lot of steps in the beginning of each run to get away from _resetpos
    for (int i = 0; i < 3; i++){
    _startpos[i] = _ppos(i) + _boxCoord[i];
    }
}



void CConfiguration::makeStep(){
    //move the particle according to the forces and record trajectory like watched by outsider
    for (int i = 0; i < 3; i++){
        _prevpos(i) = _ppos(i);
        const double disp = _timestep * _f_mob(i) + _mu_sto * _f_sto(i);
        _ppos(i) += disp;
        if (abs(disp) > 5){
            cout << "**** Way too big jump! " << endl;
        }
        else if (std::isnan(_ppos(i))){
            cout << "axis " << i << "\n_prevpos " << _prevpos(i) << "\n_ppos " << _ppos(i) << endl;
            cout << "_f_mob(i) " << _f_mob(i) << "\n_f_sto(i) " << _f_sto(i) << endl;
            abort();
        }
    }
}



void CConfiguration::checkBoxCrossing(){
    //should the particle cross the confinement of the cube, let it appear on the opposite side of the box
    int exitmarker = 0;
    for (int i = 0; i < 3; i++){
        exitmarker =0;
        if (_ppos(i) < -0.05*_b_array[i][0]){//Only create new rod config if the particle has crossed the border by a certain fraction
            _ppos(i) += _b_array[i][0];
            _boxCoord[i] -= _b_array[i][0];
            exitmarker = -1;
        }
        else if (_ppos(i) > 1.05*_b_array[i][2]){//Only create new rod config if the particle has crossed the border by a certain fraction
            _ppos(i) -= _boxsize[i];
            _boxCoord[i] += _boxsize[i];
            exitmarker = 1;
            if (_ppos(i) > 10){
                cout << "Bad _ppos after boxcrossing. Aborting!" << endl;
                cout << "_ppos b4 " << _ppos(i) + _boxsize[i] << endl;
                cout << "_ppos after " << _ppos(i) << endl;
                cout << "\naxis " << i << "\n_boxsize[ax] " << _boxsize[i] << endl;
                abort();
            }
        }
        if (exitmarker!=0){
            updateRanb(i,exitmarker);
            if (_rand){
                updateRand(i,exitmarker);
                ifdebug(
                    cout << "[["<<i<<"," << exitmarker <<"] ";
                    prinRodPos(0); // cout print rod pos!
                )
            }
            if (_ranRod){
                //TODO rel
                updateRodsArr(i, exitmarker);
                //updateRodsRel(i, exitmarker);
            }
            if (_ranU){
                for (auto & vrods : _drods[i]){
                    for (auto & rod :  vrods){
                        rod.shiftSigns(exitmarker);
                    }
                }
            }
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
    // samples from normal distribution with standard deviation 1 (later sqrt(2) is multiplied)
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > ran_gen(
            *m_igen, boost::normal_distribution<double>(0, 1));

    for (int i = 0; i < 3; i++) {
        _f_sto(i) = ran_gen();
    }
}


void CConfiguration::calcMobilityForces(){
    //calculate mobility forces from potential Epot - Unified version that includes 2ndOrder if k is larger than or equal 0.2 b , except if ranPot is activated.
    double Epot = 0;
    double rcSq = 1.25992 * _cylLJSq;
    double r_i = 0, r_k = 0;
    array<double,4> r_is, r_ks;
    std::array<double, 16> ri_arr, rk_arr, rSq_arr;
    double r_absSq;
    double utmp = 0, frtmp = 0;     //temporary "hilfsvariables"
    //reset mobility forces to zero
    _f_mob = Vector3d::Zero();


    for (int i = 0; i < 3; i++){
        unsigned int cnt=0;// counter to loop over array indices
        int k = i + 1;   //k always one direction "further", i.e. if i = 0 = x-direction, then k = 1 = y-direction
        if ( k == 3 ) k = 0;
        int plane = 3 - (i+k); //this is the current plane of the cylindrical coordinates
        int n = 0;     // reset counter for index of next rod in plane  n = 0, 1, 2, 3 -> only needed for ranPot
        double z1, z2, z1inv;
        if (_ranU){
            z1 = 0.25 * _boxsize[plane];
            z2 = _boxsize[plane] - z1;   //z is in cylindrical coordinates. This indicates above/below which value the exp potential is modifed for random signs.
            z1inv = 1./z1;
        }

        if (_ranRod){
            for (int abc=0;abc<3;abc++){
                for (int def=0;def<3;def++){
                    r_i = _ppos(i) - _rodarr[plane][abc][def].coord(i);
                    r_k = _ppos(k) - _rodarr[plane][abc][def].coord(k);
                    ri_arr[cnt]=(r_i);
                    rk_arr[cnt]=(r_k);
                    rSq_arr[cnt]=( r_i * r_i + r_k * r_k);
                    cnt++;
                }
            }
        }
        else if (_rand){
            for (int abcd=0;abcd<4;abcd++){
                for (int efgh=0;efgh<4;efgh++){
                    r_i = _ppos(i) - _drods[plane][abcd][efgh].coord(i);
                    r_k = _ppos(k) - _drods[plane][abcd][efgh].coord(k);
                    ri_arr[cnt]=(r_i);
                    rk_arr[cnt]=(r_k);
                    rSq_arr[cnt]=( r_i * r_i + r_k * r_k);
                    cnt++;
                }
            }
        }
        else{
            //This creates the distance vectors from the rods to the tracer
            r_ks[0] = _ppos(k) + _b_array[k][0];
            r_is[0] = _ppos(i) + _b_array[i][0];
            for (int rodi = 0; rodi < 3; rodi++){
                r_ks[rodi+1] = r_ks[rodi] - _b_array[k][rodi];
                r_is[rodi+1] = r_is[rodi] - _b_array[i][rodi];
            }
            for (int nk = 0; nk < r_ks.size(); nk++){
                for (int ni = 0; ni < r_is.size(); ni++){
                    r_i = r_is[ni];
                    r_k = r_ks[nk];
                    ri_arr[cnt]=(r_i);
                    rk_arr[cnt]=(r_k);
                    rSq_arr[cnt]=( r_i * r_i + r_k * r_k);
                    cnt++;
                }
            }
        }
    
        //TODO distarr[i] = rSq_arr; //store distances for writing them to a file later
        for (int j=0;j<cnt;j++){
            const double rSq = rSq_arr.at(j);
            calculateExpPotential(rSq, utmp, frtmp);

        
            if (_ranU){
                int abcd = j/4; // This could be made more efficient by replacing the j loop with two loops abcd efgh and a counter j
                int efgh = j%4;
                int sign = _drods[plane][abcd][efgh].signs[1];
                //cout << "abcd " << abcd << "efgh " << efgh << " sign: " << sign << endl;
                utmp *= sign;
                frtmp *= sign;
                if (_ppos(plane) > z2){
                    if (! _drods[plane][abcd][efgh].samesign[1]){
                        _f_mob(plane) += utmp * z1inv;              //this takes care of the derivative of the potential modification and resulting force
                        modifyPot(utmp, frtmp, (_boxsize[plane] - _ppos(plane)) * z1inv);
                    }
                }
                else if (_ppos(plane) < z1){
                    if (! _drods[plane][abcd][efgh].samesign[0]){
                        _f_mob(plane) -= utmp * z1inv;              //this takes care of the derivative of the potential modification and resulting force
                        modifyPot(utmp, frtmp, _ppos(plane) * z1inv);
                    }
                }
            }

            if (_LJPot && ( rSq < rcSq )) calcLJPot(rSq, utmp, frtmp, _cylLJSq);

            //TODO del
            if (utmp > 100){
                cout << "utmp " << utmp  <<"\nr " << sqrt(rSq) << "\nindex j "  << j << "\nplane " << plane << endl;
                cout << "ri " << ri_arr[j] << "\nrk " << rk_arr[j] << endl;
            }


            Epot += utmp;
            _f_mob(i) += frtmp * ri_arr[j];
            _f_mob(k) += frtmp * rk_arr[j];
        }
    }
    _upot = Epot;
}



void CConfiguration::saveXYZTraj(string name, int move, string flag) {
    // so far this only writes the positions of the tracer
    Eigen::Map<Eigen::Vector3d> boxCoordinates(_boxCoord);
    Vector3d rtmp;
    if(flag=="w") {    //write to new file
        /*if(m_traj_file!=NULL) {
            fclose(m_traj_file);
        }*/
        m_traj_file = fopen(name.c_str(), flag.c_str());
        if(m_traj_file==NULL) {
            cout << "error creating trajfile" << endl;
        }
    }

    fprintf(m_traj_file, "%d\n%s (%8.3f %8.3f %8.3f) t=%u \n", 1, "sim_name", _boxsize[0], _boxsize[1], _boxsize[2], move);


    rtmp = _ppos;//+boxCoordinates;
    fprintf(m_traj_file, "%3s%9.3f%9.3f%9.3f \n","H", rtmp(0), rtmp(1),  rtmp(2));


    //fflush(m_traj_file);

    if(flag=="c") {    //close file
        if(m_traj_file!=NULL) { fclose(m_traj_file); }
    }
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
    _ppos = _prevpos;
}
//****************************POTENTIALS**********************************************************//



void CConfiguration::calculateExpPotential(const double rSq, double& U, double& Fr){
    //function to calculate an exponential Potential U = U_0 * exp(-1 * r * k)
    // k is the interaction range. U_0 is the strength of the potential
    //which is attractive if direction = -1, and repulsive if direction = 1
    //The potential is weighted with kT!
    if (rSq < _cutoffExpSq && _potStrength != 0){
        const double r = sqrt(rSq);
        U = _potStrength * exp(-1.0 * r / _potRange);
        Fr = U / (_potRange * r);  //This is the force divided by the distance to the rod!
    }
    else{
        U=0;
        Fr=0;
    }
}


void CConfiguration::calculateExpHPI(const double r, double& U, double& Fr){
    cout << "\n\n....calculateExpHPI: NOTHING HERE ..." << endl;
    abort();
//	double u = _hpi_u * exp( - r / _hpi_k);
//	U += u;
//	Fr += u / (_hpi_k * r);
}


void CConfiguration::calculateExpPotentialMOD(const double r, double& U, double& Fr, int plane){
    cout << "NOT DEFINED!!!!!!!" << endl;
}

void CConfiguration::modifyPot(double& U, double& Fr, double weight){
    //function to modify the potential according to the distance along the polymer axis to the next neighbor,
    //in case the next neighboring polymer part is of opposite sign
    //the weight is determined by the distance to the point where the sign changes, divided by the boxsize: weight = 4 * dist / boxsize, such that it's 1 for dist=boxsize/4
//    cout << "NOT DEFINED" << endl;
    U *= weight;
    Fr *= weight;
}

//****************************STERIC HINDRANCE****************************************************//

bool CConfiguration::testOverlap(){//TODO relRod
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
                    if (testTracerOverlap(i, j, _rodarr[axis][abc][def].coord(i), _rodarr[axis][abc][def].coord(j))){
                        ifdebug(cout << "Overlap for\naxis" << axis << "\nabc " << abc << "\ndef " << def << endl;)
                        return true;
                    }
                }
            }
        }
    }
    else if (_rand){
        int i,j;
        for (int axis=0;axis<3;axis++){
            for (int abc=0;abc<_drods[axis].size();abc++){
                for (int def=0;def<_drods[axis].size();def++){
                    i = axis +1;
                    if (i==3) i=0;
                    j=3-(i+axis);
                    if (testTracerOverlap(i, j, _drods[axis][abc][def].coord(i), _drods[axis][abc][def].coord(j))){
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
                        r_i = _ppos(i) - ni*_boxsize[i];
                        r_k = _ppos(k) - nk*_boxsize[k];
                        r_abs = sqrt(r_i * r_i + r_k * r_k); //distance to the rods
                        if (r_abs < _pradius) overlaps = true;
                    }
                }
            }
        }
    }
    return overlaps;
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
//     int x = _ppos(0) / _boxsize * 100;     //CAREFUL: THIS CAN'T BE DONE AT A POINT WHERE X MIGHT BE ZERO!!!
//     int y = _ppos(1) / _boxsize * 100;
//     int z = _ppos(2) / _boxsize * 100;
//     if ((x < 0) || (y < 0) || (z < 0) || (x > 99) || (y > 99) || (z > 99)){
//         cout << "The Position Histogram function 'conf.addHisto()' is in a bad place, since there is a negative position _ppos()" << endl;
//         cout << "The current position is: " << _ppos(0) << " " << _ppos(1) << " " << _ppos(2) << endl;
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
        d += pow((_timestep * _f_mob(i) + _mu_sto * _f_sto(i)), 2);
    }
    return sqrt(d);
}

