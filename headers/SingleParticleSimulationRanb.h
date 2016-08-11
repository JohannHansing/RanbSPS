
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/filesystem.hpp>
#include "CAverage.h"
#include "CConfiguration.h"


using namespace std;


//Function declarations
template<typename T>
string toString(const T& value){
    ostringstream oss;
    oss << value;
    return oss.str();
}


template <typename T, size_t N>
inline
size_t sizeOfArray( const T(&)[ N ] )
{
  return N;
}


string createDataFolder(bool setPBC, string distribution, double timestep, double simtime, double potRange, double potStrength,
                        double particlesize, bool steric, bool ranRod, bool ranU, bool rand, double dvar, double polydiam, string tmp5){
    //NOTE: Maybe I can leave out dt, as soon as I settled on a timestep
    //NOTE: As soon as I create input-list with variables, I must change this function
    char range[5];
    sprintf(range, "%.3f", potRange);
    //In the definition of folder, the addition has to START WITH A STRING! for the compiler to know what to do (left to right).
    string folder = "sim_data";
    if (setPBC) folder += "/setPBC";
    if (ranU) folder = folder + "/ranU";
    if (ranRod) folder += "/ranRod";
    if (rand) folder += "/rand/d" + toString(dvar);
    folder += "/" + distribution;
    if (steric) folder = folder + "/steric";    //TODO steric2
    folder = folder
            + "/dt" + toString(timestep)
            + "/t" + toString(simtime)
            + "/a" + toString(polydiam)
            + "/p" + toString(particlesize)
            + "/k" + range
            + "/u" + toString(potStrength);
    boost::filesystem::create_directories(folder);
    boost::filesystem::create_directory(folder + "/InstantValues");
    boost::filesystem::create_directory(folder + "/Coordinates");
    return folder;
}


void settingsFile(bool setPBC, string folder, bool ranRod, double particlesize, double timestep, double runs, double steps, double potStrength, double potRange,
        bool rand, bool recordPosHisto, bool steric, bool ranU, string distribution, double dvar, double polydiam, string tmp5){
    //Creates a file where the simulation settings are stored
    //MAYBE ALSO INCLUDE TIME AND DATE!!
    ofstream settingsfile;
    settingsfile.open((folder + "/sim_Settings.txt").c_str());
    settingsfile << "Sim dir: " << folder << endl;
    settingsfile << "setPBC " << setPBC << endl;
    settingsfile << "Pore Distribution " << distribution << endl;
    settingsfile << "ranRod " << ranRod << endl;
    settingsfile << "TMP " << recordPosHisto << endl;//" (Bessel)" << endl;  //TODO Bessel!
    settingsfile << "rand " << rand << endl;
    settingsfile << "includesteric " << steric << endl;
    settingsfile << "ranU " << ranU  << endl;
    settingsfile << "p " << particlesize << endl;
    settingsfile << "dt " << timestep  << endl << "runs " << runs << endl << "steps " << steps << endl << "time: " << timestep*steps << endl;
    settingsfile << "k " << potRange << endl << "U_0 " << potStrength << endl;
    settingsfile << "dvar " << dvar << endl;
    settingsfile << "a " << polydiam << endl;

    settingsfile.close();
}
