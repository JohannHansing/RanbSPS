
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/filesystem.hpp>
#include "CAverage.h"
#include "CConfiguration.h"
#include "parameter_structs.h"

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
size_t sizeOfArray( const T(&)[ N ] ){
  return N;
}

string createDataFolder(paramstruct ps){
    //NOTE: Maybe I can leave out dt, as soon as I settled on a timestep
    //NOTE: As soon as I create input-list with variables, I must change this function
               //             cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ " <<tmp << endl;
    char range[12];
    sprintf(range, "%.3f", ps.urange);
    //In the definition of folder, the addition has to START WITH A STRING! for the compiler to know what to do (left to right).
    string folder = "sim_data";
    if (ps.setPBC) folder += "/setPBC";
    if (ps.Pointq) folder += "/pointq/drqop"+ toString(ps.drqop);
    if (ps.ranU) folder = folder + "/ranU";
    if (ps.mixU) folder += "/mixU";
    if (ps.mixU || ps.ranU) folder += "/uratio" + toString(ps.uratio) + "/Cratio" + toString(ps.Cratio);
    if (ps.ranRod) folder += "/ranRod";
    if (ps.rand) folder += "/rand/d" + toString(ps.dvar);
    folder += "/" + ps.distribution;
    if (ps.includeSteric) folder = folder + "/steric";    //TODO steric2
    folder = folder
            + "/dt" + toString(ps.timestep)
            + "/t" + toString(ps.simtime)
            + "/a" + toString(ps.polydiam)
            + "/p" + toString(ps.particlesize)
            + "/k" + range
            + "/u" + toString(ps.ustrength);
    boost::filesystem::create_directories(folder);
    boost::filesystem::create_directory(folder + "/InstantValues");
    boost::filesystem::create_directory(folder + "/Coordinates");
    return folder;
}


void settingsFile(string folder, paramstruct ps){
    //Creates a file where the simulation settings are stored
    //MAYBE ALSO INCLUDE TIME AND DATE!!
    ofstream settingsfile;
    settingsfile.open((folder + "/sim_Settings.txt").c_str());
    settingsfile << "Sim dir: " << folder << endl;
    settingsfile << "setPBC " << ps.setPBC << endl;
    settingsfile << "Pointq " << ps.Pointq << endl;
    settingsfile << "mixU " << ps.mixU << endl;
    settingsfile << "Pore Distribution " << ps.distribution << endl;
    settingsfile << "ranRod " << ps.ranRod << endl;
    settingsfile << "rand " << ps.rand << endl;
    settingsfile << "includesteric " << ps.includeSteric << endl;
    settingsfile << "ranU " << ps.ranU  << endl;
    settingsfile << "p " << ps.particlesize << endl;
    settingsfile << "dt " << ps.timestep  << endl << "runs " << ps.runs << endl << "steps " << ps.steps << endl << "time: " << ps.timestep*ps.steps << endl;
    settingsfile << "k " << ps.urange << endl << "U_0 " << ps.ustrength << endl;
    settingsfile << "dvar " << ps.dvar << endl;
    settingsfile << "a " << ps.polydiam << endl;
    settingsfile << "drqop " << ps.drqop << endl;
    if (ps.mixU || ps.ranU){
        settingsfile << "uratio " << ps.uratio << endl;
        settingsfile << "Cratio " << ps.Cratio << endl;
    }

    settingsfile.close();
}
