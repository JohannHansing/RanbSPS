
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
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



string createDataFolder(bool ranRod, double timestep, double simtime, double potRange, double potStrength,
        double boxsize, double particlesize, double rDist, bool potMod, bool steric, bool randomPot, bool hpi, double hpi_u, double hpi_k, double n_rods){
    //NOTE: Maybe I can leave out dt, as soon as I settled on a timestep
    //NOTE: As soon as I create input-list with variables, I must change this function
    char range[5];
    sprintf(range, "%.3f", potRange);
    //In the definition of folder, the addition has to START WITH A STRING! for the compiler to know what to do (left to right).
    string folder = "sim_data/noreset";
    if (ranRod) folder = folder + "/ranRod/nrods" + toString(n_rods);
    if (randomPot) folder = folder + "/ranPot";
    if (steric) folder = folder + "/steric";    //TODO steric2
    if (potMod) folder = folder +  "/potMod";   //"/potMod";  TODO!!! Bessel
    if (hpi) folder = folder + "/HPI/hpiu" + toString(hpi_u) + "/hpik" + toString(hpi_k);
    folder = folder
            + "/dt" + toString(timestep)
            + "/t" + toString(simtime)
            + "/d" + toString(rDist)
            + "/b" + toString(boxsize)
            + "/p" + toString(particlesize)
            + "/k" + range
            + "/u" + toString(potStrength);
    boost::filesystem::create_directories(folder);
    boost::filesystem::create_directory(folder + "/InstantValues");
    boost::filesystem::create_directory(folder + "/Coordinates");
    return folder;
}


void settingsFile(string folder, bool ranRod, double particlesize, double boxsize, double timestep, double runs, double steps, double potStrength, double potRange, double rDist,
        bool potMod, bool recordMFP, bool steric, bool randomPot, bool hpi, double hpi_u, double hpi_k, double n_rods){
    //Creates a file where the simulation settings are stored
    //MAYBE ALSO INCLUDE TIME AND DATE!!
    ofstream settingsfile;
    settingsfile.open((folder + "/sim_Settings.txt").c_str());
    settingsfile << "Sim dir: " << folder << endl;
    settingsfile << "ranRod " << ranRod << endl;
    settingsfile << "n_rods " << n_rods << endl;
    settingsfile << "potMod " << potMod << endl;//" (Bessel)" << endl;  //TODO Bessel!
    settingsfile << "recordMFP " << recordMFP << endl;
    settingsfile << "includesteric " << steric << endl;
	settingsfile << "ranPot " << randomPot  << endl;
	settingsfile << "HPI " << hpi  << endl;
	if (hpi == true){
		settingsfile << "hpi_u " << randomPot  << endl;
		settingsfile << "hpi_ k " << randomPot  << endl;
    }
    settingsfile << "d " << rDist << endl;
    settingsfile << "p " << particlesize << endl;
	settingsfile << "b " << boxsize << endl;
    settingsfile << "dt " << timestep  << endl << "runs " << runs << endl << "steps " << steps << endl << "time: " << timestep*steps << endl;
    settingsfile << "k " << potRange << endl << "U_0 " << potStrength << endl;
    
    settingsfile.close();
}


