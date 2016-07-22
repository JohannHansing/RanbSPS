#include "headers/SingleParticleSimulationRanb.h"

using namespace std;



int main(int argc, const char* argv[]){
    //Main includes the iteration loops for the simulation

    //NOTE: so far wallcrossings is added for all runs!!! makes it kind of incorrect, since after each run, ppos is reset.
    //NOTE: so far saving Instant Values for each tenth step!

    //TRIGGERS:
    string distribution = argv[1];    // TODO
    bool ranRod = (strcmp(argv[2] , "true") == 0 ) ;
    bool rand = (strcmp(argv[3] , "true") == 0 ) ; 
    bool writeRods = (strcmp(argv[4] , "true") == 0 ) ;  // TODO CHANGE THIS TO PBC or something
    bool recordPosHisto = (strcmp(argv[5] , "true") == 0 ) ;
    bool includeSteric = (strcmp(argv[6] , "true") == 0 ) ;  // steric 2
    bool ranU = (strcmp(argv[7] , "true") == 0 ) ;
    string tmp5 = argv[8];
    int boolpar = 8;
    ifdebug(cout << "copied bools. ";)

    // Checking for correct structure of input arguments
    for (int k= 1; k < argc; k++ ) cout << "parameter " << k << " " << argv[k] << endl;
    for (int b_i=2; b_i<boolpar; b_i++){
        if (!((strcmp(argv[b_i] , "true") == 0 )  || (strcmp(argv[b_i] , "false") == 0 ))){
            cerr << "Error; Bool parameter " << b_i << " is not either 'true' or 'false'!" << endl;
            exit(1);
        }
    }

    int runs = atoi( argv[boolpar+1] );                       // Number of Simulation runs to get mean values from
    double timestep = atof( argv[boolpar+2] );
    int simtime = atoi( argv[boolpar+3] );                   // simulation time
    int instantvalues = 200;
    unsigned int steps;

    double particlesize = atof( argv[boolpar+4] );
    double urange = atof( argv[boolpar+5] );
    double ustrength = atof( argv[boolpar+6] );
    double dvar = atof( argv[boolpar+7] );
    double polydiam = atof( argv[boolpar+8] );
    unsigned int saveInt;
    int instValIndex;                             //Counter for addInstantValue

    steps = simtime/timestep;
    saveInt = steps/instantvalues;
    const int trajout = (int)(10/timestep);
    

    ifdebug(cout << "copied  params. ";)

    cout << "distribution " << distribution << endl;

    if (distribution != "fixb" && rand){
        cout << "b needs to be fixed at this time to include rand!" << endl;
        abort();
    }
    if (ranRod  && rand){
        cout << "Cant have both rand and ranRod!" << endl;
        abort();
    }
    if (ranRod  && ranU){
        cout << "ranRod + ranU will not work properly due to definition of calcMobilityForces for ranU!" << endl;
        abort();
    }

    //Create data folders and print location as string to string "folder"
    string folder = createDataFolder(distribution, timestep, simtime, urange, ustrength, particlesize, includeSteric, ranRod, ranU, rand, 
                             dvar,polydiam, tmp5);
    ifdebug(cout << "created folder. ";)
    cout << "writing to folder " << folder << endl;


    //initialize averages
    CAverage energyU = CAverage("Upot", folder, instantvalues, runs);
    CAverage squareDisp = CAverage("squaredisp", folder, instantvalues, runs);
    ifdebug(cout << "created CAverage files. ";)

    //initialize instance of configuration
    CConfiguration conf = CConfiguration(writeRods, distribution,timestep, urange, ustrength, rand, particlesize, recordPosHisto, 
                            includeSteric, ranU, ranRod, dvar,polydiam, tmp5);
    ifdebug(cout << "created CConf conf. ";)



    unsigned int stepcount = 0;
    ofstream trajectoryfile;
    trajectoryfile.open((folder + "/Coordinates/trajectory.txt").c_str());
    
    ofstream distancesfile;
    // TODO distancefile
    distancesfile.open((folder + "/Coordinates/squareDistances.txt").c_str());
    
    settingsFile(folder, ranRod, particlesize, timestep, runs, steps, ustrength, urange, rand, recordPosHisto, includeSteric, ranU, distribution, 
                    dvar, polydiam, tmp5);
                    
    //create .xyz file to save the trajectory for VMD
    string traj_file = folder + "/Coordinates/single_traj.xyz";
    conf.saveXYZTraj(traj_file,0,"w");
    
    // write rod positions of rods are fixed throughout the simulation
    ofstream rodPosFile;
    FILE*  snapFile;
    if (writeRods){
        rodPosFile.open((folder + "/Coordinates/rodposfile.txt").c_str());
        rodPosFile << "set rad 0.1\nmol new" << endl;
        conf.writeRodForVMD(rodPosFile);
        snapFile = fopen((folder + "/Coordinates/snapFile.xyz").c_str(), "w");
        fprintf(snapFile, "%s\n%s (%8.3f %8.3f %8.3f) t=%u \n", "XXX", "sim_name", 10., 10., 10., 0);
    }
    


    //cout << "Starting Run Number: " << simcounter << " out of " << totalsims << endl;
    cout << "Starting Simulation!" << endl;


// **************START OF RUNS-LOOP
    for (int l = 0; l<runs; l++){

        conf.updateStartpos();

        instValIndex = 0;


        for (int i = 0; i < steps; i++){  //calculate stochastic force first, then mobility force!!


            conf.calcStochasticForces();

            conf.calcMobilityForces();


            if (((i+1)%saveInt) == 0){       //saving Instant Values for each saveInt'th step!
                energyU.addInstantValue(conf.getUpot(), instValIndex);
                instValIndex += 1;
            }


            conf.makeStep();    //move particle at the end of iteration

            /*
            if (includeSteric && conf.testOverlap()) conf.moveBack();   //TODO steric2
            else conf.checkBoxCrossing();
            */



                //TODO steric
            while (includeSteric && conf.testOverlap()){
                conf.moveBack();
                conf.calcStochasticForces();
                conf.makeStep();
            }
            if (conf.checkBoxCrossing()){//check if particle has crossed the confinement of the box
                // if (writeRods){
//                     rodPosFile << "#" << endl;
//                     conf.writeRodForVMD(rodPosFile);
//                 }
            }

            stepcount++;
            if (stepcount%trajout == 0) {
                std::vector<double> ppos = conf.getppos();
                trajectoryfile << fixed << stepcount * timestep << "\t" << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;
                ifdebug(cout << stepcount * timestep << "\t" << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;)
                //TODO pass distancefile to function in conf.
                if (stepcount%(100*trajout) == 0) conf.writeDistances( distancesfile, stepcount);
            }
            
            
            
            if (((i+1)%100 == 0) && (l == 0)){       //Save the first trajectory to file
                conf.saveXYZTraj(traj_file, i, "a");                    // TODO change back ((i+1)%XXX == 0) to 100

                if (((i+1)%100000 == 0) && writeRods){
                    std::vector<double> ppos = conf.getppos_rel();
                    fprintf(snapFile, "%3s%9.3f%9.3f%9.3f \n","H", ppos[0], ppos[1],  ppos[2]);
                }
            }
        }
        if (l==0){
            conf.saveXYZTraj(traj_file, steps, "c"); // Close XYZ traj_file
        }
        
        if (l%100 == 0){
            cout << "run " << l << endl;
            if (l==1000) energyU.saveAverageInstantValues(saveInt*timestep);
        }


    }//----------END OF RUNS-LOOP



    //watch out: Save average instant values at timeInterval: timestep * saveinterval saveInt!!!
    squareDisp.saveAverageInstantValues(saveInt*timestep);




    //printReport(ranRod, conf.getwallcrossings(0), conf.getwallcrossings(1), conf.getwallcrossings(2), timestep, urange, ustrength, rodDist, particlesize, runs,
    //        sizeOfArray(timestep), sizeOfArray(urange), sizeOfArray(ustrength), sizeOfArray(rodDist), sizeOfArray(particlesize), potentialMod, includeSteric);


    cout << "Simulation Finished" << endl;

    fclose(snapFile);
    trajectoryfile.close();
    rodPosFile.close();
    // TODO distancefile
    distancesfile.close();

    return 0;
}
