#include "headers/SingleParticleSimulationRanb.h"

using namespace std;



int main(int argc, const char* argv[]){
    //Main includes the iteration loops for the simulation

    //NOTE: so far wallcrossings is added for all runs!!! makes it kind of incorrect, since after each run, ppos is reset.
    //NOTE: so far saving Instant Values for each tenth step!
    
    struct paramstruct ps;

    //TRIGGERS:
    ps.distribution = argv[1];    // TODO
    ps.ranRod = (strcmp(argv[2] , "-X-") == 0 ) ;
    ps.rand = (strcmp(argv[3] , "-X-") == 0 ) ; 
    ps.setPBC = (strcmp(argv[4] , "-X-") == 0 ) ;  // TODO CHANGE THIS TO PBC or something
    ps.Tmp4 = (strcmp(argv[5] , "-X-") == 0 ) ;
    ps.includeSteric = (strcmp(argv[6] , "-X-") == 0 ) ;  // steric 2
    ps.ranU = (strcmp(argv[7] , "-X-") == 0 ) ;
    ps.Pointq = (strcmp(argv[8] , "-X-") == 0 ) ;
    int boolpar = 8;
    ifdebug(cout << "copied bools. ";)

    // Checking for correct structure of input arguments
    for (int k= 1; k < argc; k++ ) cout << "parameter " << k << " " << argv[k] << endl;
    for (int b_i=2; b_i<boolpar; b_i++){
        if (!((strcmp(argv[b_i] , "-X-") == 0 )  || (strcmp(argv[b_i] , "---") == 0 ))){
            cerr << "Error; Bool parameter " << b_i << " is not either '-X-' or '---'!" << endl;
            exit(1);
        }
    }

    ps.runs = atoi( argv[boolpar+1] );                       // Number of Simulation runs to get mean values from
    ps.timestep = atof( argv[boolpar+2] );
    ps.simtime = atoi( argv[boolpar+3] );                   // simulation time
    ps.instantvalues = 200;
    ps.steps = ps.simtime/ps.timestep;

    ps.particlesize = atof( argv[boolpar+4] );
    ps.urange = atof( argv[boolpar+5] );
    ps.ustrength = atof( argv[boolpar+6] );
    ps.dvar = atof( argv[boolpar+7] );
    ps.polydiam = atof( argv[boolpar+8] );
    ps.drqop = atof( argv[boolpar+9] );

    // main loop parameters
    unsigned int saveInt= ps.steps/ps.instantvalues;
    const int trajout = (int)(10/ps.timestep);
    int instValIndex;                             //Counter for addInstantValue



    ifdebug(cout << "copied  params. ";)

    cout << "distribution " << ps.distribution << endl;

    if (ps.distribution != "fixb" && ps.rand){
        cout << "b needs to be fixed at this time to include rand!" << endl;
        abort();
    }
    if (ps.ranRod  && ps.rand){
        cout << "Cant have both rand and ranRod!" << endl;
        abort();
    }
    if ((ps.ranRod  && (ps.ranU || ps.Pointq)) || (ps.setPBC && ps.Pointq) || (ps.ranRod && ps.Pointq)){
        cout << "ranRod + ranU will not work properly due to definition of calcMobilityForces for ranU!" << endl;
        cout << "Or some other bad combination of bool parameters!" << endl;
        abort();
    }

    //Create data folders and print location as string to string "folder"
    string folder = createDataFolder( ps);
    ifdebug(cout << "created folder. ";)
    cout << "writing to folder " << folder << endl;

    //initialize averages
    CAverage energyU = CAverage("Upot", folder, ps.instantvalues, ps.runs);
    CAverage squareDisp = CAverage("squaredisp", folder, ps.instantvalues, ps.runs);
    ifdebug(cout << "created CAverage files. ";)

    //initialize instance of configuration
    CConfiguration conf = CConfiguration(ps);
    ifdebug(cout << "created CConf conf. ";)



    unsigned long long stepcount = 0;
    ofstream trajectoryfile;
    trajectoryfile.open((folder + "/Coordinates/trajectory.txt").c_str());
    
    ofstream distancesfile;
    // TODO distancefile
    distancesfile.open((folder + "/Coordinates/squareDistances.txt").c_str());
    
    settingsFile(folder, ps);
                    
    //create .xyz file to save the trajectory for VMD
    //string traj_file = folder + "/Coordinates/single_traj.xyz";
    //conf.saveXYZTraj(traj_file,0,"w");
    
    // write rod positions of rods are fixed throughout the simulation
    ofstream rodPosFileVMD;
    FILE*  snapFile;
    if (ps.setPBC){
        rodPosFileVMD.open((folder + "/Coordinates/rodPosFileVMD.txt").c_str());
        rodPosFileVMD << "set rad 0.15\nmol new" << endl;
        if (ps.ranU) conf.writeRodForVMDranU(rodPosFileVMD);
        else conf.writeRodForVMD(rodPosFileVMD);
        snapFile = fopen((folder + "/Coordinates/snapFile.xyz").c_str(), "w");
        fprintf(snapFile, "%s\n%s (%8.3f %8.3f %8.3f) t=%u \n", "XXX", "sim_name", 10., 10., 10., 0);
    }
    
    ofstream rodsfile;
    rodsfile.open((folder + "/Coordinates/rodsfile.txt").c_str());


    //cout << "Starting Run Number: " << simcounter << " out of " << totalsims << endl;
    cout << "Starting Simulation!" << endl;


// **************START OF RUNS-LOOP
    for (int l = 0; l<ps.runs; l++){

        conf.updateStartpos();

        instValIndex = 0;


        for (int i = 0; i < ps.steps; i++){  //calculate stochastic force first, then mobility force!!


            conf.calcStochasticForces();

            conf.calcMobilityForces();


            if (((i+1)%saveInt) == 0){       //saving Instant Values for each saveInt'th step!
                energyU.addInstantValue(conf.getUpot(), instValIndex);
                instValIndex += 1;
            }


            conf.makeStep();    //move particle at the end of iteration

            /*
            if (ps.includeSteric && conf.testOverlap()) conf.moveBack();   //TODO steric2
            else conf.checkBoxCrossing();
            */



                //TODO steric
            while (ps.includeSteric && conf.testOverlap()){
                conf.moveBack();
                conf.calcStochasticForces();
                conf.makeStep();
            }
            if (conf.checkBoxCrossing()){//check if particle has crossed the confinement of the box
                // if (ps.setPBC){
//                     rodPosFileVMD << "#" << endl;
//                     conf.writeRodForVMD(rodPosFileVMD);
//                 }
            }

            stepcount++;
            if (stepcount%trajout == 0) {
                std::vector<double> ppos = conf.getppos();
                trajectoryfile << fixed << stepcount * ps.timestep << "\t" << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;
                ifdebug(cout << stepcount * ps.timestep << "\t" << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;)
                //TODO pass distancefile to function in conf.
                if (stepcount%(10*trajout) == 0){
                    conf.writeDistances( distancesfile, stepcount);
                }
                if (ps.setPBC && (stepcount%(10*trajout) == 0)){ //should be 10
                    std::vector<double> ppos = conf.getppos_rel();
                    fprintf(snapFile, "%3s%9.3f%9.3f%9.3f \n","H", ppos[0], ppos[1],  ppos[2]);
                }
            }
            
            
            
            if ( (l == 0) && ((i+1)%100 == 0)){       //Save the first trajectory to file
                //conf.saveXYZTraj(traj_file, i, "a");                    // TODO change back ((i+1)%XXX == 0) to 100
            }
        }
        if (l==0){
            //conf.saveXYZTraj(traj_file, ps.steps, "c"); // Close XYZ traj_file
        }
        
        if (l%50 == 0){
            cout << "run " << l << endl;
            if (ps.setPBC) fflush(snapFile); // flush the snapFile data to file
            energyU.saveAverageInstantValues(saveInt*ps.timestep);
            conf.writeRodsToFile(rodsfile, stepcount*ps.timestep);
        }
    }//----------END OF RUNS-LOOP



    //watch out: Save average instant values at timeInterval: ps.timestep * saveinterval saveInt!!!
    squareDisp.saveAverageInstantValues(saveInt*ps.timestep);




    //printReport(ranRod, conf.getwallcrossings(0), conf.getwallcrossings(1), conf.getwallcrossings(2), ps.timestep, ps.urange, ps.ustrength, rodDist, ps.particlesize, ps.runs,
    //        sizeOfArray(ps.timestep), sizeOfArray(ps.urange), sizeOfArray(ps.ustrength), sizeOfArray(rodDist), sizeOfArray(ps.particlesize), potentialMod, ps.includeSteric);


    cout << "Simulation Finished" << endl;

    if (ps.setPBC)  fclose(snapFile);
    trajectoryfile.close();
    rodPosFileVMD.close();
    rodsfile.close();
    // TODO distancefile
    distancesfile.close();

    return 0;
}
