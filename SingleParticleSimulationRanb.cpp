#include "headers/SingleParticleSimulationRanb.h"


using namespace std;


int main(int argc, const char* argv[]){
    //Main includes the iteration loops for the simulation

    //NOTE: so far wallcrossings is added for all runs!!! makes it kind of incorrect, since after each run, ppos is reset.
    //NOTE: so far saving Instant Values for each tenth step!

    //TRIGGERS:
    bool writeTrajectory = (strcmp(argv[1] , "true") == 0 ) ;    // relative position TODO
    bool ranRod = (strcmp(argv[2] , "true") == 0 ) ;
    bool potentialMod = (strcmp(argv[3] , "true") == 0 ) ;       //BESSEL TODO
    bool recordMFP = (strcmp(argv[4] , "true") == 0 ) ;
    bool recordPosHisto = (strcmp(argv[5] , "true") == 0 ) ;
    bool includeSteric = (strcmp(argv[6] , "true") == 0 ) ;  // steric 2
    bool ranPot = (strcmp(argv[7] , "true") == 0 ) ;
    bool hpi = (strcmp(argv[8] , "true") == 0 ) ;          // hpi exp
    int boolpar = 8;
	
	// Checking for correct structure of input arguments
	for (int k= 0; k < argc; k++ ) cout << "parameter " << k << " " << argv[k] << endl;
	for (int b_i=1; b_i<=boolpar; b_i++){
		if ((strcmp(argv[b_i] , "true") == 1 )  && (strcmp(argv[b_i] , "false") == 1 )){
			cerr << "Error; Bool parameter " << b_i << " is not either 'true' or 'false'!" << endl;
			exit(1);
		}
	}


    int runs = atoi( argv[boolpar+1] );                       // Number of Simulation runs to get mean values from
    double timestep = atof( argv[boolpar+2] );
    int simtime = atoi( argv[boolpar+3] );                   // simulation time
    int instantvalues = 200;
    unsigned int steps;
    
    double rodDist = atof( argv[boolpar+4] );                //Distance of the rods depending on boxsize. for zero do rodDist[]={0.0}
    double boxsize = atof( argv[boolpar+5] );
    double particlesize = atof( argv[boolpar+6] );
    double urange = atof( argv[boolpar+7] );
    double ustrength = atof( argv[boolpar+8] );
    double hpi_u = atof( argv[boolpar+9] );
    double hpi_k = atof( argv[boolpar+10] );
    double n_rods = atof( argv[boolpar+11] );
    unsigned int saveInt;
    int instValIndex;                             //Counter for addInstantValue
    
    if (ranRod && (ranPot || includeSteric)){
        cout << "ERROR !!!!!!!!!!!!!!!!!!! This is not yet supported!" << endl;
        return 3;
    }

	
	
    //MFP
    double fpInt = boxsize/10;

    //printReport(ranRod, conf.getwallcrossings(0), conf.getwallcrossings(1), conf.getwallcrossings(2), timestep, urange, ustrength, rodDist, particlesize, runs,
    //            sizeOfArray(timestep), sizeOfArray(urange), sizeOfArray(ustrength), sizeOfArray(rodDist), sizeOfArray(particlesize), potentialMod, includeSteric);


    
    steps = simtime/timestep;
    saveInt = steps/instantvalues;
    const int trajout = (int)(10/timestep);
        
    //Create data folders and print location as string to string "folder"
    string folder = createDataFolder(ranRod, timestep, simtime, urange, ustrength, boxsize, particlesize, rodDist, potentialMod, includeSteric, ranPot, hpi, hpi_u, hpi_k, n_rods);


    //initialize averages
    CAverage energyU = CAverage("Upot", folder, instantvalues, runs);
    CAverage squareDisp = CAverage("squaredisp", folder, instantvalues, runs);
    //CAverage displacem = CAverage("displacement", folder, instantvalues, runs);
    //CAverage squareDisp_x = CAverage("squaredisp_x", folder, instantvalues, runs);
    //CAverage squareDisp_y = CAverage("squaredisp_y", folder, instantvalues, runs);
    //CAverage squareDisp_z = CAverage("squaredisp_z", folder, instantvalues, runs);
    //CAverage mfp_x = CAverage("mfp_x", folder, 1, 1);
    CAverage mfp_xyz;
    if ( recordMFP ) mfp_xyz = CAverage("mfp_xyz", folder, 1, 1);


    //initialize instance of configuration
    CConfiguration conf = CConfiguration(timestep, urange, ustrength, boxsize, rodDist, potentialMod, particlesize, recordPosHisto, includeSteric,
    ranPot, hpi , hpi_u, hpi_k, ranRod, n_rods);


    //create file to save the trajectory
    string traj_file = folder + "/Coordinates/single_traj.xyz";
    if (writeTrajectory) conf.saveXYZTraj(traj_file,0,"w");
    
    unsigned int stepcount = 0;
    ofstream trajectoryfile;
    trajectoryfile.open((folder + "/Coordinates/trajectory.txt").c_str());


    //cout << "Starting Run Number: " << simcounter << " out of " << totalsims << endl;
    cout << "Starting Simulation!" << endl;


// **************START OF RUNS-LOOP
    for (int l = 0; l<runs; l++){

        conf.updateStartpos();

        instValIndex = 0;
        int fpCounter[3] = {0};                  //counter for first passage times (next point to pass first: fpCounter*fpInt


        for (int i = 0; i < steps; i++){  //calculate stochastic force first, then mobility force!!


            conf.calcStochasticForces();

            conf.calcMobilityForces();


            if (((i+1)%100 == 0) && (l == 0) && writeTrajectory){       //Save the first trajectory to file
                conf.saveXYZTraj(traj_file, i, "a");                    // TODO change back ((i+1)%XXX == 0) to 100
            }




            if (((i+1)%saveInt) == 0){       //saving Instant Values for each saveInt'th step!
                energyU.addInstantValue(conf.getUpot(), instValIndex);
                squareDisp.addInstantValue(conf.getPosVariance(), instValIndex);
            //    displacem.addInstantValue(conf.getDisplacement(), instValIndex);
            //    squareDisp_x.addInstantValue(conf.get1DPosVariance(0), instValIndex);
            //    squareDisp_y.addInstantValue(conf.get1DPosVariance(1), instValIndex);
            //    squareDisp_z.addInstantValue(conf.get1DPosVariance(2), instValIndex);

                instValIndex += 1;
            }

            /*if ((conf.checkFirstPassage(fpInt*(fpCounter+1))) && recordMFP) {
                fpCounter+=1;
                mfp_x.addFPValue(i*timestep, fpCounter);
            }
            */
            if (recordMFP){
                for (int a=0; a < 3; a++){
                    if (conf.checkFirstPassage(fpInt*(fpCounter[a]+1), a)) {
                        fpCounter[a]+=1;
                        mfp_xyz.addFPValue(i*timestep, fpCounter[a]);
                    }
                }
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
            conf.checkBoxCrossing(); //check if particle has crossed the confinement of the box
            
            stepcount++;
            if (stepcount%trajout == 0) {
                std::vector<double> ppos = conf.getppos();
                trajectoryfile << fixed << stepcount * timestep << "\t" << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;
            }



            if (((i % 5) == 0) && recordPosHisto) conf.addHistoValue();

        }
        if ( recordPosHisto ) conf.printHistoMatrix(folder);
        
        
    }//----------END OF RUNS-LOOP
    
    conf.printAvRods();




    //watch out: Save average instant values at timeInterval: timestep * saveinterval saveInt!!!
    energyU.saveAverageInstantValues(saveInt*timestep);
    squareDisp.saveAverageInstantValues(saveInt*timestep);
    //displacem.saveAverageInstantValues(saveInt*timestep);
    //squareDisp_x.saveAverageInstantValues(saveInt*timestep);
    //squareDisp_y.saveAverageInstantValues(saveInt*timestep);
    //squareDisp_z.saveAverageInstantValues(saveInt*timestep);
    //if (recordMFP) mfp_x.saveAverageFPValue(fpInt);
    if (recordMFP) mfp_xyz.saveAverageFPValue(fpInt);

    


    //printReport(ranRod, conf.getwallcrossings(0), conf.getwallcrossings(1), conf.getwallcrossings(2), timestep, urange, ustrength, rodDist, particlesize, runs,
    //        sizeOfArray(timestep), sizeOfArray(urange), sizeOfArray(ustrength), sizeOfArray(rodDist), sizeOfArray(particlesize), potentialMod, includeSteric);

    
    cout << "Simulation Finished" << endl;
    
    //If settingsFile is saved, then the simulation was successfull
    settingsFile(folder, ranRod, particlesize, boxsize, timestep, runs, steps, ustrength, urange, rodDist, potentialMod, recordMFP, includeSteric, ranPot ,hpi, hpi_u, hpi_k, n_rods);
	
    trajectoryfile.close();
	
    return 0;
}
