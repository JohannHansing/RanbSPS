#include "headers/SingleParticleSimulationRanb.h"

using namespace std;



int main(int argc, const char* argv[]){
    //Main includes the iteration loops for the simulation

    //NOTE: so far wallcrossings is added for all runs!!! makes it kind of incorrect, since after each run, ppos is reset.
    //NOTE: so far saving Instant Values for each tenth step!

    //TRIGGERS:
    string distribution = argv[1];    // TODO
    bool ranRod = (strcmp(argv[2] , "true") == 0 ) ;
    bool potentialMod = (strcmp(argv[3] , "true") == 0 ) ;       //BESSEL TODO
    bool recordMFP = (strcmp(argv[4] , "true") == 0 ) ;
    bool recordPosHisto = (strcmp(argv[5] , "true") == 0 ) ;
    bool includeSteric = (strcmp(argv[6] , "true") == 0 ) ;  // steric 2
    bool ranPot = (strcmp(argv[7] , "true") == 0 ) ;
    bool hpi = (strcmp(argv[8] , "true") == 0 ) ;          // hpi exp
    int boolpar = 8;
    ifdebug(cout << "copied bools. ";)
    
    // Checking for correct structure of input arguments
    for (int k= 1; k < argc; k++ ) cout << "parameter " << k << " " << argv[k] << endl;
    for (int b_i=2; b_i<=boolpar; b_i++){
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
    unsigned int saveInt;
    int instValIndex;                             //Counter for addInstantValue
    
    ifdebug(cout << "copied  params. ";)
    
    cout << "distribution " << distribution << endl;


    
    steps = simtime/timestep;
    saveInt = steps/instantvalues;
    const int trajout = (int)(10/timestep);
        
    //Create data folders and print location as string to string "folder"
    string folder = createDataFolder(distribution, timestep, simtime, urange, ustrength, particlesize, includeSteric, ranPot, ranRod);
    ifdebug(cout << "created folder. ";)


    //initialize averages
    CAverage energyU = CAverage("Upot", folder, instantvalues, runs);
    CAverage squareDisp = CAverage("squaredisp", folder, instantvalues, runs);
    ifdebug(cout << "created CAverage files. ";)

    //initialize instance of configuration
    CConfiguration conf = CConfiguration(distribution,timestep, urange, ustrength, potentialMod, particlesize, recordPosHisto, includeSteric, ranPot, hpi, ranRod);
    ifdebug(cout << "created CConf conf. ";)

    
    unsigned int stepcount = 0;
    ofstream trajectoryfile;
    trajectoryfile.open((folder + "/Coordinates/trajectory.txt").c_str());


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
            conf.checkBoxCrossing(); //check if particle has crossed the confinement of the box
            
            stepcount++;
            if (stepcount%trajout == 0) {
                std::vector<double> ppos = conf.getppos();
                trajectoryfile << fixed << stepcount * timestep << "\t" << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;
                ifdebug(cout << stepcount * timestep << "\t" << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;)
            }
        }
        
        
    }//----------END OF RUNS-LOOP
    


    //watch out: Save average instant values at timeInterval: timestep * saveinterval saveInt!!!
    energyU.saveAverageInstantValues(saveInt*timestep);
    squareDisp.saveAverageInstantValues(saveInt*timestep);

    


    //printReport(ranRod, conf.getwallcrossings(0), conf.getwallcrossings(1), conf.getwallcrossings(2), timestep, urange, ustrength, rodDist, particlesize, runs,
    //        sizeOfArray(timestep), sizeOfArray(urange), sizeOfArray(ustrength), sizeOfArray(rodDist), sizeOfArray(particlesize), potentialMod, includeSteric);

    
    cout << "Simulation Finished" << endl;
    
    //If settingsFile is saved, then the simulation was successfull
    settingsFile(folder, ranRod, particlesize, timestep, runs, steps, ustrength, urange, potentialMod, recordMFP, includeSteric, ranPot, hpi, distribution);
	
    trajectoryfile.close();
	
    return 0;
}
