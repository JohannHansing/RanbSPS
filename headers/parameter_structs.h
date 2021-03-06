//
//  paramater_structs.h
//
//
//
//

#ifndef ____paramater_structs__ //checks if already defined by other header file
#define ____paramater_structs__


#include <string.h>

struct paramstruct {
    string distribution;
    bool ranRod;
    bool rand;
    bool setPBC;
    bool Tmp4;
    bool includeSteric;
    bool ranU;
    bool Pointq;
    bool mixU;

    int runs;
    double timestep;
    int simtime;
    int instantvalues;
    unsigned int steps;

    double particlesize;
    double urange;
    double ustrength;
    double dvar;
    double polydiam;
    double drqop;
    double Cratio; // = Catt/Crep
    double uratio; //  = Uatt/Urep
    int N_patches; // Number of patches along a "full" fiber (i.e. from end to end of simulation box)
    int bol; // the boxsize over the patchlength b/l
};



#endif /* defined(____paramater_structs__) */ //end check
