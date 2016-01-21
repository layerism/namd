/**
 ***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
 ***  The Board of Trustees of the University of Illinois.
 ***  All rights reserved.
 **/

#ifndef PATCHTYPES_H
#define PATCHTYPES_H

#include "NamdTypes.h"
#include "Lattice.h"

class Flags {
public:
	int step;			// timestep number reported to user
	// Same number may appear multiple times!
	int sequence;			// sequence number of compute call
	// changes by 1 every time!
	int doEnergy;
	int doNonbonded;
	int doFullElectrostatics;
	int doMolly;
	// BEGIN LA
	int doLoweAndersen;
	// END LA
	int doGBIS;				// gbis
	int doLCPO;				// LCPO
	int submitLoadStats;
	int maxForceUsed;		// may ignore slower force classes
	int maxForceMerged;		// add this and faster to normal

	int usePairlists;
	int savePairlists;
	BigReal pairlistTolerance;
	BigReal maxAtomMovement;
	BigReal maxGroupRadius;

	Lattice lattice;		// rather than shipping around separately

	inline Flags& operator=(const Flags &flags) {
		CmiMemcpy(this, &flags, sizeof(Flags));
		return *this;
	}

};

class Results {
public:

	enum {
		normal = 0,
		nbond = 1,
		slow = 2,
		amdf = 3,
		vdw11 = 4,
		vdw12 = 5,
		vdw22 = 6,
		vdw11_rep = 7,
		vdw12_rep = 8,
		vdw22_rep = 9,
		elect11 = 10,
		elect12 = 11,
		elect22 = 12,
		pme = 13,
		maxNumForces = 14 //lyk
	};
	//enum { normal=0, nbond=1, slow=2, amdf=3,maxNumForces=4 };
	Force *f[maxNumForces];
};

#endif

