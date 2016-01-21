/**
 ***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
 ***  The Board of Trustees of the University of Illinois.
 ***  All rights reserved.
 **/

#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "converse.h"
#include "Node.h"
#include "common.h"
#include <fstream>

class ControllerBroadcasts;
class NamdState;
class SimParameters;
class RequireReduction;

#ifdef MEM_OPT_VERSION
class CollectionMasterHandler;
#else
class CollectionMaster;
#endif

class Random;
class PressureProfileReduction;

class Controller {
public:
	Controller(NamdState *s);
	virtual ~Controller(void);
	void run(void);             // spawn thread, etc.
	void awaken(void) {
		CthAwaken(thread);
	}
	;
	void resumeAfterTraceBarrier(int);
#ifdef MEASURE_NAMD_WITH_PAPI
	void resumeAfterPapiMeasureBarrier(int step);
#endif

protected:
	friend class ScriptTcl;
	virtual void algorithm(void);	// subclasses redefine this method

	void integrate(); // Verlet integrator
	void minimize(); // CG minimizer
	RequireReduction *min_reduction;

	void receivePressure(int step, int minimize = 0);
	Tensor pressure_normal;
	Tensor pressure_nbond;
	Tensor pressure_slow;
	Tensor virial_amd;
	//lyk
	Tensor virial_acg;
	Tensor virial_its;
	//lyk
	Tensor groupPressure_normal;
	Tensor groupPressure_nbond;
	Tensor groupPressure_slow;
	Tensor controlPressure_normal;
	Tensor controlPressure_nbond;
	Tensor controlPressure_slow;
	int nbondFreq;
	int slowFreq;
	BigReal temp_avg;
	BigReal pressure_avg;
	BigReal groupPressure_avg;
	int avg_count;
	Tensor pressure_tavg;
	Tensor groupPressure_tavg;
	int tavg_count;
	void compareChecksums(int, int = 0);
	int computeChecksum;
	int marginViolations;
	int pairlistWarnings;
	void printTiming(int);
	void printMinimizeEnergies(int);
	BigReal min_energy;
	BigReal min_f_dot_f;
	BigReal min_f_dot_v;
	BigReal min_v_dot_v;
	int min_huge_count;
	void printDynamicsEnergies(int);
	void printEnergies(int step, int minimize);
	int numDegFreedom;
	BigReal totalEnergy;
	BigReal electEnergy;
	BigReal electEnergySlow;
	BigReal ljEnergy;
	BigReal groLJEnergy;
	BigReal groGaussEnergy;
	BigReal goNativeEnergy;
	BigReal goNonnativeEnergy;
	BigReal goTotalEnergy;
//fepb
	BigReal electEnergy_f;
	BigReal electEnergySlow_f;
	BigReal ljEnergy_f;
	BigReal ljEnergy_f_left;	// used by WCA repulsive, [s1,s2]
	BigReal exp_dE_ByRT;
	BigReal net_dE;
	BigReal dG;
	int FepNo;
	void printFepMessage(int);
	BigReal fepSum;
//fepe

	BigReal electEnergy_ti_1;
	BigReal electEnergySlow_ti_1;
	BigReal ljEnergy_ti_1;
	BigReal electEnergy_ti_2;
	BigReal electEnergySlow_ti_2;
	BigReal ljEnergy_ti_2;
	BigReal net_dEdl_elec_1;
	BigReal net_dEdl_elec_2;
	BigReal net_dEdl_lj_1;
	BigReal net_dEdl_lj_2;
	BigReal electEnergyPME_ti_1;
	BigReal electEnergyPME_ti_2;
	int TiNo;
	BigReal recent_dEdl_elec_1;
	BigReal recent_dEdl_elec_2;
	BigReal recent_dEdl_lj_1;
	BigReal recent_dEdl_lj_2;
	int recent_TiNo;
	void printTiMessage(int);

	BigReal drudeComTemp;  // temperature of Drude centers of mass
	BigReal drudeBondTemp; // temperature of Drude bonds
	BigReal drudeComTempAvg;
	BigReal drudeBondTempAvg;

	BigReal kineticEnergy;
	BigReal kineticEnergyHalfstep;
	BigReal kineticEnergyCentered;
	BigReal temperature;
	BigReal smooth2_avg;
	BigReal smooth2_avg2;  // avoid internal compiler error
	Tensor pressure;
	Tensor groupPressure;
	int controlNumDegFreedom;
	Tensor controlPressure;
	void enqueueCollections(int);
	void correctMomentum(int step);
	void rescaleVelocities(int);
	BigReal rescaleVelocities_sumTemps;
	int rescaleVelocities_numTemps;
	void reassignVelocities(int);
	void tcoupleVelocities(int);
	void berendsenPressure(int);
	Tensor berendsenPressure_avg;
	int berendsenPressure_count;
	void langevinPiston1(int);
	void langevinPiston2(int);
	Tensor langevinPiston_strainRate;
	Tensor strainRate_old;  // for langevinPistonBarrier no
	Tensor positionRescaleFactor;  // for langevinPistonBarrier no

	int ldbSteps;
	void rebalanceLoad(int);
	int fflush_count;
	void cycleBarrier(int, int);

	void traceBarrier(int, int);

#ifdef MEASURE_NAMD_WITH_PAPI
	void papiMeasureBarrier(int, int);
#endif

	// void suspend(void) { CthSuspend(); };
	void terminate(void);

	Random *random;
	SimParameters * const simParams;	// for convenience

	NamdState * const state;		// access data in state
	RequireReduction *reduction;
	RequireReduction *amd_reduction;

	// data for pressure profile reductions and output
	PressureProfileReduction *ppbonded;
	PressureProfileReduction *ppnonbonded;
	PressureProfileReduction *ppint;
	int pressureProfileSlabs;
	int pressureProfileCount;
	BigReal *pressureProfileAverage;

	CollectionMaster * const collection;

	ControllerBroadcasts * broadcast;
	std::ofstream xstFile;
	void outputExtendedSystem(int step);
	void writeExtendedSystemLabels(std::ofstream &file);
	void writeExtendedSystemData(int step, std::ofstream &file);

//fepb
	std::ofstream fepFile;
	void outputFepEnergy(int step);
	void writeFepEnergyData(int step, std::ofstream &file);
//fepe
	std::ofstream tiFile;
	void outputTiEnergy(int step);
	void writeTiEnergyData(int step, std::ofstream &file);

	// for checkpoint/revert
	int checkpoint_stored;
	Lattice checkpoint_lattice;
	Tensor checkpoint_langevinPiston_strainRate;
	Tensor checkpoint_berendsenPressure_avg;
	int checkpoint_berendsenPressure_count;
	BigReal checkpoint_smooth2_avg;

//for accelMD
	void rescaleaccelMD(int step, int minimize = 0);
	BigReal accelMDdVAverage;

	//lyk

#define MAXLOOP 1000
#define WHAMFLAG
	BigReal DENSITYALAPHA;
	int SKIP;
	BigReal LOGZERROR;

	//for ITS
	BigReal its_lj;
	BigReal its_replj;
	BigReal its_elect;

	BigReal *its_replicaWeight;
	int its_replicaNum;
	int its_samplesNum;
	BigReal *its_logZ;
	BigReal *its_meanCV;
	BigReal *its_sgmCV;
	BigReal *its_energySamples;
	BigReal *its_effCvSamples;
	BigReal its_lbdMin;
	BigReal its_lbdMax;
	BigReal *its_lbdList;
	Vector its_rescalingFactor;
#ifdef WHAMFLAG
	//wham
	std::ofstream its_whamLog;
	BigReal **its_cvEff;
	BigReal **its_logZHis;
	BigReal *its_logZEffHis;
	BigReal *its_whamLogZ;
	BigReal *its_energySamples_wham;
	int its_samplesNum_wham;
	void its_whamEstimate(int step, BigReal **U_eff, BigReal *est_logZ, BigReal *mean, BigReal *sgm);
	BigReal its_whamEstimateLogZEff(int step, BigReal **U_eff, BigReal *logZEffHis);
	void its_cvEffHis(int step, BigReal **U_eff);
#endif
	BigReal *its_lbdHistogram;
	BigReal *its_dynamicWeight;
	BigReal *its_dynamicWeightSum;
	int its_intType;
	std::ofstream its_updateLog;
	std::ifstream its_inFile;
	Bool its_inFileOpen;
	void itsMD(int step, int minimize = 0);
	void its_updateVirial(Vector scale);
	void its_selectEnergy(BigReal &energy);
	void its_init(void);
	BigReal its_lbd(int i);
	BigReal its_cv(BigReal energy, BigReal lbd);
	BigReal its_eff_cv(BigReal energy, BigReal *logZ);
	Vector its_eff_lbd(BigReal energy, BigReal *logZ);
	void its_estimate(BigReal *mean, BigReal *sgm, BigReal *logz);
	void its_eflatUpdateLogZ(int step, BigReal *logz, BigReal *weight);

	// for adaptcg
	BigReal aa_lj;
	BigReal aa_replj;
	BigReal aa_elect;
	BigReal cg_lj;
	BigReal cg_replj;
	BigReal cg_elect;
	BigReal *acg_replicaWeight;
	BigReal *acg_logZ;
	BigReal *acg_lbdHistogram;
	int acg_replicaNum;
	int acg_samplesNum;
	int acg_intType;
	BigReal *acg_aaEnergySamples;
	BigReal *acg_cgEnergySamples;
	BigReal acg_lbdMin;
	BigReal acg_lbdMax;
	BigReal *acg_lbdList;
	Vector acg_rescalingFactor;
	std::ofstream acg_updateLog;
	std::ifstream acg_inFile;
	Bool acg_inFileOpen;

#ifdef WHAMFLAG
	//wham
	std::ofstream acg_whamLog;
	BigReal **acg_cvEff;
	BigReal **acg_logZHis;
	BigReal *acg_logZEffHis;
	BigReal *acg_whamLogZ;
	BigReal *acg_aaEnergySamples_wham;
	BigReal *acg_cgEnergySamples_wham;
	int acg_samplesNum_wham;
	void acg_whamEstimate(int step, BigReal **U_eff, BigReal *est_logZ, BigReal *mean, BigReal *sgm);
	BigReal acg_whamEstimateLogZEff(int step, BigReal **U_eff, BigReal *logZEffHis);
	void acg_cvEffHis(int step, BigReal **U_eff);
#endif
	//non wham
	BigReal *acg_effCvSamples;
	BigReal *acg_meanCV;
	BigReal *acg_sgmCV;
	BigReal *acg_dynamicWeight;
	BigReal *acg_dynamicWeightSum;
	void adaptCgMD(int step, int minimize = 0);
	void acg_init(void);
	void acg_updateVirial(Vector scale);
	void acg_selectEnergy(BigReal &aaEnergy, BigReal &cgEnergy);
	BigReal acg_lbd(int i);
	BigReal acg_lbdaa(BigReal lbd);
	BigReal acg_lbdcg(BigReal lbd);
	BigReal acg_cv(BigReal aa_energy, BigReal cg_energy, BigReal lbd);
	BigReal acg_eff_cv(BigReal aa_energy, BigReal cg_energy, BigReal *logZ);
	Vector acg_eff_lbd(BigReal aa_energy, BigReal cg_energy, BigReal *logZ);
	void acg_estimate(BigReal *mean, BigReal *sgm, BigReal *logz);
	void acg_eflatUpdateLogZ(int step, BigReal *logz, BigReal *weight);
	void acg_purifyLogZ(int step, BigReal *logZ);
	//lyk

//JS for adaptive temperature sampling
	void adaptTempInit(int step);
	void adaptTempUpdate(int step, int minimize = 0);
	void adaptTempWriteRestart(int step);

	BigReal *adaptTempPotEnergyAveNum;
	BigReal *adaptTempPotEnergyAveDen;
	BigReal *adaptTempPotEnergyVarNum;
	BigReal *adaptTempPotEnergyAve;
	BigReal *adaptTempPotEnergyVar;
	int *adaptTempPotEnergySamples;
	BigReal *adaptTempBetaN;
	BigReal adaptTempT;
	BigReal adaptTempDTave;
	BigReal adaptTempDTavenum;
	BigReal adaptTempBetaMin;
	BigReal adaptTempBetaMax;
	int adaptTempBin;
	int adaptTempBins;
	BigReal adaptTempDBeta;
	BigReal adaptTempCg;
	BigReal adaptTempDt;
	Bool adaptTempAutoDt;
	BigReal adaptTempDtMin;
	BigReal adaptTempDtMax;
	std::ofstream adaptTempRestartFile;

private:
	CthThread thread;
	static void threadRun(Controller*);

	double startCTime;
	double startWTime;
	double firstCTime;
	double firstWTime;
	double startBenchTime;

	int computesPartitioned;
};

#endif // CONTROLLER_H

