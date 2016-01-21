/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEMGR_H
#define COMPUTEMGR_H

#include "charm++.h"
#include "main.h"
#include <new>

#include "NamdTypes.h"
#include "BOCgroup.h"

#include "ResizeArray.h"

#include "GlobalMaster.h"
#include "GlobalMasterServer.h"
#include "ComputeMgr.decl.h"

class Compute;
class ComputeMap;
class CkQdMsg;

class ComputeGlobal;
class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;

class ComputeDPME;
class ComputeDPMEDataMsg;
class ComputeDPMEResultsMsg;
class ComputeConsForceMsg;

class ComputeEwald;
class ComputeEwaldMsg;

class ComputeNonbondedCUDA;
class NonbondedCUDASlaveMsg;

class ComputeNonbondedMIC;
class NonbondedMICSlaveMsg;

class ComputeNonbondedWorkArrays;

class ComputeMgr : public CBase_ComputeMgr
{
public:

  ComputeMgr();
  ~ComputeMgr();
  void createComputes(ComputeMap *map);
  void updateComputes(int,CkGroupID);
  void updateComputes2(CkQdMsg *);
  void updateComputes3();
  void splitComputes();
  void splitComputes2(CkQdMsg *);
  void updateLocalComputes();
  void updateLocalComputes2(CkQdMsg *);
  void updateLocalComputes3();
  void updateLocalComputes4(CkQdMsg *);
  void updateLocalComputes5();
  void doneUpdateLocalComputes();

  void sendComputeGlobalConfig(ComputeGlobalConfigMsg *);
  void recvComputeGlobalConfig(ComputeGlobalConfigMsg *);
  void sendComputeGlobalData(ComputeGlobalDataMsg *);
  void recvComputeGlobalData(ComputeGlobalDataMsg *);
  void sendComputeGlobalResults(ComputeGlobalResultsMsg *);
  void recvComputeGlobalResults(ComputeGlobalResultsMsg *);
  void enableComputeGlobalResults();

  void sendComputeDPMEData(ComputeDPMEDataMsg *);
  void recvComputeDPMEData(ComputeDPMEDataMsg *);
  void sendComputeDPMEResults(ComputeDPMEResultsMsg *, int);
  void recvComputeDPMEResults(ComputeDPMEResultsMsg *);

  void sendComputeEwaldData(ComputeEwaldMsg *);
  void recvComputeEwaldData(ComputeEwaldMsg *);
  void sendComputeEwaldResults(ComputeEwaldMsg *);
  void recvComputeEwaldResults(ComputeEwaldMsg *);

  void recvComputeConsForceMsg(ComputeConsForceMsg *);

  // Made public in order to access the ComputeGlobal on the node
  ComputeGlobal *computeGlobalObject; /* node part of global computes */
  ResizeArray<ComputeGlobalResultsMsg*> computeGlobalResultsMsgs;
  int computeGlobalResultsMsgSeq;
  int computeGlobalResultsMsgMasterSeq;

  void sendYieldDevice(int pe);
  void recvYieldDevice(int pe);
  void sendBuildCudaForceTable();
  void recvBuildCudaForceTable();

  // DMK
  void sendBuildMICForceTable();
  void recvBuildMICForceTable();

  void sendCreateNonbondedCUDASlave(int,int);
  void recvCreateNonbondedCUDASlave(NonbondedCUDASlaveMsg *);
  void sendNonbondedCUDASlaveReady(int,int,int,int);
  void recvNonbondedCUDASlaveReady(int,int,int);
  void sendNonbondedCUDASlaveEnqueue(ComputeNonbondedCUDA *c, int,int,int,int);
  void sendCreateNonbondedMICSlave(int,int);
  void recvCreateNonbondedMICSlave(NonbondedMICSlaveMsg *);
  void sendNonbondedMICSlaveReady(int,int,int,int);
  void recvNonbondedMICSlaveReady(int,int,int);
  void sendNonbondedMICSlaveEnqueue(ComputeNonbondedMIC *c, int,int,int,int);
  void sendMICPEData(int,int);
  void recvMICPEData(int,int);
  int isMICProcessor(int);
  
private:
  void createCompute(ComputeID, ComputeMap *);

  GlobalMasterServer *masterServerObject; /* master part of global computes */
  ComputeDPME *computeDPMEObject;

  ComputeEwald *computeEwaldObject;

  ComputeNonbondedCUDA *computeNonbondedCUDAObject;

  ComputeNonbondedMIC *computeNonbondedMICObject;

  ComputeNonbondedWorkArrays *computeNonbondedWorkArrays;

  int skipSplitting;
  int updateComputesCount;
  int updateComputesReturnEP;
  CkGroupID updateComputesReturnChareID;

  ResizeArray<int> computeFlag;

  int* micPEData;  // DMK : An array (1 bit per PE) which will hold a flag indicating if a given PE is driving a MIC card or not
};


#endif /* COMPUTEMGR_H */
