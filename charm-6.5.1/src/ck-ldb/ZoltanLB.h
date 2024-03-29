/** \file ZoltanLB.h
 *
 * Load balancer using Zoltan hypergraph partitioner. This is a multicast aware
 * load balancer
 * Harshitha, 2012/02/21
 */

/**
 * \addtogroup CkLdb
*/
/*@{*/

#ifndef _ZOLTANLB_H_
#define _ZOLTANLB_H_

#include "CentralLB.h"
#include "ZoltanLB.decl.h"

void CreateZoltanLB();
BaseLB * AllocateZoltanLB();

class ZoltanLB : public CentralLB {
public:
  ZoltanLB(const CkLBOptions &);
  ZoltanLB(CkMigrateMessage *m):CentralLB(m) { lbname = "ZoltanLB"; }
private:
  CmiBool QueryBalanceNow(int step) { return CmiTrue; }
  void work(LDStats* stats);
};

#endif /* _ZOLTANLB_H_ */

/*@}*/
