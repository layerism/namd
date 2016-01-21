#include "SimParameters.h"
#include "Parameters.h"
#include "Molecule.h"
#include "BackEnd.h"
/*
class Molecule;
class Parameters;
class SimParameters;

class HREMD {
public:
	HREMD(void) {
		SimParameters *simParam = Node::Object()->simParameters;
		Molecule *mol = Node::Object()->molecule;
		Parameters *param = Node::Object()->parameters;
	}
	~HREMD(void) {

	}
private:
	SimParameters *simParam;
	Molecule *mol;
	Parameters *param;	
};*/

#ifdef NAMD_TCL

static int Tcl_hremdinit(ClientData, Tcl_Interp *interp, int argc, char *argv[]) {
	SimParameters *simParam = Node::Object()->simParameters;
	Molecule *mol = Node::Object()->molecule;
	Parameters *param = Node::Object()->parameters;
	//simParam->dt=0.0;
	return TCL_OK;
}

static int Tcl_gettimestep(ClientData, Tcl_Interp *interp, int argc, char *argv[]) {
	SimParameters *simParam = Node::Object()->simParameters;
	char s[1024];
	simParam->dt=0.0;
	sprintf(s,"%g", simParam->dt);
	Tcl_SetResult(interp,s,TCL_VOLATILE);
	return TCL_OK;
}

#endif

