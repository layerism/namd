
#ifdef NAMD_TCL

int ScriptTcl::Tcl_loadAdaptCgParam(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[]) {
	ScriptTcl *script = (ScriptTcl *)clientData;
	script->initcheck();
	if (argc != 2) {
		Tcl_AppendResult(interp, "usage: loadAdaptCgParam <filename>", NULL);
		return TCL_ERROR;
	}

	Node::Object()->loadAdaptCgParam(argv[1]);

	script->runController(SCRIPT_LOADCGPARAM);

	return TCL_OK;

}

#endif

