inline void Sequencer::init_lyk_HREMD(void) {
	reduction_rev = ReductionMgr::Object()->willRequire(REDUCTIONS_BASIC);
	bond_array_origin = new BondValue[param->NumBondParams];
	angle_array_origin = new AngleValue[param->NumAngleParams];

	memcpy(bond_array_origin,param->bond_array,sizeof(BondValue)*(param->NumBondParams));
	memcpy(angle_array_origin,param->angle_array,sizeof(AngleValue)*(param->NumAngleParams));

}
inline void Sequencer::end_lyk_HREMD(void) {
	delete []bond_array_origin;
	delete []angle_array_origin;
	delete reduction_rev;
}
void Sequencer::getRescaleBondIndex() {

}

