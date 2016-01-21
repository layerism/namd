inline BigReal exp_cut(BigReal x) {

	if (x >= -300 && x <= 300) {
		return (exp(x));
	} else if (x >= 300) {
		return (exp(300));
	} else {
		return (0.0);
	}
}

inline BigReal get_max(BigReal *x, int size) {
	BigReal max = -9999999999999999.0;
	for (int i = 0; i < size; ++i) {
		if (max < x[i])
			max = x[i];
	}
	return (max);
}

#define WHAMFLAG

//dym build A[N][M]
#define REALARRAY2D(A,N,M) \
	A = new BigReal *[N]; \
	for (int i=0; i<N; ++i) { \
		A[i] = new BigReal [M]; \
	} \
	for (int i=0; i<N; ++i) { \
		for (int j=0; j<M; ++j) { \
			A[i][j] = 0.0; \
		} \
	}
#define DELETE2D(A,N,M) \
	for (int i=0; i<N; ++i) { \
		delete[] A[i]; \
	} \
	delete[] A

#define REALARRAY1D(A,N) \
	A = new BigReal [N]; \
	for (int i=0; i<N; ++i) { \
		A[i] = 0.0; \
	}

//==============================================================
//
//			MODULE FOR ACG
//
//==============================================================
inline void Controller::acg_updateVirial(Vector scale) {
	Tensor virial_aa;
	Tensor virial_cg;

	if (simParams->aaGroup == -1) {
		//Tensor virial_dihe;
		Tensor virial_nbond;
		Tensor virial_slow;
		//GET_TENSOR(virial_dihe, amd_reduction, REDUCTION_VIRIAL_AMD_DIHE);
		GET_TENSOR(virial_nbond, amd_reduction, REDUCTION_VIRIAL_NBOND);
		GET_TENSOR(virial_slow, amd_reduction, REDUCTION_VIRIAL_SLOW);

		//virial_aa = (virial_dihe + virial_nbond + virial_slow);
		virial_aa = (virial_nbond + virial_slow);
	} else {
		Tensor virial_11_vdw;
		Tensor virial_11_repvdw;
		Tensor virial_11_elect;
		GET_TENSOR(virial_11_vdw, amd_reduction, REDUCTION_VIRIAL_11_VDW);
		GET_TENSOR(virial_11_repvdw, amd_reduction, REDUCTION_VIRIAL_11_REPVDW);
		GET_TENSOR(virial_11_elect, amd_reduction, REDUCTION_VIRIAL_11_ELECT);
		if (acg_intType == 0) {
			virial_aa = virial_11_vdw;
		} else if (acg_intType == 1) {
			virial_aa = virial_11_elect;
		} else if (acg_intType == 2) {
			virial_aa = (virial_11_vdw + virial_11_elect);
		} else if (acg_intType == 3) {
			Tensor virial_11_attvdw = virial_11_vdw - virial_11_repvdw;
			virial_aa = virial_11_attvdw;
		} else if (acg_intType == 4) {
			Tensor virial_11_attvdw = virial_11_vdw - virial_11_repvdw;
			virial_aa = (virial_11_attvdw + virial_11_elect);
		}
	}

	Tensor virial_22_vdw;
	Tensor virial_22_repvdw;
	Tensor virial_22_elect;
	GET_TENSOR(virial_22_vdw, amd_reduction, REDUCTION_VIRIAL_22_VDW);
	GET_TENSOR(virial_22_repvdw, amd_reduction, REDUCTION_VIRIAL_22_REPVDW);
	GET_TENSOR(virial_22_elect, amd_reduction, REDUCTION_VIRIAL_22_ELECT);

	if (acg_intType == 0) {
		virial_cg = virial_22_vdw;
	} else if (acg_intType == 1) {
		virial_cg = virial_22_elect;
	} else if (acg_intType == 2) {
		virial_cg = (virial_22_vdw + virial_22_elect);
	} else if (acg_intType == 3) {
		Tensor virial_22_attvdw = virial_22_vdw - virial_22_repvdw;
		virial_cg = virial_22_attvdw;
	} else if (acg_intType == 4) {
		Tensor virial_22_attvdw = virial_22_vdw - virial_22_repvdw;
		virial_cg = (virial_22_attvdw + virial_22_elect);
	}

	virial_acg = scale[0] * virial_aa + scale[1] * virial_cg;
}
inline void Controller::acg_selectEnergy(BigReal &aaEnergy, BigReal &cgEnergy) {

	if (acg_intType == 0) {
		aaEnergy = aa_lj;
		cgEnergy = cg_lj;

	} else if (acg_intType == 1) {
		aaEnergy = aa_elect;
		cgEnergy = cg_elect;

	} else if (acg_intType == 2) {
		aaEnergy = aa_lj + aa_elect;
		cgEnergy = cg_lj + cg_elect;

	} else if (acg_intType == 3) {
		aaEnergy = aa_lj - aa_replj;
		cgEnergy = cg_lj - cg_replj;

	} else if (acg_intType == 4) {
		aaEnergy = aa_lj - aa_replj + aa_elect;
		cgEnergy = cg_lj - cg_replj + cg_elect;
	}
}
inline BigReal Controller::acg_lbd(int i) {
	//BigReal lbd_i = BigReal(i) * (acg_lbdMax - acg_lbdMin) / BigReal(acg_replicaNum - 1) + acg_lbdMin;
	//return (lbd_i);
	return (acg_lbdList[i]);
}
inline BigReal Controller::acg_lbdaa(BigReal lbd) {
	BigReal value;
	value = -lbd;
	return (value);
}
inline BigReal Controller::acg_lbdcg(BigReal lbd) {
	BigReal value;
	value = lbd;
	return (value);
}
inline BigReal Controller::acg_cv(BigReal aa_energy, BigReal cg_energy, BigReal lbd) {
	BigReal cv = 0.0;
	cv = acg_lbdaa(lbd) * aa_energy + acg_lbdcg(lbd) * cg_energy;
	return (cv);
}
BigReal Controller::acg_eff_cv(BigReal aa_energy, BigReal cg_energy, BigReal *logZ) {
	BigReal kbT = BOLTZMANN * (simParams->adaptCgTemp);

	BigReal *logExp = new BigReal[acg_replicaNum];
	BigReal max = -9999999999999.0;
	for (int i = 0; i < acg_replicaNum; ++i) {
		BigReal lbd_i = acg_lbd(i);
		BigReal cv_i = acg_cv(aa_energy, cg_energy, lbd_i);
		BigReal c_i = acg_replicaWeight[i];
		logExp[i] = -cv_i * DENSITYALAPHA / kbT - logZ[i] * DENSITYALAPHA + log(c_i);
		if (max < logExp[i])
			max = logExp[i];
	}
	BigReal sum = 0.0;
	for (int i = 0; i < acg_replicaNum; ++i) {
		sum += exp_cut(logExp[i] - max);
	}
	delete[] logExp;
	return (-(log(sum) + max) * kbT / DENSITYALAPHA);
}
Vector Controller::acg_eff_lbd(BigReal aa_energy, BigReal cg_energy, BigReal *logZ) {
	BigReal kbT = BOLTZMANN * (simParams->adaptCgTemp);

	BigReal *logExp = new BigReal[acg_replicaNum];
	BigReal max = -9999999999.0;
	Vector scale;
	for (int i = 0; i < acg_replicaNum; ++i) {
		BigReal lbd_i = acg_lbd(i);
		BigReal cv_i = acg_cv(aa_energy, cg_energy, lbd_i);
		BigReal c_i = acg_replicaWeight[i];
		logExp[i] = -cv_i * DENSITYALAPHA / kbT - logZ[i] * DENSITYALAPHA + log(c_i);
		if (max < logExp[i])
			max = logExp[i];
	}

	BigReal sum_up_aa = 0.0;
	BigReal sum_up_cg = 0.0;
	BigReal sum_down = 0.0;
	for (int i = 0; i < acg_replicaNum; ++i) {
		BigReal lbd_i = acg_lbd(i);
		BigReal preFi = exp_cut(logExp[i] - max);
		sum_up_aa += preFi * acg_lbdaa(lbd_i);
		sum_up_cg += preFi * acg_lbdcg(lbd_i);
		sum_down += preFi;
	}
	scale[0] = sum_up_aa / sum_down;
	scale[1] = sum_up_cg / sum_down;
	scale[2] = 0.0;
	delete[] logExp;
	return (scale);
}
void Controller::acg_estimate(BigReal *mean, BigReal *sgm, BigReal *logz) {

	BigReal kbT = BOLTZMANN * (simParams->adaptCgTemp);
	BigReal *logExp = new BigReal[acg_samplesNum];

	for (int i = 0; i < acg_replicaNum; ++i) {
		BigReal lbd_i = acg_lbd(i);

		//calc mean and sgm
		BigReal max = -9999999999.0;
		BigReal eff_cv_k = 0.0;
		BigReal cv_k = 0.0;
		for (int k = 0; k < acg_samplesNum; ++k) {
			eff_cv_k = acg_effCvSamples[k];
			cv_k = acg_cv(acg_aaEnergySamples[k], acg_cgEnergySamples[k], lbd_i);
			logExp[k] = -(cv_k - eff_cv_k) / kbT;
			if (max <= logExp[k])
				max = logExp[k];
		}

		BigReal sum_mean = 0.0;
		BigReal sum_sgm = 0.0;
		BigReal sum = 0.0;

		for (int k = 0; k < acg_samplesNum; ++k) {
			BigReal PreF = exp_cut(logExp[k] - max);
			cv_k = acg_cv(acg_aaEnergySamples[k], acg_cgEnergySamples[k], lbd_i);
			sum_mean += PreF * cv_k;
			sum_sgm += PreF * cv_k * cv_k;
			sum += PreF;
		}

		mean[i] = sum_mean / sum;
		BigReal sgm_ = sum_sgm / sum - mean[i] * mean[i];
		sgm[i] = sgm_ >= 0 ? sqrt(sgm_) : sqrt(-sgm_);
		logz[i] = log(sum) + max + log(1.0 / BigReal(acg_samplesNum));
		//done calc log_Z[i]
	}

	delete[] logExp;

}
void Controller::acg_eflatUpdateLogZ(int step, BigReal *logz, BigReal *weight) {
	//estimate weight
	int updateTimes = int(step / simParams->adaptCgUpdateFreq);
	BigReal est_logWeight[acg_replicaNum];
	BigReal wsum = 0.0;
	BigReal wmax = -99999999999.0;

	for (int i = 1; i < acg_replicaNum; ++i) {
		logz[i] -= logz[0];
	}
	logz[0] = 0.0;
	for (int i = 0; i < acg_replicaNum; ++i) {
		est_logWeight[i] = logz[i] - acg_logZ[i] + log(acg_replicaWeight[i]);
		if (wmax <= est_logWeight[i])
			wmax = est_logWeight[i];

	}
	for (int i = 0; i < acg_replicaNum; ++i) {
		wsum += exp_cut(est_logWeight[i] - wmax);
	}
	for (int i = 0; i < acg_replicaNum; ++i) {
		est_logWeight[i] -= log(wsum) + wmax;
		weight[i] = exp_cut(est_logWeight[i]);
		if ((updateTimes == 1) && !acg_inFileOpen) {
			acg_logZ[i] = logz[i];
		} else {
			BigReal dw = est_logWeight[i] - log(acg_replicaWeight[i]);
			BigReal WCAP = simParams->WCAP;
			if (dw > WCAP)
				dw = WCAP;
			else if (dw < -WCAP)
				dw = -WCAP;
			acg_logZ[i] += dw;
		}
	}
	for (int i = 1; i < acg_replicaNum; ++i) {
		acg_logZ[i] -= acg_logZ[0];
	}
	acg_logZ[0] = 0.0;

}
void Controller::acg_purifyLogZ(int step, BigReal *logZ) {
	BigReal kbT = BOLTZMANN * (simParams->adaptCgTemp);

	BigReal *eff_cv = new BigReal[acg_samplesNum];
	BigReal *logZ_ = new BigReal[acg_replicaNum];
	BigReal d_error_max = -9999999999999.0;

	int loop = 0;
	for (loop = 0; loop < MAXLOOP; ++loop) {
		for (int k = 0; k < acg_samplesNum; ++k) {
			BigReal aa_ener_k = acg_aaEnergySamples[k];
			BigReal cg_ener_k = acg_cgEnergySamples[k];
			eff_cv[k] = acg_eff_cv(aa_ener_k, cg_ener_k, logZ);
		}

		for (int i = 0; i < acg_replicaNum; ++i) {
			BigReal lbd_i = acg_lbd(i);

			//get <exp(-beta*(Q_i-Qeff))>_eff=Z_i
			BigReal *logExp = new BigReal[acg_samplesNum];
			BigReal max = -99999999999.0;
			for (int k = 0; k < acg_samplesNum; ++k) {
				BigReal aa_ener_k = acg_aaEnergySamples[k];
				BigReal cg_ener_k = acg_cgEnergySamples[k];
				BigReal eff_cv_k = eff_cv[k];
				BigReal cv_k = acg_cv(aa_ener_k, cg_ener_k, lbd_i);
				logExp[k] = -(cv_k - eff_cv_k) / kbT;
				if (max < logExp[k])
					max = logExp[k];
			}
			BigReal sum = 0.0;
			for (int k = 0; k < acg_samplesNum; ++k) {
				sum += exp_cut(logExp[k] - max);
			}
			delete[] logExp;
			logZ_[i] = log(sum) + max + log(1.0 / BigReal(acg_samplesNum));
		}

		d_error_max = -9999999999999.0;
		for (int i = 0; i < acg_replicaNum; ++i) {
			BigReal d_error = (logZ_[i] - logZ[i]) * (logZ_[i] - logZ[i]);
			if (d_error_max < d_error)
				d_error_max = d_error;
			logZ[i] = logZ_[i];
		}

		if (d_error_max <= LOGZERROR)
			break;
	}
	iout << "PURIFY: STEP: " << step << " STOP-AT: " << loop + 1 << "/" << MAXLOOP << "  MAX-ERROR: " << d_error_max
			<< "\n" << endi;
	delete[] logZ_;
	delete[] eff_cv;
}

#ifdef WHAMFLAG
void Controller::acg_cvEffHis(int step, BigReal **U_eff) {

	int current_l = (int(step / simParams->adaptCgUpdateFreq) - 1) % simParams->adaptCgWhamTimes;

	int updateTimes = int(step / simParams->adaptCgUpdateFreq);
	if (updateTimes >= simParams->adaptCgWhamTimes)
		updateTimes = simParams->adaptCgWhamTimes;

	int sampleNum = int(step / simParams->adaptCgSamplingFreq);
	if (sampleNum >= acg_samplesNum_wham)
		sampleNum = acg_samplesNum_wham;

	BigReal aa_energy_t = 0.0;
	BigReal cg_energy_t = 0.0;

	for (int t = 0; t < sampleNum; ++t) {
		aa_energy_t = acg_aaEnergySamples_wham[t];
		cg_energy_t = acg_cgEnergySamples_wham[t];
		U_eff[t][current_l] = acg_eff_cv(aa_energy_t, cg_energy_t, acg_logZHis[current_l]);
	}

	for (int l = 0; l < updateTimes; ++l) {
		for (int t = 0; t < acg_samplesNum; ++t) {
			aa_energy_t = acg_aaEnergySamples[t];
			cg_energy_t = acg_cgEnergySamples[t];
			int t_ = t + current_l * acg_samplesNum;
			U_eff[t_][l] = acg_eff_cv(aa_energy_t, cg_energy_t, acg_logZHis[l]);
		}
	}

}
BigReal Controller::acg_whamEstimateLogZEff(int step, BigReal **U_eff, BigReal *logZEffHis) {
	BigReal kbT = BOLTZMANN * (simParams->adaptCgTemp);

	int updateTimes = int(step / simParams->adaptCgUpdateFreq);
	if (updateTimes >= simParams->adaptCgWhamTimes)
		updateTimes = simParams->adaptCgWhamTimes;

	int sampleNum = int(step / simParams->adaptCgSamplingFreq / SKIP);
	if (sampleNum >= int(acg_samplesNum_wham / SKIP))
		sampleNum = int(acg_samplesNum_wham / SKIP);

	BigReal *logZEffHis_ = new BigReal[updateTimes];
	for (int l = 0; l < updateTimes; ++l) {
		logZEffHis_[l] = logZEffHis[l] - logZEffHis[0];
	}

	BigReal *logExp_k = new BigReal[sampleNum];
	BigReal *logExp_t = new BigReal[updateTimes];
	for (int k = 0; k < updateTimes; ++k) {

		BigReal max_k = -9999999999.0;
		for (int t = 0; t < sampleNum; ++t) {
			int t_ = t * SKIP;
			//get Sum(l=1 to n) exp(-beta*(Ueff_l-Ueff_k)-logZ_l)=exp(-logExp_t[l])
			BigReal max_t = -999999999999.0;
			for (int l = 0; l < updateTimes; ++l) {
				logExp_t[l] = -(U_eff[t_][l] - U_eff[t_][k]) / kbT - logZEffHis_[l];
				if (max_t <= logExp_t[l])
					max_t = logExp_t[l];
			}
			BigReal expSum_t = 0.0;
			for (int l = 0; l < updateTimes; ++l) {
				expSum_t += exp_cut(logExp_t[l] - max_t);
			}
			//finish get

			logExp_k[t] = -log(expSum_t) - max_t;
			if (max_k <= logExp_k[t])
				max_k = logExp_k[t];
		}

		BigReal sum_k = 0.0;
		for (int t = 0; t < sampleNum; ++t) {
			sum_k += exp_cut(logExp_k[t] - max_k);
		}

		logZEffHis[k] = log(sum_k) + max_k - log(BigReal(sampleNum) / BigReal(updateTimes));
		if (k >= 1)
			logZEffHis[k] -= logZEffHis[0];
	}
	logZEffHis[0] = 0.0;

	BigReal d_error_max = -9999999999.0;
	for (int l = 0; l < updateTimes; ++l) {
		BigReal d_error = sqrt((logZEffHis_[l] - logZEffHis[l]) * (logZEffHis_[l] - logZEffHis[l]));
		if (d_error_max < d_error)
			d_error_max = d_error;
	}

	delete[] logZEffHis_;
	delete[] logExp_k;
	delete[] logExp_t;
	return (d_error_max);

}
void Controller::acg_whamEstimate(int step, BigReal **U_eff, BigReal *est_logZ, BigReal *mean, BigReal *sgm) {
	BigReal kbT = BOLTZMANN * (simParams->adaptCgTemp);

	int updateTimes = int(step / simParams->adaptCgUpdateFreq);
	if (updateTimes >= simParams->adaptCgWhamTimes)
		updateTimes = simParams->adaptCgWhamTimes;

	int sampleNum = int(step / simParams->adaptCgSamplingFreq);
	if (sampleNum >= acg_samplesNum_wham)
		sampleNum = acg_samplesNum_wham;

	BigReal *logExp_i = new BigReal[sampleNum];
	BigReal *logExp_t = new BigReal[updateTimes];
	for (int i = 0; i < acg_replicaNum; ++i) {

		BigReal lbd_i = acg_lbd(i);
		BigReal max_i = -9999999999.0;
		BigReal aa_energy_t = 0.0;
		BigReal cg_energy_t = 0.0;
		BigReal cv_t = 0.0;
		for (int t = 0; t < sampleNum; ++t) {

			//get Sum(l=1 to n) exp(-beta*(Ueff_l-Ueff_k)-logZ_l)=exp(-logExp_t[l])
			BigReal max_t = -999999999999.0;
			aa_energy_t = acg_aaEnergySamples_wham[t];
			cg_energy_t = acg_cgEnergySamples_wham[t];
			cv_t = acg_cv(aa_energy_t, cg_energy_t, lbd_i);
			for (int l = 0; l < updateTimes; ++l) {
				logExp_t[l] = -(U_eff[t][l] - cv_t) / kbT - acg_logZEffHis[l];
				if (max_t <= logExp_t[l])
					max_t = logExp_t[l];
			}
			BigReal expSum_t = 0.0;
			for (int l = 0; l < updateTimes; ++l) {
				expSum_t += exp_cut(logExp_t[l] - max_t);
			}
			//finish get

			logExp_i[t] = -log(expSum_t) - max_t;
			if (max_i <= logExp_i[t])
				max_i = logExp_i[t];
		}

		BigReal sum_i = 0.0;
		BigReal sum_mean_i = 0.0;
		BigReal sum_sgm_i = 0.0;
		for (int t = 0; t < sampleNum; ++t) {
			aa_energy_t = acg_aaEnergySamples_wham[t];
			cg_energy_t = acg_cgEnergySamples_wham[t];
			cv_t = acg_cv(aa_energy_t, cg_energy_t, lbd_i);
			BigReal preF = exp_cut(logExp_i[t] - max_i);
			sum_mean_i += preF * cv_t;
			sum_sgm_i += preF * cv_t * cv_t;
			sum_i += preF;
		}
		est_logZ[i] = log(sum_i) + max_i - log(BigReal(sampleNum) / BigReal(updateTimes));
		mean[i] = sum_mean_i / sum_i;
		sgm[i] = sum_sgm_i / sum_i;
		sgm[i] -= mean[i] * mean[i];
		sgm[i] = (sgm[i] >= 0.0) ? sqrt(sgm[i]) : sqrt(-sgm[i]);
		if (i != 0) {
			est_logZ[i] -= est_logZ[0];
		}
	}
	est_logZ[0] = 0.0;
	delete[] logExp_i;
	delete[] logExp_t;
}
#endif

void Controller::acg_init(void) {

	acg_rescalingFactor[0] = 0.0;
	acg_rescalingFactor[1] = 0.0;
	acg_rescalingFactor[2] = 0.0;
	acg_replicaNum = simParams->adaptCgReplicaNum;
	acg_samplesNum = int(simParams->adaptCgUpdateFreq / simParams->adaptCgSamplingFreq);

	acg_replicaWeight = new BigReal[acg_replicaNum];
	acg_logZ = new BigReal[acg_replicaNum];
	acg_lbdList = new BigReal[acg_replicaNum];
#ifdef WHAMFLAG
	acg_samplesNum_wham = acg_samplesNum * simParams->adaptCgWhamTimes;
	REALARRAY2D(acg_cvEff, acg_samplesNum_wham, simParams->adaptCgWhamTimes);
	REALARRAY2D(acg_logZHis, simParams->adaptCgWhamTimes, acg_replicaNum);
	REALARRAY1D(acg_logZEffHis, simParams->adaptCgWhamTimes);
	REALARRAY1D(acg_whamLogZ, acg_replicaNum);
#endif
	acg_meanCV = new BigReal[acg_replicaNum];
	acg_sgmCV = new BigReal[acg_replicaNum];
	acg_lbdHistogram = new BigReal[acg_replicaNum];
	acg_dynamicWeight = new BigReal[acg_replicaNum];
	acg_dynamicWeightSum = new BigReal[acg_replicaNum];
	acg_inFile.open(simParams->adaptCgInFile, std::ios::in);
	acg_inFileOpen = acg_inFile.is_open();

	if (simParams->adaptCgUpdateOn) {
		BigReal kbT = BOLTZMANN * (simParams->adaptCgTemp);
		acg_updateLog.open(simParams->adaptCgUpdateLog, std::ios::out);
		acg_updateLog.close();
		acg_updateLog.open(simParams->adaptCgUpdateLog, std::ios::app);
#ifdef WHAMFLAG
		acg_whamLog.open("ACGWHAM.LOG", std::ios::out);
		acg_whamLog.close();
		acg_whamLog.open("ACGWHAM.LOG", std::ios::app);
#endif
		if (acg_inFileOpen) {
			int index;
			BigReal lbd;
			BigReal radio;
			for (int i = 0; i < acg_replicaNum; ++i) {
				acg_inFile >> index >> acg_lbdList[i] >> acg_dynamicWeightSum[i] >> acg_replicaWeight[i]
						>> acg_meanCV[i] >> acg_sgmCV[i] >> acg_lbdHistogram[i] >> acg_logZ[i] >> acg_logZ[i];
				acg_lbdHistogram[i] = 0;
			}
			//normalize the acg_replicaWeight
			BigReal sum = 0.0;
			for (int i = 0; i < acg_replicaNum; ++i) {
				sum += acg_replicaWeight[i];
			}
			for (int i = 0; i < acg_replicaNum; ++i) {
				acg_replicaWeight[i] /= sum;
				acg_whamLogZ[i] = acg_logZ[i];
			}
			acg_inFile >> acg_lbdMin >> acg_lbdMax;
			acg_inFile.close();
		} else {
			acg_lbdMax = simParams->adaptCgLbdMax;
			acg_lbdMin = simParams->adaptCgLbdMin;
			BigReal aa_energy = 0.0;
			BigReal cg_energy = 0.0;
			//acg_selectEnergy(aa_energy, cg_energy);

			for (int i = 0; i < acg_replicaNum; ++i) {
				//BigReal lbd_i = acg_lbd(i);
				BigReal lbd_i = BigReal(i) * (acg_lbdMax - acg_lbdMin) / BigReal(acg_replicaNum - 1) + acg_lbdMin;
				acg_lbdList[i] = lbd_i;
				BigReal cv_i = acg_cv(aa_energy, cg_energy, lbd_i);
				BigReal p_eta_i = 1.0 / BigReal(acg_replicaNum);
				acg_replicaWeight[i] = p_eta_i;
				acg_logZ[i] = -cv_i / kbT;
#ifdef WHAMFLAG
				acg_logZHis[0][i] = acg_logZ[i];
				acg_whamLogZ[i] = acg_logZ[i];
#endif
				acg_lbdHistogram[i] = 0;
				acg_meanCV[i] = 0.0;
				acg_sgmCV[i] = 0.0;
				acg_dynamicWeight[i] = 0.0;
				acg_dynamicWeightSum[i] = 1.0 / BigReal(acg_replicaNum);
			}
#ifdef WHAMFLAG
			acg_logZEffHis[0] = 0.0;
#endif
		}

		acg_aaEnergySamples = new BigReal[acg_samplesNum];
		acg_cgEnergySamples = new BigReal[acg_samplesNum];
		acg_effCvSamples = new BigReal[acg_samplesNum];
#ifdef WHAMFLAG
		acg_aaEnergySamples_wham = new BigReal[acg_samplesNum_wham];
		acg_cgEnergySamples_wham = new BigReal[acg_samplesNum_wham];
#endif

	} else {

		if (acg_inFileOpen) {
			int index;
			BigReal lbd;
			BigReal radio;
			for (int i = 0; i < acg_replicaNum; ++i) {
				acg_inFile >> index >> acg_lbdList[i] >> acg_dynamicWeightSum[i] >> acg_replicaWeight[i]
						>> acg_meanCV[i] >> acg_sgmCV[i] >> acg_lbdHistogram[i] >> acg_logZ[i] >> acg_logZ[i];
				acg_lbdHistogram[i] = 0;
			}

			//normalize the acg_replicaWeight
			BigReal sum = 0.0;
			for (int i = 0; i < acg_replicaNum; ++i) {
				sum += acg_replicaWeight[i];
			}
			for (int i = 0; i < acg_replicaNum; ++i) {
				acg_replicaWeight[i] /= sum;
			}

			acg_inFile >> acg_lbdMin >> acg_lbdMax;
			acg_inFile.close();
		} else {
			acg_lbdMax = simParams->adaptCgLbdMax;
			acg_lbdMin = simParams->adaptCgLbdMin;
			for (int i = 0; i < acg_replicaNum; ++i) {
				//BigReal lbd_i = acg_lbd(i);
				BigReal lbd_i = BigReal(i) * (acg_lbdMax - acg_lbdMin) / BigReal(acg_replicaNum - 1) + acg_lbdMin;
				BigReal p_eta_i = 1.0 / BigReal(acg_replicaNum);
				acg_replicaWeight[i] = p_eta_i;
				acg_logZ[i] = 0.0;
				acg_lbdHistogram[i] = 0;
				acg_meanCV[i] = 0.0;
				acg_sgmCV[i] = 0.0;
			}
		}

	}

}

//==============================================================
//
//			MODULE FOR ITS
//
//==============================================================
inline void Controller::its_updateVirial(Vector scale) {
	Tensor virial_dihe;
	Tensor virial_nbond;
	Tensor virial_slow;

	Tensor virial_11_vdw;
	Tensor virial_12_vdw;
	Tensor virial_22_vdw;

	Tensor virial_11_repvdw;
	Tensor virial_12_repvdw;
	Tensor virial_22_repvdw;

	Tensor virial_11_elect;
	Tensor virial_12_elect;
	Tensor virial_22_elect;

	GET_TENSOR(virial_dihe, amd_reduction, REDUCTION_VIRIAL_AMD_DIHE);
	GET_TENSOR(virial_nbond, amd_reduction, REDUCTION_VIRIAL_NBOND);
	GET_TENSOR(virial_slow, amd_reduction, REDUCTION_VIRIAL_SLOW);

	GET_TENSOR(virial_11_vdw, amd_reduction, REDUCTION_VIRIAL_11_VDW);
	GET_TENSOR(virial_12_vdw, amd_reduction, REDUCTION_VIRIAL_12_VDW);
	GET_TENSOR(virial_22_vdw, amd_reduction, REDUCTION_VIRIAL_22_VDW);

	GET_TENSOR(virial_11_repvdw, amd_reduction, REDUCTION_VIRIAL_11_REPVDW);
	GET_TENSOR(virial_12_repvdw, amd_reduction, REDUCTION_VIRIAL_12_REPVDW);
	GET_TENSOR(virial_22_repvdw, amd_reduction, REDUCTION_VIRIAL_22_REPVDW);

	GET_TENSOR(virial_11_elect, amd_reduction, REDUCTION_VIRIAL_11_ELECT);
	GET_TENSOR(virial_12_elect, amd_reduction, REDUCTION_VIRIAL_12_ELECT);
	GET_TENSOR(virial_22_elect, amd_reduction, REDUCTION_VIRIAL_22_ELECT);

	if (simParams->onlyIntOn) {
		if (its_intType == 0) {
			virial_its = scale[0] * virial_12_vdw;

		} else if (its_intType == 1) {
			virial_its = scale[0] * virial_12_elect;

		} else if (its_intType == 2) {
			virial_its = scale[0] * (virial_12_vdw + virial_12_elect);

		} else if (its_intType == 3) {
			Tensor virial_12_attvdw = virial_12_vdw - virial_12_repvdw;
			virial_its = scale[0] * virial_12_attvdw;

		} else if (its_intType == 4) {
			Tensor virial_12_attvdw = virial_12_vdw - virial_12_repvdw;
			virial_its = scale[0] * (virial_12_attvdw + virial_12_elect);
		}
	} else {
		if (its_intType == 0) {
			virial_its = scale[0] * (virial_11_vdw + virial_12_vdw + virial_22_vdw);

		} else if (its_intType == 1) {
			virial_its = scale[0] * (virial_11_elect + virial_12_elect + virial_22_elect);

		} else if (its_intType == 2) {
			virial_its = scale[0] * (virial_11_vdw + virial_12_vdw + virial_22_vdw);
			virial_its += scale[0] * (virial_11_elect + virial_12_elect + virial_22_elect);

		} else if (its_intType == 3) {
			Tensor virial_11_attvdw = virial_11_vdw - virial_11_repvdw;
			Tensor virial_12_attvdw = virial_12_vdw - virial_12_repvdw;
			Tensor virial_22_attvdw = virial_22_vdw - virial_22_repvdw;
			virial_its = scale[0] * (virial_11_attvdw + virial_12_attvdw + virial_22_attvdw);

		} else if (its_intType == 4) {
			Tensor virial_11_attvdw = virial_11_vdw - virial_11_repvdw;
			Tensor virial_12_attvdw = virial_12_vdw - virial_12_repvdw;
			Tensor virial_22_attvdw = virial_22_vdw - virial_22_repvdw;
			virial_its = scale[0] * (virial_11_attvdw + virial_12_attvdw + virial_22_attvdw);
			virial_its += scale[0] * (virial_11_elect + virial_12_elect + virial_22_elect);
		}
	}
	if (simParams->group1 == -1) {
		virial_its = scale[0] * (virial_dihe + virial_nbond + virial_slow);
	}
}
inline void Controller::its_selectEnergy(BigReal &energy) {
	if (its_intType == 0) {
		energy = its_lj;

	} else if (its_intType == 1) {
		energy = its_elect;

	} else if (its_intType == 2) {
		energy = its_lj + its_elect;

	} else if (its_intType == 3) {
		energy = its_lj - its_replj;

	} else if (its_intType == 4) {
		energy = its_lj - its_replj + its_elect;

	}
}
inline BigReal Controller::its_lbd(int i) {
	return (its_lbdList[i]);
}
inline BigReal Controller::its_cv(BigReal energy, BigReal lbd) {
	BigReal cv = 0.0;
	cv = lbd * energy;
	return (cv);
}
BigReal Controller::its_eff_cv(BigReal energy, BigReal *logZ) {
	BigReal kbT = BOLTZMANN * (simParams->itsTemp);

	BigReal *logExp = new BigReal[its_replicaNum];
	BigReal max = -9999999999999.0;
	for (int i = 0; i < its_replicaNum; ++i) {
		BigReal lbd_i = its_lbd(i);
		BigReal cv_i = its_cv(energy, lbd_i);
		BigReal c_i = its_replicaWeight[i];
		logExp[i] = -cv_i * DENSITYALAPHA / kbT - logZ[i] * DENSITYALAPHA + log(c_i);
		if (max < logExp[i])
			max = logExp[i];
	}
	BigReal sum = 0.0;
	for (int i = 0; i < its_replicaNum; ++i) {
		sum += exp_cut(logExp[i] - max);
	}
	delete[] logExp;
	return (-(log(sum) + max) * kbT / DENSITYALAPHA);
}
Vector Controller::its_eff_lbd(BigReal energy, BigReal *logZ) {
	BigReal kbT = BOLTZMANN * (simParams->itsTemp);

	BigReal *logExp = new BigReal[its_replicaNum];
	BigReal max = -9999999999.0;
	Vector scale;
	for (int i = 0; i < its_replicaNum; ++i) {
		BigReal lbd_i = its_lbd(i);
		BigReal cv_i = its_cv(energy, lbd_i);
		BigReal c_i = its_replicaWeight[i];
		logExp[i] = -cv_i * DENSITYALAPHA / kbT - logZ[i] * DENSITYALAPHA + log(c_i);
		if (max < logExp[i])
			max = logExp[i];
	}

	BigReal sum_up = 0.0;
	BigReal sum_down = 0.0;
	for (int i = 0; i < its_replicaNum; ++i) {
		BigReal lbd_i = its_lbd(i);
		BigReal preFi = exp_cut(logExp[i] - max);
		sum_up += preFi * lbd_i;
		sum_down += preFi;
	}
	scale[0] = sum_up / sum_down;
	scale[1] = 0.0;
	scale[2] = 0.0;
	delete[] logExp;
	return (scale);
}
void Controller::its_estimate(BigReal *mean, BigReal *sgm, BigReal *logz) {

	BigReal kbT = BOLTZMANN * (simParams->itsTemp);
	BigReal *logExp = new BigReal[its_samplesNum];

	for (int i = 0; i < its_replicaNum; ++i) {
		BigReal lbd_i = its_lbd(i);

		//calc mean and sgm
		BigReal max = -9999999999.0;
		BigReal eff_cv_k = 0.0;
		BigReal cv_k = 0.0;
		for (int k = 0; k < its_samplesNum; ++k) {
			eff_cv_k = its_effCvSamples[k];
			cv_k = its_cv(its_energySamples[k], lbd_i);
			logExp[k] = -(cv_k - eff_cv_k) / kbT;
			if (max <= logExp[k])
				max = logExp[k];
		}

		BigReal sum_mean = 0.0;
		BigReal sum_sgm = 0.0;
		BigReal sum = 0.0;

		for (int k = 0; k < its_samplesNum; ++k) {
			BigReal PreF = exp_cut(logExp[k] - max);
			cv_k = its_cv(its_energySamples[k], lbd_i);
			sum_mean += PreF * cv_k;
			sum_sgm += PreF * cv_k * cv_k;
			sum += PreF;
		}

		mean[i] = sum_mean / sum;
		BigReal sgm_ = sum_sgm / sum - mean[i] * mean[i];
		sgm[i] = sgm_ >= 0 ? sqrt(sgm_) : sqrt(-sgm_);
		logz[i] = log(sum) + max + log(1.0 / BigReal(its_samplesNum));
		//done calc log_Z[i]
	}

	delete[] logExp;

}
void Controller::its_eflatUpdateLogZ(int step, BigReal *logz, BigReal *weight) {
	//estimate weight
	BigReal est_logWeight[its_replicaNum];
	BigReal wsum = 0.0;
	BigReal wmax = -99999999999.0;

	for (int i = 1; i < its_replicaNum; ++i) {
		logz[i] -= logz[0];
	}
	logz[0] = 0.0;
	for (int i = 0; i < its_replicaNum; ++i) {
		est_logWeight[i] = logz[i] - its_logZ[i] + log(its_replicaWeight[i]);
		if (wmax <= est_logWeight[i])
			wmax = est_logWeight[i];

	}
	for (int i = 0; i < its_replicaNum; ++i) {
		wsum += exp_cut(est_logWeight[i] - wmax);
	}
	for (int i = 0; i < its_replicaNum; ++i) {
		est_logWeight[i] -= log(wsum) + wmax;
		weight[i] = exp_cut(est_logWeight[i]);
		BigReal dw = est_logWeight[i] - log(its_replicaWeight[i]);
		BigReal WCAP = simParams->WCAP;
		if (dw > WCAP)
			dw = WCAP;
		else if (dw < -WCAP)
			dw = -WCAP;
		its_logZ[i] += dw;
	}
	for (int i = 1; i < its_replicaNum; ++i) {
		its_logZ[i] -= its_logZ[0];
	}
	its_logZ[0] = 0.0;

}

#ifdef WHAMFLAG
void Controller::its_cvEffHis(int step, BigReal **U_eff) {

	int current_l = (int(step / simParams->itsUpdateFreq) - 1) % simParams->itsWhamTimes;

	int updateTimes = int(step / simParams->itsUpdateFreq);
	if (updateTimes >= simParams->itsWhamTimes)
		updateTimes = simParams->itsWhamTimes;

	int sampleNum = int(step / simParams->itsSamplingFreq);
	if (sampleNum >= its_samplesNum_wham)
		sampleNum = its_samplesNum_wham;

	BigReal energy_t = 0.0;

	for (int t = 0; t < sampleNum; ++t) {
		energy_t = its_energySamples_wham[t];
		U_eff[t][current_l] = its_eff_cv(energy_t, its_logZHis[current_l]);
	}

	for (int l = 0; l < updateTimes; ++l) {
		for (int t = 0; t < its_samplesNum; ++t) {
			energy_t = its_energySamples[t];
			int t_ = t + current_l * its_samplesNum;
			U_eff[t_][l] = its_eff_cv(energy_t, its_logZHis[l]);
		}
	}

}
BigReal Controller::its_whamEstimateLogZEff(int step, BigReal **U_eff, BigReal *logZEffHis) {
	BigReal kbT = BOLTZMANN * (simParams->itsTemp);

	int updateTimes = int(step / simParams->itsUpdateFreq);
	if (updateTimes >= simParams->itsWhamTimes)
		updateTimes = simParams->itsWhamTimes;

	int sampleNum = int(step / simParams->itsSamplingFreq / SKIP);
	if (sampleNum >= int(its_samplesNum_wham / SKIP))
		sampleNum = int(its_samplesNum_wham / SKIP);

	BigReal *logZEffHis_ = new BigReal[updateTimes];
	for (int l = 0; l < updateTimes; ++l) {
		logZEffHis_[l] = logZEffHis[l] - logZEffHis[0];
	}

	BigReal *logExp_k = new BigReal[sampleNum];
	BigReal *logExp_t = new BigReal[updateTimes];
	for (int k = 0; k < updateTimes; ++k) {

		BigReal max_k = -9999999999.0;
		for (int t = 0; t < sampleNum; ++t) {
			int t_ = t * SKIP;
			//get Sum(l=1 to n) exp(-beta*(Ueff_l-Ueff_k)-logZ_l)=exp(-logExp_t[l])
			BigReal max_t = -999999999999.0;
			for (int l = 0; l < updateTimes; ++l) {
				logExp_t[l] = -(U_eff[t_][l] - U_eff[t_][k]) / kbT - logZEffHis_[l];
				if (max_t <= logExp_t[l])
					max_t = logExp_t[l];
			}
			BigReal expSum_t = 0.0;
			for (int l = 0; l < updateTimes; ++l) {
				expSum_t += exp_cut(logExp_t[l] - max_t);
			}
			//finish get

			logExp_k[t] = -log(expSum_t) - max_t;
			if (max_k <= logExp_k[t])
				max_k = logExp_k[t];
		}

		BigReal sum_k = 0.0;
		for (int t = 0; t < sampleNum; ++t) {
			sum_k += exp_cut(logExp_k[t] - max_k);
		}

		logZEffHis[k] = log(sum_k) + max_k - log(BigReal(sampleNum) / BigReal(updateTimes));
		if (k >= 1)
			logZEffHis[k] -= logZEffHis[0];
	}
	logZEffHis[0] = 0.0;

	BigReal d_error_max = -9999999999.0;
	for (int l = 0; l < updateTimes; ++l) {
		BigReal d_error = sqrt((logZEffHis_[l] - logZEffHis[l]) * (logZEffHis_[l] - logZEffHis[l]));
		if (d_error_max < d_error)
			d_error_max = d_error;
	}

	delete[] logZEffHis_;
	delete[] logExp_k;
	delete[] logExp_t;
	return (d_error_max);

}
void Controller::its_whamEstimate(int step, BigReal **U_eff, BigReal *est_logZ, BigReal *mean, BigReal *sgm) {
	BigReal kbT = BOLTZMANN * (simParams->itsTemp);

	int updateTimes = int(step / simParams->itsUpdateFreq);
	if (updateTimes >= simParams->itsWhamTimes)
		updateTimes = simParams->itsWhamTimes;

	int sampleNum = int(step / simParams->itsSamplingFreq);
	if (sampleNum >= its_samplesNum_wham)
		sampleNum = its_samplesNum_wham;

	BigReal *logExp_i = new BigReal[sampleNum];
	BigReal *logExp_t = new BigReal[updateTimes];
	for (int i = 0; i < its_replicaNum; ++i) {

		BigReal lbd_i = its_lbd(i);
		BigReal max_i = -9999999999.0;
		BigReal energy_t = 0.0;
		BigReal cv_t = 0.0;
		for (int t = 0; t < sampleNum; ++t) {

			//get Sum(l=1 to n) exp(-beta*(Ueff_l-Ueff_k)-logZ_l)=exp(-logExp_t[l])
			BigReal max_t = -999999999999.0;
			energy_t = its_energySamples_wham[t];
			cv_t = its_cv(energy_t, lbd_i);
			for (int l = 0; l < updateTimes; ++l) {
				logExp_t[l] = -(U_eff[t][l] - cv_t) / kbT - its_logZEffHis[l];
				if (max_t <= logExp_t[l])
					max_t = logExp_t[l];
			}
			BigReal expSum_t = 0.0;
			for (int l = 0; l < updateTimes; ++l) {
				expSum_t += exp_cut(logExp_t[l] - max_t);
			}
			//finish get

			logExp_i[t] = -log(expSum_t) - max_t;
			if (max_i <= logExp_i[t])
				max_i = logExp_i[t];
		}

		BigReal sum_i = 0.0;
		BigReal sum_mean_i = 0.0;
		BigReal sum_sgm_i = 0.0;
		for (int t = 0; t < sampleNum; ++t) {
			energy_t = its_energySamples_wham[t];
			cv_t = its_cv(energy_t, lbd_i);
			BigReal preF = exp_cut(logExp_i[t] - max_i);
			sum_mean_i += preF * cv_t;
			sum_sgm_i += preF * cv_t * cv_t;
			sum_i += preF;
		}
		est_logZ[i] = log(sum_i) + max_i - log(BigReal(sampleNum) / BigReal(updateTimes));
		mean[i] = sum_mean_i / sum_i;
		sgm[i] = sum_sgm_i / sum_i;
		sgm[i] -= mean[i] * mean[i];
		sgm[i] = (sgm[i] >= 0.0) ? sqrt(sgm[i]) : sqrt(-sgm[i]);
		if (i != 0) {
			est_logZ[i] -= est_logZ[0];
		}
	}
	est_logZ[0] = 0.0;
	delete[] logExp_i;
	delete[] logExp_t;
}
#endif

void Controller::its_init(void) {

	its_rescalingFactor[0] = 0.0;
	its_rescalingFactor[1] = 0.0;
	its_rescalingFactor[2] = 0.0;
	its_replicaNum = simParams->itsReplicaNum;
	its_samplesNum = int(simParams->itsUpdateFreq / simParams->itsSamplingFreq);

	its_replicaWeight = new BigReal[its_replicaNum];
	its_logZ = new BigReal[its_replicaNum];
	its_lbdList = new BigReal[its_replicaNum];
#ifdef WHAMFLAG
	its_samplesNum_wham = its_samplesNum * simParams->itsWhamTimes;
	REALARRAY2D(its_cvEff, its_samplesNum_wham, simParams->itsWhamTimes);
	REALARRAY2D(its_logZHis, simParams->itsWhamTimes, its_replicaNum);
	REALARRAY1D(its_logZEffHis, simParams->itsWhamTimes);
	REALARRAY1D(its_whamLogZ, its_replicaNum);
#endif
	its_meanCV = new BigReal[its_replicaNum];
	its_sgmCV = new BigReal[its_replicaNum];
	its_lbdHistogram = new BigReal[its_replicaNum];
	its_dynamicWeight = new BigReal[its_replicaNum];
	its_dynamicWeightSum = new BigReal[its_replicaNum];
	its_inFile.open(simParams->itsInFile, std::ios::in);
	its_inFileOpen = its_inFile.is_open();

	if (simParams->itsUpdateOn) {
		BigReal kbT = BOLTZMANN * (simParams->itsTemp);
		its_updateLog.open(simParams->itsUpdateLog, std::ios::out);
		its_updateLog.close();
		its_updateLog.open(simParams->itsUpdateLog, std::ios::app);
#ifdef WHAMFLAG
		its_whamLog.open("ITSWHAM.LOG", std::ios::out);
		its_whamLog.close();
		its_whamLog.open("ITSWHAM.LOG", std::ios::app);
#endif
		if (its_inFileOpen) {
			int index;
			BigReal lbd;
			BigReal radio;
			for (int i = 0; i < its_replicaNum; ++i) {
				its_inFile >> index >> its_lbdList[i] >> its_dynamicWeightSum[i] >> its_replicaWeight[i]
						>> its_meanCV[i] >> its_sgmCV[i] >> its_lbdHistogram[i] >> its_logZ[i] >> its_logZ[i];
				its_lbdHistogram[i] = 0;
			}
			//normalize the its_replicaWeight
			BigReal sum = 0.0;
			for (int i = 0; i < its_replicaNum; ++i) {
				sum += its_replicaWeight[i];
			}
			for (int i = 0; i < its_replicaNum; ++i) {
				its_replicaWeight[i] /= sum;
				its_whamLogZ[i] = its_logZ[i];
			}
			its_inFile >> its_lbdMin >> its_lbdMax;
			its_inFile.close();
		} else {
			its_lbdMax = simParams->itsLbdMax;
			its_lbdMin = simParams->itsLbdMin;
			BigReal energy = 0.0;
			//its_selectEnergy(energy);

			for (int i = 0; i < its_replicaNum; ++i) {
				//BigReal lbd_i = its_lbd(i);
				BigReal lbd_i = BigReal(i) * (its_lbdMax - its_lbdMin) / BigReal(its_replicaNum - 1) + its_lbdMin;
				its_lbdList[i] = lbd_i;
				BigReal cv_i = its_cv(energy, lbd_i);
				BigReal p_eta_i = 1.0 / BigReal(its_replicaNum);
				its_replicaWeight[i] = p_eta_i;
				its_logZ[i] = -cv_i / kbT;
#ifdef WHAMFLAG
				its_logZHis[0][i] = its_logZ[i];
				its_whamLogZ[i] = its_logZ[i] - its_logZ[0];
#endif
				its_lbdHistogram[i] = 0;
				its_meanCV[i] = 0.0;
				its_sgmCV[i] = 0.0;
				its_dynamicWeight[i] = 0.0;
				its_dynamicWeightSum[i] = 1.0 / BigReal(its_replicaNum);
			}
#ifdef WHAMFLAG
			its_logZEffHis[0] = 0.0;
#endif
		}

		its_energySamples = new BigReal[its_samplesNum];
		its_effCvSamples = new BigReal[its_samplesNum];
#ifdef WHAMFLAG
		its_energySamples_wham = new BigReal[its_samplesNum_wham];
#endif

	} else {

		if (its_inFileOpen) {
			int index;
			BigReal lbd;
			BigReal radio;
			for (int i = 0; i < its_replicaNum; ++i) {
				its_inFile >> index >> its_lbdList[i] >> its_dynamicWeightSum[i] >> its_replicaWeight[i]
						>> its_meanCV[i] >> its_sgmCV[i] >> its_lbdHistogram[i] >> its_logZ[i] >> its_logZ[i];
				its_lbdHistogram[i] = 0;
			}

			//normalize the its_replicaWeight
			BigReal sum = 0.0;
			for (int i = 0; i < its_replicaNum; ++i) {
				sum += its_replicaWeight[i];
			}
			for (int i = 0; i < its_replicaNum; ++i) {
				its_replicaWeight[i] /= sum;
			}

			its_inFile >> its_lbdMin >> its_lbdMax;
			its_inFile.close();
		} else {
			its_lbdMax = simParams->itsLbdMax;
			its_lbdMin = simParams->itsLbdMin;
			for (int i = 0; i < its_replicaNum; ++i) {
				//BigReal lbd_i = its_lbd(i);
				BigReal lbd_i = BigReal(i) * (its_lbdMax - its_lbdMin) / BigReal(its_replicaNum - 1) + its_lbdMin;
				BigReal p_eta_i = 1.0 / BigReal(its_replicaNum);
				its_replicaWeight[i] = p_eta_i;
				its_logZ[i] = 0.0;
				its_lbdHistogram[i] = 0;
				its_meanCV[i] = 0.0;
				its_sgmCV[i] = 0.0;
			}
		}

	}

}

//
//update mean and sgm
//
/*
 int avgtimes = int(step / simParams->adaptCgUpdateFreq) - 1;
 BigReal eff_samplesNum = 0.0;
 for (int i = 0; i < acg_replicaNum; ++i) {
 eff_samplesNum += acg_dynamicWeight[i];
 }

 for (int i = 0; i < acg_replicaNum; ++i) {
 if ((avgtimes == 0) && !acg_inFileOpen) {
 //acg_logZ[i] = logz[i];
 acg_meanCV[i] = mean[i];
 acg_sgmCV[i] = sgm[i];
 }
 if ((avgtimes > 0) || acg_inFileOpen) {
 BigReal wi = acg_dynamicWeight[i] / eff_samplesNum;
 //acg_logZ[i] = acg_logZ[i] * acg_dynamicWeightSum[i] + logz[i] * wi;
 //acg_logZ[i] /= acg_dynamicWeightSum[i] + wi;
 acg_meanCV[i] = acg_meanCV[i] * acg_dynamicWeightSum[i] + mean[i] * wi;
 acg_meanCV[i] /= acg_dynamicWeightSum[i] + wi;
 acg_sgmCV[i] = acg_sgmCV[i] * acg_dynamicWeightSum[i] + sgm[i] * wi;
 acg_sgmCV[i] /= acg_dynamicWeightSum[i] + wi;
 acg_dynamicWeightSum[i] += wi;
 }
 }
 */

/*
 BigReal sgm2 = dlbd * simParams->adaptCgUpdateScaling;
 sgm2 = sgm2 * sgm2;

 for (int i = 0; i < acg_replicaNum; ++i) {
 BigReal lbd_i = acg_lbd(i);
 acg_dynamicWeight[i] += exp(-0.5 * (lbd_i - current_lbd) * (lbd_i - current_lbd) / sgm2);
 }
 */
