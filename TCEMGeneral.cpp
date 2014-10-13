/*
 * TCEM.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#include "TCEMGeneral.h"

/**
 * update [z_ij]
 */
void TCEMGeneral::EStep(){
	// estimate g0
	double norm;
	for (int j=1; j<=M; j++){
		if (gH(j)==0) continue;
		norm = 0;
		for (int i=j; i<=W; i++){
			ZV[Idx(i,j)] = PV[Idx(i,j)]*ThV[i];
			norm += ZV[Idx(i,j)];
			ZV[Idx(i,j)] *= gH(j);
		}
		for (int i=j; i<=W; i++)
			ZV[Idx(i,j)] = norm<1E-9 ? 0 : ZV[Idx(i,j)]/norm;
	}
	Assert(ZV[Idx(0,0)]>=0 && ZV[Idx(0,0)]<=1E-9);
}

/**
 * update ThV and Nt
 */
bool TCEMGeneral::MStep(const double Eps){
	// store the old parameters before we update them
	for (int i=1; i<=W; i++) {
		ThV_pre[i] = ThV[i];
		ThV[i] = 0;
	}
	// now update
	double norm = 0;
	for (int j=1; j<=M; j++){
		for(int i=j; i<=W; i++){
			ThV[i] += ZV[Idx(i,j)];
			norm += ZV[Idx(i,j)];
		}
	}
	// normalize Thv and calculate diff
	double diff = 0; //TMath::Abs(N-N_pre);
	for (int i=1; i<=W; i++) {
		ThV[i] /= norm;
		diff += TMath::Abs(ThV[i]-ThV_pre[i]);
	}
//	printf("%.6f\n", diff);
	return diff <= Eps;
}

/**
 * do EM iterations
 */
bool TCEMGeneral::Run(const int max_iter){
	for (int iter=0; iter<max_iter; iter++){
		EStep();
		if(MStep()) {
			Scale();
			return true;
		}
	}
	return false;
}

void TCEMGeneral::Scale(){
	ThV[0] = 0;
	for (int i=1; i<=W; i++) {
		ThV[i] /= (1-TSpecFunc::Binomial(0, i, Pd));
		ThV[0] += ThV[i];
	}
	double qth = 0; // q_theta
	for (int i=1; i<=W; i++) {
		ThV[i] /= ThV[0];
		qth += ThV[i]*TSpecFunc::Binomial(0, i, Pd);
	}
	N = g/(1-qth);
	ThV[0] = N; // store N to ThV[0]
//	printf("N =  %.2f, q = %.2f\n", N, qth);
}