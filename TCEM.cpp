/*
 * TCEM.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#include "TCEM.h"

/**
 * update [z_ij]
 */
void TCEM::EStep(){
	for (int j=0; j<=M; j++){
		if (gH(j)==0) continue;
		double norm = 0;
		for (int i=j; i<=W; i++){
			ZV[Idx(i,j)] = BV[Idx(i,j)]*ThV[i];
			norm += ZV[Idx(i,j)];
			ZV[Idx(i,j)] *= gH(j);
		}
		for (int i=j; i<=W; i++)
			ZV[Idx(i,j)] = norm<1E-9 ? 0 : ZV[Idx(i,j)]/norm;
	}
}

/**
 * update ThV and Nt
 */
bool TCEM::MStep(const double Eps){
	// store the old parameters before we update them
	for (int i=0; i<=W; i++) {
		ThV_pre[i] = ThV[i];
		ThV[i] = 0;
	}
	// now update
	double norm = 0;
	for (int j=0; j<=M; j++)
		for(int i=j; i<=W; i++){
			ThV[i] += ZV[Idx(i,j)];
			norm += ZV[Idx(i,j)];
		}
	// normalize Thv and calculate diff
	double diff = 0;
	for (int i=0; i<=W; i++) {
		ThV[i] /= norm;
		diff += TMath::Abs(ThV[i]-ThV_pre[i]);
	}
	return diff <= Eps;
}

/**
 * do EM iterations
 */
bool TCEM::Run(const int max_iter){
	for (int iter=0; iter<max_iter; iter++){
		EStep();
		if(MStep()) return true;
	}
	return false;
}

