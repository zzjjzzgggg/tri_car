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
	for (int id=0; id<gV.Len(); id++){
		const int j = gV[id].Val1, freq = gV[id].Val2;
		double norm = 0;
		for (int i=j; i<=W; i++){
			const int k = Idx(i,j);
			ZV[k] = ThV[i]*TSpecFunc::Binomial(j, i, p);
			norm += ZV[k];
			ZV[k] *= freq;
		}
		for (int i=j; i<=W; i++){
			const int k = Idx(i,j);
			ZV[k] = norm<1E-10 ? 0 : ZV[k]/norm;
		}
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
	for (int id=0; id<gV.Len(); id++){
		const int j = gV[id].Val1;
		for(int i=j; i<=W; i++){
			const int k = Idx(i,j);
			ThV[i] += ZV[k];
			norm += ZV[k];
		}
	}
	// normalize Thv and calculate diff
	double diff = 0;
	for (int i=0; i<=W; i++) {
		ThV[i] /= norm;
		diff += TMath::Abs(ThV[i] - ThV_pre[i]);
	}
	return diff <= Eps;
}

/**
 * do EM iterations
 */
bool TCEM::Run(const int max_iter){
	for (Iters=0; Iters<max_iter; Iters++){
		EStep();
		if(MStep()) return true;
	}
	return false;
}

