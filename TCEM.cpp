/*
 * TCEM.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#include "TCEM.h"

void TCEM::Init(){
	// init B
	BV.Gen((W+1)*(M+1));
	for(int j=0; j<=M; j++){
		for (int i=j; i<=W; i++){
			BV[Idx(i,j)] = TSpecFunc::Binormal(j, i, p);
			IAssertR(BV[Idx(i,j)]>=0, TStr::Fmt("ER: %.6e", BV[Idx(i,j)].Val));
		}
	}
	// init Theta
	ThV.Gen(W+1); ThV_pre.Gen(W+1);
	TRandom::InitUni(ThV);
	// space for Z
	ZV.Gen((W+1)*(M+1));
}

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
void TCEM::MStep(){
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
	// normalize Thv
	for (int i=0; i<=W; i++) ThV[i] /= norm;
//	printf("N: %.2f  ", norm);
}

/**
 * do EM iterations
 */
bool TCEM::Run(const int max_iter){
	Init();
	for (int iter=0; iter<max_iter; iter++){
		EStep();
		MStep();
		if (IsConverged()) return true;
	}
	return false;
}

void TCEM::Save(TStr ofnm){
	FILE* fw = fopen(ofnm.CStr(), "w");
	for (int i=0; i<=W; i++) fprintf(fw, "%d\t%.6e\n", i, ThV[i].Val);
	fclose(fw);
}

