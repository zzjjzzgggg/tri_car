/*
 * TCEM.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#include "TCEMBetaBinomGeneral.h"

/**
 * update [z_ij]
 */
void TCEMBetaBinomGeneral::EStep(){
	TInt i, j, k, gj;
	double norm;
	NonZeroZV.Clr();
	for (int id=0; id<gV.Len(); id++){
		gV[id].GetVal(j, gj);
		norm = 0;
		for (i=j; i<=W; i++){
			k = Idx(i,j);
			ZV[k] = ThV[i]*GetA(i,j);
			norm += ZV[k];
		}
		for (i=j; i<=W; i++){
			k = Idx(i,j);
			ZV[k] = norm<1E-9 ? 0 : gj*ZV[k]/norm;
			if (i<=BoundW && ZV[k]>1E-9) NonZeroZV.Add(TIntIntFltTr(i,j,ZV[k]));
		}
	}
}

/**
 * update ThV and Nt
 */
bool TCEMBetaBinomGeneral::MStep_theta(const double Eps){
	// store the old parameters before we update them
	for (int i=1; i<=W; i++) {
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
	for (int i=1; i<=W; i++) {
		ThV[i] /= norm;
		diff += fabs(ThV[i] - ThV_pre[i]);
	}
	return diff <= Eps;
}

/**
 * to speed up the estimation of alpha, we do not need handle nodes with i>MinW.
 */
bool TCEMBetaBinomGeneral::MStep_alpha(const double Eps, const int MxNewtonIters) {
	TInt i, j, cnt=0;
	TFlt dQ, ddQ, e1, e2, e, z, a, da, q, dq, ddq, alpha_pre = alpha+1;
	// Newton iterations
	while (fabs(alpha_pre - alpha) > Eps && cnt < MxNewtonIters) {
		dQ = ddQ = 0;
		for (int id=0; id < NonZeroZV.Len(); id++) {
			NonZeroZV[id].GetVal(i, j, z);
			e1 = e2 = a = da = 0; q = 1;
			for (int s = 0; s < i; s++) {
				if (s < j) e = s/(s*alpha+Pd);
				else e = (s-j)/((s-j)*alpha+1-Pd);
				e1 += e - s/(s*alpha+1);
				e2 += -e*e + pow(s/(s*alpha+1), 2);
				q *= 1 - Pd/(s*alpha+1);
				a += s*Pd / (s*alpha+1) / (s*alpha+1-Pd);
				da += -s*s*Pd*(2*s*alpha+2-Pd) / pow((s*alpha+1)*(s*alpha+1-Pd), 2);
			}
			dq = a*q;
			ddq = da*q + a*dq;
			dQ += z * ( e1 + dq/(1-q) );
			ddQ += z * ( e2 + pow(dq/(1-q), 2) + ddq/(1-q) );
		}
		alpha_pre = alpha;
		alpha -= dQ/ddQ;
		cnt++;
	}
	return cnt<MxNewtonIters;
}

/**
 * EM iterations
 */
bool TCEMBetaBinomGeneral::Run(const int MxEMIters){
	for (int iters=0; iters<MxEMIters; iters++){
		EStep();
		if (MStep_theta()) break;
	}
	for (int iters=0; iters<MxEMIters; iters++){
		EStep();
		if (MStep_alpha() && MStep_theta()) {
			Scale();
			return true;
		}
	}
	return false;
}

void TCEMBetaBinomGeneral::Scale(){
	ThV[0] = 0;
	for (int i=1; i<=W; i++) {
		ThV[i] /= 1- GetQ(i);
		ThV[0] += ThV[i];
	}
	double qth = 0; // q_theta
	for (int i=1; i<=W; i++) {
		ThV[i] /= ThV[0];
		qth += ThV[i]*GetQ(i);
	}
	ThV[0] = g/(1-qth);// store N to ThV[0]
}
