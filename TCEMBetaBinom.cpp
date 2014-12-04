/*
 * TCEM.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#include "TCEMBetaBinom.h"

/**
 * update [z_ij]
 */
void TCEMBetaBinom::EStep(){
	TInt j, gj;
	NonZeroZV.Clr();
	for (int id=0; id<gV.Len(); id++){
		gV[id].GetVal(j, gj);
		double norm = 0;
		for (int i=j; i<=W; i++){
			const int k = Idx(i,j);
			ZV[k] = ThV[i]*TSpecFunc::BetaBinomial(j, i, Pd/alpha, (1-Pd)/alpha);
			norm += ZV[k];
			ZV[k] *= gj;
		}
		for (int i=j; i<=W; i++){
			const int k = Idx(i,j);
			ZV[k] = norm<1E-9 ? 0 : ZV[k]/norm;
			if (i<=BoundW && ZV[k]>1E-9) NonZeroZV.Add(TIntIntFltTr(i,j,ZV[k]));
		}
	}
}

/**
 * update ThV and Nt
 */
bool TCEMBetaBinom::MStep_theta(const double Eps){
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
		diff += fabs(ThV[i] - ThV_pre[i]);
	}
	return diff <= Eps;
}

/**
 * to speed up the estimation of alpha, we do not need handle nodes with i>MinW.
 */
bool TCEMBetaBinom::MStep_alpha(const double Eps, const int MxNewtonIters) {
	TInt i, j, cnt=0;
	TFlt dQ, ddQ, z, e1, e2, e, alpha_pre = alpha+1;
	// Newton iterations
	while (fabs(alpha_pre - alpha) > Eps && cnt < MxNewtonIters) {
		dQ = ddQ = 0;
		for (int id=0; id < NonZeroZV.Len(); id++){
			NonZeroZV[id].GetVal(i, j, z);
			e1 = e2 = 0;
			for (int s = 0; s < i; s++) {
				if (s < j) e = s/(s*alpha+Pd);
				else e = (s-j)/((s-j)*alpha+1-Pd);
				e1 += e - s/(s*alpha+1);
				e2 += -e*e + pow(s/(s*alpha+1), 2);
			}
			dQ += z*e1;
			ddQ += z*e2;
		}
		alpha_pre = alpha;
		alpha -= dQ/ddQ;
		cnt++;
	}
	return cnt<MxNewtonIters;
}

/**
 * do EM iterations
 */
bool TCEMBetaBinom::Run(const int MxEMIters){
	for (int iters=0; iters<MxEMIters; iters++){
		EStep();
		if (MStep_theta()) break;
	}
	for (int iters=0; iters<MxEMIters; iters++){
		EStep();
		if (MStep_alpha() && MStep_theta()) return true;
	}
	return false;
}

