/*
 * TCEM.h
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#ifndef TCEM_BETABINOM_GENERAL_H_
#define TCEM_BETABINOM_GENERAL_H_

#include "stdafx.h"

/**
 * For the case N is not given
 */
class TCEMBetaBinomGeneral {
private:
	int W, M, g, BoundW;
	double Pd;
	TFltV ZV, ThV_pre;
	TIntIntFltTrV NonZeroZV;
	TIntPrV gV;
public:
	double alpha;
	TFltV ThV;
private:
	int Idx(const int i, const int j) const { return (i<=M) ? (i*(i+1)/2+j) : ((M+1)*(M+2)/2+(i-M-1)*(M+1)+j); }
	double GetA(const int i, const int j) const { return TSpecFunc::BetaBinomial(j, i, Pd/alpha, (1-Pd)/alpha) / (1-TSpecFunc::BetaBinomial(0, i, Pd/alpha, (1-Pd)/alpha)); }
	double GetQ(const int i) const { return TSpecFunc::BetaBinomial(0, i, Pd/alpha, (1-Pd)/alpha); }
	void EStep();
	bool MStep_theta(const double Eps=2e-3);
	bool MStep_alpha(const double Eps=1e-6, const int MxNewtonIters = 100);
	void Scale();
public:
	TCEMBetaBinomGeneral(const int W, const int BW, const double Pd, const double alpha, const TIntPrV& igV):
		W(W), BoundW(BW), Pd(Pd), alpha(alpha) {
		M = g = 0;
		for (int id=0; id<igV.Len(); id++) {
			const int card = igV[id].Val1, freq = igV[id].Val2;
			if(card>=1 && card<=W) {
				gV.Add(TIntPr(card, freq));
				g += freq;
				if (card > M) M = card;
			}
		}

		// init ThV
		ThV.Gen(W+1); ThV_pre.Gen(W+1);
		TRandom::InitUni(ThV, 1);

		// space for Z
		ZV.Gen((M+1)*(2*W-M+2)/2);
	};
	bool Run(const int MxEMIters = 500);
};

#endif /* TCEM_H_ */
