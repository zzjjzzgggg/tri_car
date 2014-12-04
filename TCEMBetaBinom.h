/*
 * TCEM.h
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#ifndef TCEM_BETABINOM_H_
#define TCEM_BETABINOM_H_

#include "stdafx.h"

/**
 * For the case N is given, and P(Y=j|X=i) is modeled by a beta binomial distribution.
 */
class TCEMBetaBinom {
private:
	int W, M, N, g, BoundW;
	double Pd;
	TFltV ZV, ThV_pre;
	TIntIntFltTrV NonZeroZV;
	TIntPrV gV;
public:
	double alpha;
	TFltV ThV;
private:
	int Idx(const int i, const int j) const { return (i<=M) ? (i*(i+1)/2+j) : ((M+1)*(M+2)/2+(i-M-1)*(M+1)+j); }
	void EStep();
	bool MStep_theta(const double Eps=0.005);
	bool MStep_alpha(const double Eps=0.0001, const int MxNewtonIters = 100);
public:
	TCEMBetaBinom(const int W, const int N, const int BW, const double Pd, const double alpha, const TIntPrV& igV):
		W(W), N(N), BoundW(BW), Pd(Pd), alpha(alpha) {
		M = g = 0;
		for (int i=0; i<igV.Len(); i++) {
			const int card = igV[i].Val1, freq = igV[i].Val2;
			if(card>=1 && card<=W){
				gV.Add(TIntPr(card, freq));
				g += freq;
				if (card > M) M = card;
			}
		}
		gV.Add(TIntPr(0, N-g));

		// init ThV
		ThV.Gen(W+1);ThV_pre.Gen(W+1);
//		TRandom::InitUni(ThV);

		TRnd Rnd;
		double norm = 0;
		for (int i=0; i<=W; i++) {
			ThV[i] = pow(i+1, -(Rnd.GetUniDev()*0.4+0.8));
			norm += ThV[i];
		}
		for (int i=0; i<=W; i++) ThV[i] /= norm;

		// space for Z
		ZV.Gen((M+1)*(2*W-M+2)/2);
	};

	bool Run(const int MxEMIters = 500);
};

#endif /* TCEM_H_ */
