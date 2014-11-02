/*
 * TCEM.h
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#ifndef TCEM_H_
#define TCEM_H_

#include "stdafx.h"

/**
 * For the case N is given
 */
class TCEM {
private:
	int W, N, g;
	double Pd;  /// p == p_delta
	TFltV BV, ZV, ThV_pre;
	TIntPrV gV;
public:
	int M, Iters;
	TFltV ThV;
public:
	TCEM(const int W, const int N, const double p_delta, const TIntPrV& igV): W(W), N(N), Pd(p_delta), Iters(0) {
		M = g = 0;
		for (int i=0; i<igV.Len(); i++) {
			const int card = igV[i].Val1, freq = igV[i].Val2;
			if(card>0){
				gV.Add(TIntPr(card, freq));
				g += freq;
				if (card > M) M = card;
			}
		}
		gV.Add(TIntPr(0, N-g));
		// init B
		BV.Gen((M+1)*(2*W-M+2)/2);
		for (int id=0; id<gV.Len(); id++){
			int j = gV[id].Val1;
			for (int i=j; i<=W; i++)
				BV[Idx(i,j)] = TSpecFunc::Binomial(j, i, Pd);
		}
		// init Theta
		ThV.Gen(W+1); ThV_pre.Gen(W+1);
		TRandom::InitUni(ThV);
		// space for Z
		ZV.Gen((M+1)*(2*W-M+2)/2);
	};
	bool Run(const int max_iter = 500);
private:
	int Idx(const int i, const int j) const { return (i<=M) ? (i*(i+1)/2+j) : ((M+1)*(M+2)/2+(i-M-1)*(M+1)+j); }
	void EStep();
	bool MStep(const double Eps=0.005);
};

#endif /* TCEM_H_ */
