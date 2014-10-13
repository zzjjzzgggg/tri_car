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
	int W, N, g, M;
	double p;  /// p == p_delta
	TFltV BV, ZV, ThV_pre;
	TIntH gH;
public:
	TFltV ThV;
public:
	TCEM(const int W, const int N, const double p_delta, const TIntPrV& TridCnt): W(W), N(N), p(p_delta) {
		M = g = 0;
		for (int i=0; i<TridCnt.Len(); i++) {
			const int card = TridCnt[i].Val1, freq = TridCnt[i].Val2;
			gH(card) = freq;
			g += freq;
			if (card > M) M = card;
		}
		g -= gH(0);
		gH(0) = N-g;
		// init B
		BV.Gen((W+1)*(M+1));
		for(int j=0; j<=M; j++)
			for (int i=j; i<=W; i++)
				BV[Idx(i,j)] = TSpecFunc::Binomial(j, i, p);
		// init Theta
		ThV.Gen(W+1); ThV_pre.Gen(W+1);
		TRandom::InitUni(ThV);
		// space for Z
		ZV.Gen((W+1)*(M+1));
	};
	bool Run(const int max_iter = 1000);
private:
	const int Idx(const int i, const int j){ return i*M+j; }
	void EStep();
	bool MStep(const double Eps=0.005);
};

#endif /* TCEM_H_ */
