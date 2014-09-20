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
	TCEM(const int W, const int N, const double p_delta, const TIntH& freqH): W(W), N(N), p(p_delta) {
		M = g = 0;
		for (int i=0; i<freqH.Len(); i++) {
			gH(freqH.GetKey(i)) = freqH[i];
			g += freqH[i];
			if (freqH.GetKey(i)>M) M = freqH.GetKey(i);
		}
		g -= gH(0);
		gH(0) = N-g;
		// init B
		BV.Gen((W+1)*(M+1));
		for(int j=0; j<=M; j++){
			for (int i=j; i<=W; i++){
				BV[Idx(i,j)] = TSpecFunc::Binomial(j, i, p);
				IAssertR(BV[Idx(i,j)]>=0, TStr::Fmt("ER: %.6e", BV[Idx(i,j)].Val));
			}
		}
		// init Theta
		ThV.Gen(W+1); ThV_pre.Gen(W+1);
		// space for Z
		ZV.Gen((W+1)*(M+1));
	};
	bool Run(const int max_iter = 1000);
private:
	const int Idx(const int i, const int j){ return i*M+j; }
	void Init();
	void EStep();
	bool MStep(const double Eps=0.005);
};

#endif /* TCEM_H_ */
