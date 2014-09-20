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
	int W, M, N, g;
	double p;  /// p == p_delta
	TFltV BV, ZV, ThV_pre;
	TIntH gH;
public:
	TFltV ThV;
public:
	TCEM(const int W, const int M, const int N, const double p_delta, const TIntH& freqH): W(W), M(M), N(N), p(p_delta) {
		g = 0;
		for (int i=0; i<freqH.Len(); i++) {
			gH(freqH.GetKey(i)) = freqH[i];
			g += freqH[i];
		}
		g -= gH(0);
		gH(0) = N-g;
	};
	bool Run(const int max_iter = 1000);
	void Save(TStr ofnm);
private:
	const int Idx(const int i, const int j){ return i*M+j; }
	bool IsConverged(const double eps=0.001) {
		double diff = 0;
		for(int i=0; i<=W; i++) diff += TMath::Abs(ThV[i]-ThV_pre[i]);
		double sum = 0;
		for(int i=0; i<=W; i++) sum += ThV[i];
		IAssertR(TMath::Abs(sum-1)<1E-9, TStr::Fmt("UN: %.6e", sum));
		return diff<=eps;
	}

	void Init();
	void EStep();
	void MStep();
};

#endif /* TCEM_H_ */
