/*
 * TCEM.h
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#ifndef TCEMGENERAL_H_
#define TCEMGENERAL_H_

#include "stdafx.h"

/**
 * For the case N is not given
 */
class TCEMGeneral {
private:
	int W, g;
	double Pd;  /// p == p_delta
	TFltV PV, ZV, ThV_pre;
	TIntH gH;
public:
	int M;
	TFltV ThV;
public:
	TCEMGeneral(const int W, const double Pdelta, const TIntPrV& TridCnt): W(W), Pd(Pdelta) {
		M = g = 0;
		for (int i=0; i<TridCnt.Len(); i++) {
			const int card = TridCnt[i].Val1, freq = TridCnt[i].Val2;
			if (card>0){
				gH(card) = freq;
				if (card > M) M = card;
				g += freq;
			}
		}
		// init A
		PV.Gen((W+1)*(M+1));
		for (int j=1; j<=M; j++)
			for (int i=j; i<=W; i++)
				PV[Idx(i,j)] = TSpecFunc::Binomial(j, i, Pd)/(1-TSpecFunc::Binomial(0, i, Pd));
		// init Theta
		ThV.Gen(W+1); ThV_pre.Gen(W+1);
		TRandom::InitUni(ThV, 1);
		// space for Z
		ZV.Gen((W+1)*(M+1));
	};
	bool Run(const int max_iter = 500);
private:
	const int Idx(const int i, const int j){ return i*M+j; }
	void EStep();
	bool MStep(const double Eps=0.002);
	void Scale();
};

#endif /* TCEMGENERAL_H_ */
