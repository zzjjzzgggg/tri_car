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
	TFltV AV, ZV, ThV_pre;
	TIntPrV gV;
public:
	int M, Iters;
	TFltV ThV;
public:
	TCEMGeneral(const int W, const double Pdelta, const TIntPrV& igV): W(W), Pd(Pdelta), Iters(0) {
		M = g = 0;
		for (int i=0; i<igV.Len(); i++) {
			const int card = igV[i].Val1, freq = igV[i].Val2;
			if(card>0){
				gV.Add(TIntPr(card, freq));
				g += freq;
				if (card > M) M = card;
			}
		}
		// init A
		AV.Gen((M+1)*(2*W-M+2)/2);
		for (int id=0; id<gV.Len(); id++){
			int j = gV[id].Val1;
			for (int i=j; i<=W; i++)
				AV[Idx(i,j)] = (TSpecFunc::Binomial(j, i, Pd)/(1-TSpecFunc::Binomial(0, i, Pd)));
		}
		// init Theta
		ThV.Gen(W+1); ThV_pre.Gen(W+1);
		TRandom::InitUni(ThV, 1);
		// space for Z
		ZV.Gen((M+1)*(2*W-M+2)/2);
	};
	bool Run(const int max_iter = 200);
private:
	int Idx(const int i, const int j) const { return (i<=M) ? (i*(i+1)/2+j) : ((M+1)*(M+2)/2+(i-M-1)*(M+1)+j); }
	void EStep();
	bool MStep(const double Eps=0.005);
	void Scale();
	void ScaleTail();
};

#endif /* TCEMGENERAL_H_ */
