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
	TFltV ZV, ThV_pre;
	TIntH gH;
public:
	int M;
	TFltV ThV;
public:
	TCEMGeneral(const int W, const double Pdelta, const TIntH& igH): W(W), Pd(Pdelta) {
		TInt card, freq;
		M = g = 0;
		for (int i=0; i<igH.Len(); i++) {
			igH.GetKeyDat(i, card, freq);
			Assert(card<=W);
			if(card>0){
				gH(card) = freq;
				g += freq;
				if (card > M) M = card;
			}
		}
		// init Theta
		ThV.Gen(W+1); ThV_pre.Gen(W+1);
		TRandom::InitUni(ThV, 1);
		// space for Z
		ZV.Gen((M+1)*(2*W-M+2)/2);
	};
	bool Run(const int max_iter = 500);
private:
	int Idx(const int i, const int j) const { return (i<=M) ? (i*(i+1)/2+j) : ((M+1)*(M+2)/2+(i-M-1)*(M+1)+j); }
	void EStep();
	bool MStep(const double Eps=0.005);
	void Scale();
	void ScaleTail();
};

#endif /* TCEMGENERAL_H_ */
