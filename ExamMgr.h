/*
 * ExamMgr.h
 *
 *  Created on: Sep 20, 2014
 *      Author: jzzhao
 */

#ifndef EXAMMGR_H_
#define EXAMMGR_H_

#include "stdafx.h"

class ExamMgr {
public:
	typedef PNEGraph::TObj::TEdgeI EI;
	typedef PNEGraph::TObj::TNodeI NI;
	TStr GFNm;
	PNEGraph GFull;
	int N, W, CPU, Rpt;
	double PEdge;
	bool TrimTail;
public:
	ExamMgr(const TStr& GFnm, const int NCPU=8, const int MX_TC=10000, const double Pe=0.1, const int NRpt=12, const bool Tail=false);
	void GetSampledGraph(PNEGraph& G, const double Pe=-1);
	void Sample(TIntPrV& gV, const double Pe=-1);
	TStr GetTHFNm() const{
		return GFNm.GetFPath() + TStr::Fmt("th_%s_W%dK_p%g.dist", GFNm.GetFMid().CStr(), W/1000, PEdge);
	}
	TStr GetNTHFNm() const{
		return GFNm.GetFPath() + TStr::Fmt("nth_%s_W%dK_p%g.dist", GFNm.GetFMid().CStr(), W/1000, PEdge);
	}
	TStr GetGTFNm() const {
		return GFNm.GetFPath() + TStr::Fmt("gndtruth_%s.dat", GFNm.GetFMid().CStr());
	}
	TStr GetNTFNm() const {
		return GFNm.GetFPath() + TStr::Fmt("ndtriads_%s.gz", GFNm.GetFMid().CStr());
	}
	int GetRpt() const {return CPU*Rpt;}
};

#endif /* EXAMMGR_H_ */
