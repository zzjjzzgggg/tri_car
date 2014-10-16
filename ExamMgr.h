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
public:
	ExamMgr(const TStr& GFnm, const int NCPU=8, const int MX_TC=10000, const double Pe=0.1, const int NRpt=12);
	void GetSampledGraph(PNEGraph& G, const double Pe=-1);
	TStr GetTHFNm() const{
		return GFNm.GetFPath() + TStr::Fmt("th_%s_W%dK_p%g.dist", GFNm.GetFMid().CStr(), W/1000, PEdge);
	}
	TStr GetNTHFNm() const{
		return GFNm.GetFPath() + TStr::Fmt("nth_%s_W%dK_p%g.dist", GFNm.GetFMid().CStr(), W/1000, PEdge);
	}
};

#endif /* EXAMMGR_H_ */
