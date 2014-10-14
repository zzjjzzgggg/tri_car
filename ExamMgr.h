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
	PNEGraph GFull;
	int N, W;
	double PEdge;
public:
	ExamMgr(const TStr& GFnm, const int MX_TC=1000, const double Pe=0.1);
	void GetSampledGraph(PNEGraph& G, const double Pe=-1);
};

#endif /* EXAMMGR_H_ */
