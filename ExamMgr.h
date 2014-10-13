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
	ExamMgr(const TStr& GFnm, const int MX_TC, const double Pe);
	void GetSampledGraph(PNEGraph& G);
};

#endif /* EXAMMGR_H_ */
