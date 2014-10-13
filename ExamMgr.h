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
	int N;
public:
	ExamMgr(const TStr& GFnm);
	void GetSampledGraph(PNEGraph& G, const double PEdge);
};

#endif /* EXAMMGR_H_ */
