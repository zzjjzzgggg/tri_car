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
private:
	typedef PNEGraph::TObj::TEdgeI EI;
	PNEGraph GFull;
	double PEdge;
public:
	ExamMgr(const TStr& GFnm, const double PEdge);
	void GetSampledGraph(PNEGraph& G);
};

#endif /* EXAMMGR_H_ */
