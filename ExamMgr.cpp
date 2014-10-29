/*
 * ExamMgr.cpp
 *
 *  Created on: Sep 20, 2014
 *      Author: jzzhao
 */

#include "ExamMgr.h"

ExamMgr::ExamMgr(const TStr& GFnm, const int NCPU, const int MX_TC, const double Pe, const int NRpt) {
	// TODO Auto-generated constructor stub
	GFNm = GFnm;
	GFull = TSnap::LoadEdgeList<PNEGraph>(GFnm);
	N = GFull->GetNodes();
	W = MX_TC;
	PEdge = Pe;
	CPU = NCPU;
	Rpt = NRpt;
}

void ExamMgr::GetSampledGraph(PNEGraph& G, const double Pe){
	const double p = Pe<0 ? PEdge : Pe;
	TRnd rnd;
	G->Clr();
	for(EI ei=GFull->BegEI(); ei!=GFull->EndEI(); ei++){
		if(rnd.GetUniDev() > p) continue;
		const int src = ei.GetSrcNId(), dst = ei.GetDstNId();
		if(!G->IsNode(src)) G->AddNode(src);
		if(!G->IsNode(dst)) G->AddNode(dst);
		G->AddEdge(src, dst);
	}
	G->Defrag();
}
