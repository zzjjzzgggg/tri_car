/*
 * ExamMgr.cpp
 *
 *  Created on: Sep 20, 2014
 *      Author: jzzhao
 */

#include "ExamMgr.h"

ExamMgr::ExamMgr(const TStr& GFnm, const double PEdge) {
	// TODO Auto-generated constructor stub
	this->PEdge = PEdge;
	GFull = TSnap::LoadEdgeList<PNEGraph>(GFnm);
}

void ExamMgr::GetSampledGraph(PNEGraph& G){
	TRnd rnd;
	for(EI ei=GFull->BegEI(); ei!=GFull->EndEI(); ei++){
		if(rnd.GetUniDev() > PEdge) continue;
		const int src = ei.GetSrcNId(), dst = ei.GetDstNId();
		if(!G->IsNode(src)) G->AddNode(src);
		if(!G->IsNode(dst)) G->AddNode(dst);
		G->AddEdge(src, dst);
	}
	printf("Sampled: \n\t nodes: %d, edges: %d\n", G->GetNodes(), G->GetEdges());
}
