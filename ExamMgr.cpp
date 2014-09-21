/*
 * ExamMgr.cpp
 *
 *  Created on: Sep 20, 2014
 *      Author: jzzhao
 */

#include "ExamMgr.h"

ExamMgr::ExamMgr(const TStr& GFnm) {
	// TODO Auto-generated constructor stub
	GFull = TSnap::LoadEdgeList<PNEGraph>(GFnm);
}

void ExamMgr::GetSampledGraph(PNEGraph& G, const double PEdge){
	TRnd rnd;
	G->Clr();
	for(EI ei=GFull->BegEI(); ei!=GFull->EndEI(); ei++){
		if(rnd.GetUniDev() > PEdge) continue;
		const int src = ei.GetSrcNId(), dst = ei.GetDstNId();
		if(!G->IsNode(src)) G->AddNode(src);
		if(!G->IsNode(dst)) G->AddNode(dst);
		G->AddEdge(src, dst);
	}
	printf("Sampled: \n\t nodes: %d, edges: %d\n", G->GetNodes(), G->GetEdges());
}
