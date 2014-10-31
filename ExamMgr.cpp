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


void ExamMgr::Sample(TIntH& gH, const double Pe){
	const double p = Pe<0 ? PEdge : Pe;
	PNEGraph G = PNEGraph::TObj::New();
	TRnd rnd;
	for(EI ei=GFull->BegEI(); ei!=GFull->EndEI(); ei++){
		if(rnd.GetUniDev() > p) continue;
		const int src = ei.GetSrcNId(), dst = ei.GetDstNId();
		if(!G->IsNode(src)) G->AddNode(src);
		if(!G->IsNode(dst)) G->AddNode(dst);
		G->AddEdge(src, dst);
	}
	TIntPrV TridCnt;
	TSnap::GetTriadParticipAll(G, TridCnt);
	gH.Clr();
	for (int i=0; i<TridCnt.Len(); i++) {
		const int card = TridCnt[i].Val1, freq = TridCnt[i].Val2;
		if (card<0||card>W) continue;
		gH(card) = freq;
	}
}
