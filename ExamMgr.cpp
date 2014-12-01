/*
 * ExamMgr.cpp
 *
 *  Created on: Sep 20, 2014
 *      Author: jzzhao
 */

#include "ExamMgr.h"

void ExamMgr::GetSampledGraph(PNEGraph& G, const double Pe){
	const double p = Pe<0 ? PEdge : Pe;
	TRnd rnd;
	G->Clr();
	for(EI ei=GFull->BegEI(); ei!=GFull->EndEI(); ei++){
		if(rnd.GetUniDev() <= p) {
			const int src = ei.GetSrcNId(), dst = ei.GetDstNId();
			if(!G->IsNode(src)) G->AddNode(src);
			if(!G->IsNode(dst)) G->AddNode(dst);
			G->AddEdge(src, dst);
		}
	}
	G->Defrag();
}

void ExamMgr::Sample(TIntPrV& gV, const double Pe){
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
	TSnap::GetTriadParticipAll<PNEGraph>(G, TridCnt);
	gV.Clr();
	for (int i=0; i<TridCnt.Len(); i++) {
		const int card = TridCnt[i].Val1, freq = TridCnt[i].Val2;
		if (card>=0 && card<=W) gV.Add(TIntPr(card, freq));
	}
}

void ExamMgr::SampleUC(TIntPrV& gV){
	PNEGraph G = PNEGraph::New();
	PBNEGraph UCGSampled = PBNEGraph::TObj::New();
	TRnd rnd;
	// obtain sampled graph, i.e., UCG and partial G
	for(TBNEGraph::TEdgeI ei=UCGFull->BegEI(); ei!=UCGFull->EndEI(); ei++){
		if(rnd.GetUniDev() > PEdge) continue;
		const int u = ei.GetSrcNId(), c = ei.GetDstNId();
		if(!G->IsNode(u)) G->AddNode(u);
		if(!G->IsNode(c)) G->AddNode(c);
		G->AddEdge(u, c);
		if(!UCGSampled->IsSrcNode(u)) UCGSampled->AddSrcNode(u);
		if(!UCGSampled->IsDstNode(c)) UCGSampled->AddDstNode(c);
		UCGSampled->AddEdge(ei);
	}
	TIntPrV TmUsrV;
	for(TBNEGraph::TNodeI ni = UCGSampled->BegDstNI(); ni < UCGSampled->EndDstNI(); ni++){
		if (ni.GetDeg()<2) continue;
		for (int d=0; d<ni.GetDeg(); d++) {
			TBNEGraph::TEdgeI ei = UCGSampled->GetEI(ni.GetEId(d));
			TmUsrV.AddSorted(TIntPr(ei.GetDat(), ei.GetSrcNId()));
		}
		CheckSocialRelation(TmUsrV, G);
		TmUsrV.Clr();
	}
	TIntH TridCntH;
	for(TBNEGraph::TNodeI ni = UCGSampled->BegDstNI(); ni < UCGSampled->EndDstNI(); ni++){
		const int card = TSnap::GetNodeTriadsAll<PNEGraph>(G, ni.GetId());
		TridCntH(card) ++;
	}
	gV.Clr();
	for (int id = TridCntH.FFirstKeyId(); TridCntH.FNextKeyId(id); ){
		const int card = TridCntH.GetKey(id), freq = TridCntH[id];
		gV.Add(TIntPr(card, freq));
	}
}


void ExamMgr::CheckSocialRelation(const TIntPrV& TmUsrs, PNEGraph& G){
	TRnd rnd;
	const int L = TmUsrs.Len();
	for (int i=0; i<L; i++){
		const int u = TmUsrs[i].Val2;
		for (int j=i+1; j<L; j++){
			const int v = TmUsrs[j].Val2;
			// check whether u <--- v with probability PRelation
			if (u!=v && rnd.GetUniDev()<=PRelation && FGFull->IsEdge(v, u))
				G->AddEdge(v, u);
		}
	}
}
