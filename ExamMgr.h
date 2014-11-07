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
	PNEGraph GFull;
	PBNEGraph UCGFull;
	PNGraph FGFull;
	TStr GFNm, FGFNm;
	int N, W, CPU, Rpt;
	double PEdge, PRelation;
	bool TrimTail;
private:
	void CheckSocialRelation(const TIntPrV& Users, PNEGraph& G);
public:
	ExamMgr(const TStr& GFnm, const int NCPU=8, const int MX_TC=10000,
			const double Pe=0.1, const int NRpt=12, const bool Tail=false): GFNm(GFnm),
			FGFNm(""), W(MX_TC), CPU(NCPU), Rpt(NRpt), PEdge(Pe), PRelation(0), TrimTail(Tail){
		GFull = TSnap::LoadEdgeList<PNEGraph>(GFnm);
		N = GFull->GetNodes();
	}
	/**
	 * NC is the number of OSN content
	 */
	ExamMgr(const TStr& UCGFnm, const TStr& FGFnm, const int NCPU=8, const int MX_TC=10000,
			const double Pe=0.1, const double Pr=0.1, const int NRpt=12, const bool Tail=false): GFNm(UCGFnm),
			FGFNm(FGFnm), W(MX_TC), CPU(NCPU), Rpt(NRpt), PEdge(Pe), PRelation(Pr), TrimTail(Tail){
		UCGFull = TSnap::LoadBNEGraph(UCGFnm);
		FGFull = TSnap::LoadEdgeList<PNGraph>(FGFnm);
		N = UCGFull->GetDstNodes();
	}
	void GetSampledGraph(PNEGraph& G, const double Pe=-1);

	void Sample(TIntPrV& gV, const double Pe=-1);
	void SampleUC(TIntPrV& gV);

	TStr GetTHFNm() const{
		return GFNm.GetFPath() + TStr::Fmt("th_%s_W%dK_p%g.dist", GFNm.GetFMid().CStr(), W/1000, PEdge);
	}
	TStr GetNTHFNm() const{
		return GFNm.GetFPath() + TStr::Fmt("nth_%s_W%dK_p%g.dist", GFNm.GetFMid().CStr(), W/1000, PEdge);
	}
	TStr GetGTFNm() const {
		return GFNm.GetFPath() + TStr::Fmt("gndtruth_%s.dat", GFNm.GetFMid().CStr());
	}
	TStr GetNTFNm() const {
		return GFNm.GetFPath() + TStr::Fmt("ndtriads_%s.gz", GFNm.GetFMid().CStr());
	}
	int GetRpt() const { return CPU*Rpt; }
	double GetPdUU() const { return pow(PEdge, 3); }
	double GetPdUC() const { return pow(PEdge, 2)*PRelation; }
};

#endif /* EXAMMGR_H_ */
