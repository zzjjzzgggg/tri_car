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
	int N, W, CPU, Rpt, BoundW;
	double PEdge, PSocial, alpha;
	bool TrimTail;
private:
	int CountTrids(const TIntPrV& Users);
public:
	ExamMgr(): GFNm(""), FGFNm(""), N(0), W(10000), CPU(1), Rpt(12), BoundW(100),
		PEdge(0.1), PSocial(0.1), alpha(0.0001),
		TrimTail(false) {}

	// method chaining
	// for user-content interaction, set social graph first.
	ExamMgr& SetActionGraph(const TStr& GFNm) {
		this->GFNm = GFNm;
		if (FGFNm.Len()==0) { // if social graph is not set
			GFull = TSnap::LoadEdgeList<PNEGraph>(GFNm);
			N = GFull->GetNodes();
		} else {
			UCGFull = TSnap::LoadBNEGraph(GFNm);
			N = UCGFull->GetDstNodes();
		}
		return *this;
	}
	ExamMgr& SetSocialGraph(const TStr& GFNm) { FGFNm=GFNm; FGFull = TSnap::LoadEdgeList<PNGraph>(FGFNm); return *this; }
	ExamMgr& SetCPU(const int& nCPU) { this->CPU = nCPU; return *this; }
	ExamMgr& SetW(const int& W) { this->W = W; return *this; }
	ExamMgr& SetRepeat(const int& Rpt) { this->Rpt = Rpt; return *this; }
	ExamMgr& SetPEdge(const double& Pe) { this->PEdge = Pe; return *this; }
	ExamMgr& SetPSocial(const double& Pr) { this->PSocial = Pr; return *this; }
	ExamMgr& IsTrimTail(const bool& TrimTail) { this->TrimTail = TrimTail; return *this; }
	ExamMgr& SetBoundW(const int& BW) { this->BoundW = BW; return *this; }
	ExamMgr& SetAlpha(const double& alpha) { this->alpha = alpha; return *this; }

	void GetSampledGraph(PNEGraph& G, const double Pe=-1);
	void Sample(TIntPrV& gV, const double Pe=-1);
	void SampleUC(TIntPrV& gV);

	void TrimTailTh(TFltV& ThV);
	void TrimTailNTh(TFltV& ThV, const double Alpha);

	TStr GetTHFNm() const { return GFNm.GetFPath()+TStr::Fmt("th_%s_W%dK_p%g.dist", GFNm.GetFMid().CStr(), W/1000, PEdge); }
	TStr GetBTHFNm() const { return GFNm.GetFPath()+TStr::Fmt("bth_%s_W%dK_p%g.dist", GFNm.GetFMid().CStr(), W/1000, PEdge); }
	TStr GetNTHFNm() const { return GFNm.GetFPath()+TStr::Fmt("nth_%s_W%dK_p%g.dist", GFNm.GetFMid().CStr(), W/1000, PEdge); }
	TStr GetBNTHFNm() const { return GFNm.GetFPath()+TStr::Fmt("bnth_%s_W%dK_p%g.dist", GFNm.GetFMid().CStr(), W/1000, PEdge); }
	TStr GetGTFNm() const { return GFNm.GetFPath()+TStr::Fmt("gndtruth_%s.dat", GFNm.GetFMid().CStr()); }
	TStr GetNTFNm() const { return GFNm.GetFPath()+TStr::Fmt("ndtrids_%s.gz", GFNm.GetFMid().CStr()); }
	TStr GetSGFNm() const { return GFNm.GetFPath()+TStr::Fmt("%s_p%g.gz", GFNm.GetFMid().CStr(), PEdge); }
	TStr GetEfFNm() const { return GFNm.GetFPath()+TStr::Fmt("eff_%s.dat", GFNm.GetFMid().CStr()); }
	int GetRpt() const { return CPU*Rpt; }
	double GetPdUU() const { return pow(PEdge, 3); }
	double GetPdUC() const { return pow(PEdge, 2)*PSocial; }
};

#endif /* EXAMMGR_H_ */
