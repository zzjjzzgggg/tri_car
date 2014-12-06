/*
 * main.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#include <queue>
#include <condition_variable>
#include <chrono>

#include "stdafx.h"
#include "Queue.h"
#include "ExamMgr.h"


/**
 * count the number of triangles in int-pair NIdCnt.
 */
void sub_gt(const int id, Queue<int>& Que, ExamMgr& ExM, TIntPrV& NIdCnt) {
	TExeTm2 tm;
	int nid, ntrids;
	printf("[%d] is ready.\n", id);
	while (Que.TryPop(nid)) {
		ntrids = TSnap::GetNodeTriadsAll<PNEGraph>(ExM.GFull, nid);
		NIdCnt.Add(TIntPr(nid, ntrids));
	}
	printf("[%d] costs time: %s\n", id, tm.GetStr());
}


void multi_groundtruth(ExamMgr& ExM){
	TExeTm tm;
	// assign jobs
	Queue<int> Qu;
	for (ExamMgr::NI NI=ExM.GFull->BegNI(); NI<ExM.GFull->EndNI(); NI++) Qu.Push(NI.GetId());
	TIntPrV Vs[ExM.CPU];
	for (int n=0; n<ExM.CPU; n++) Vs[n] = TIntPrV();

	// assign threads
	std::vector<std::thread> threads;
	for (int n=0; n<ExM.CPU; n++) threads.emplace_back([n, &Qu, &ExM, &Vs] { sub_gt(n, Qu, ExM, Vs[n]); });
	for(std::thread& t: threads) t.join();

	// collect results
	PSOut tz = TZipOut::New(ExM.GetNTFNm());
	tz->PutStrLn(TStr::Fmt("# Nodes: %d", ExM.N));
	TIntH TriadCntH;
	for (int n=0; n<ExM.CPU; n++) {
		for (int i=0; i<Vs[n].Len(); i++) {
			TriadCntH(Vs[n][i].Val2) ++;
			tz->PutStrLn(TStr::Fmt("%d\t%d", Vs[n][i].Val1.Val, Vs[n][i].Val2.Val));
		}
	}

	TIntPrV TriadCntV;
	TriadCntH.GetKeyDatPrV(TriadCntV);
	TriadCntV.Sort();
	FILE* fw=fopen(ExM.GetGTFNm().CStr(), "w");
	fprintf(fw, "# Time cost: %.2f seconds\n", tm.GetSecs());
	fprintf(fw, "# Nodes: %d\n", ExM.N);
	for (int i=0; i<TriadCntV.Len(); i++) {
		const int card = TriadCntV[i].Val1, freq = TriadCntV[i].Val2;
		const double prob = freq/(double)ExM.N;
		fprintf(fw, "%d\t%d\t%.6e\n", card, freq, prob);
	}
	fclose(fw);
}


void UC_groundtruth(ExamMgr& ExM){
	TExeTm tm;
	TIntPrV gV;
	ExM.SampleUC(gV);
	gV.Sort();
	FILE* fw=fopen(ExM.GetGTFNm().CStr(), "w");
	fprintf(fw, "# Time cost: %.2f seconds\n", tm.GetSecs());
	fprintf(fw, "# Nodes: %d\n", ExM.N);
	for (int i=0; i<gV.Len(); i++) {
		const int card = gV[i].Val1, freq = gV[i].Val2;
		const double prob = freq/(double)ExM.N;
		fprintf(fw, "%d\t%d\t%.6e\n", card, freq, prob);
	}
	fclose(fw);
}


void samle_graph(ExamMgr& ExM){
	PNEGraph G = PNEGraph::TObj::New();
	ExM.GetSampledGraph(G);
	TSnap::SaveEdgeList<PNEGraph>(G, ExM.GetSGFNm());
}


void sub_ef(const int id, Queue<TIntFltKd>& Que, ExamMgr& ExM, TIntFltKdV& TmV) {
	TIntFltKd Item;
	TIntPrV tridCnt;
	PNEGraph G = PNEGraph::TObj::New();
	while (Que.TryPop(Item)) {
		printf("[%d] %g\n", id, Item.Dat.Val);
		ExM.GetSampledGraph(G, Item.Dat);
		auto start = std::chrono::steady_clock::now();
		TSnap::GetTriadParticipAll(G, tridCnt);
		auto end = std::chrono::steady_clock::now();
		double millsecs = std::chrono::duration<double, std::milli>(end-start).count();
		TmV.Add(TIntFltKd(Item.Key, millsecs));
	}
}

void multi_ef(ExamMgr& ExM, const TFltV& PeV){
	Queue<TIntFltKd> Qu;
	for (int key=0; key<PeV.Len(); key++)
		for (int rpt=0; rpt<ExM.Rpt; rpt++)
			Qu.Push(TIntFltKd(key, PeV[key]));
	TIntFltKdV TmVs[ExM.CPU];
	for (int n=0; n<ExM.CPU; n++) TmVs[n] = TIntFltKdV();

	// assign threads
	std::vector<std::thread> threads;
	for (int n=0; n<ExM.CPU; n++) threads.emplace_back([n, &Qu, &ExM, &TmVs] { sub_ef(n, Qu, ExM, TmVs[n]); });
	for(std::thread& t: threads) t.join();

	// collect results
	TIntFltH KeyTmH;
	for (int n=0; n<ExM.CPU; n++)
		for (int i=0; i<TmVs[n].Len(); i++)
			KeyTmH(TmVs[n][i].Key) += TmVs[n][i].Dat;
	TFltPrV AvgV; TInt Key; TFlt SumTm;
	for (int kid = KeyTmH.FFirstKeyId(); KeyTmH.FNextKeyId(kid);) {
		KeyTmH.GetKeyDat(kid, Key, SumTm);
		AvgV.Add(TFltPr(PeV[Key], SumTm/ExM.Rpt/1000));
	}
	BIO::SaveFltPrV(AvgV, ExM.GetEfFNm(), "%.2f\t%.2f");
}


int main(int argc, char* argv[]){
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	const TStr GFNm = Env.GetIfArgPrefixStr("-i:", "test.graph", "Input graph");
	const TStr FGFNm = Env.GetIfArgPrefixStr("-f:", "fg.graph", "Follower graph");
	const int Rpt = Env.GetIfArgPrefixInt("-r:", 10, "Repeat times");
	const int CPU = Env.GetIfArgPrefixInt("-cpu:", std::thread::hardware_concurrency(), "# of CPUs");
	const double Pe = Env.GetIfArgPrefixFlt("-p:", 0.1, "Edge sampling rate");
	TFltV PeDefV(6);
	for (int i=0; i<PeDefV.Len(); i++) PeDefV[i] = 0.05*(i+1);
	const TFltV PeV = Env.GetIfArgPrefixFltV("-ps:", PeDefV, "Edge sampling rates");
	const TStr Fmts = Env.GetIfArgPrefixStr("-c:", "", "What to compute:"
				"\n\tg: get groundtruth (multi-core)"
				"\n\tc: get UC groundtruth (simple)"
				"\n\ts: get sampled graph"
				"\n\te: compare efficiency");
	if (Env.IsEndOfRun()) return 0;
	ExamMgr ExM;
	TExeTm2 tm;
	if (Fmts.SearchCh('g') != -1) {
		ExM.SetActionGraph(GFNm).SetCPU(CPU);
		multi_groundtruth(ExM);
		printf("Saved to\n  %s\n  %s\n", ExM.GetGTFNm().CStr(), ExM.GetNTFNm().CStr());
	} else if (Fmts.SearchCh('e') != -1) {
		ExM.SetActionGraph(GFNm).SetCPU(CPU).SetRepeat(Rpt);
		multi_ef(ExM, PeV);
		printf("Saved to %s\n", ExM.GetEfFNm().CStr());
	} else if (Fmts.SearchCh('c') != -1) {
		ExM.SetSocialGraph(FGFNm).SetActionGraph(GFNm).SetPEdge(1).SetPSocial(1);
		UC_groundtruth(ExM);
		printf("Saved to %s\n", ExM.GetGTFNm().CStr());
	} else if (Fmts.SearchCh('s') != -1) {
		ExM.SetActionGraph(GFNm).SetPEdge(Pe);
		samle_graph(ExM);
		printf("Saved to %s\n", ExM.GetSGFNm().CStr());
	}
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

