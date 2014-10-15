/*
 * main.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#include <thread>
#include <mutex>
#include <vector>

#include "stdafx.h"
#include "TCEM.h"
#include "TCEMGeneral.h"
#include "ExamMgr.h"

using namespace std;

/**
 * count triangles per week
 */
void count_triangles(){
	PNEGraph G = PNEGraph::TObj::New();
	TStr Root = "/media/WinE/workspace/triadic_cardinality/enron/";
	int pre_week = 0, ntrid = 0;
	TIntPrV TridCnt; TIntV NTridsPerWeek;
	TSsParser Ss(Root+"enron_add_week_no.gz");
	while(Ss.Next()){
		const int week = Ss.GetInt(0), fid = Ss.GetInt(1), tid = Ss.GetInt(2);
		if(week!=pre_week && G->GetEdges()!=0) {
			TSnap::GetTriadParticipAll(G, TridCnt);
			for (int i=0; i<TridCnt.Len(); i++) ntrid += TridCnt[i].Val1*TridCnt[i].Val2;
			NTridsPerWeek.Add(ntrid/3);
			G->Clr(); TridCnt.Clr(); pre_week = week; ntrid = 0;
		}
		if (!G->IsNode(fid)) G->AddNode(fid);
		if (!G->IsNode(tid)) G->AddNode(tid);
		G->AddEdge(fid, tid);
	}
	BIO::SaveIntsWithIdx(NTridsPerWeek, Root+"NTridsPerWeek.dat");
}

/**
 * triadic cardinality distribution
 */
void dist_triangles(const int WK){
	PNEGraph G = PNEGraph::TObj::New();
	TStr Root = "/media/WinE/workspace/triadic_cardinality/enron/";
	TSsParser Ss(Root+"enron_add_week_no.gz");
	while(Ss.Next()){
		const int week = Ss.GetInt(0), fid = Ss.GetInt(1), tid = Ss.GetInt(2);
		if(week==WK) {
			if (!G->IsNode(fid)) G->AddNode(fid);
			if (!G->IsNode(tid)) G->AddNode(tid);
			G->AddEdge(fid, tid);
		}else if(week>WK) break;
	}
	TIntPrV TridCnt;
	TSnap::GetTriadParticipAll(G, TridCnt);
	int ntrid =0;
	for (int i=0; i<TridCnt.Len(); i++) ntrid += TridCnt[i].Val1*TridCnt[i].Val2;
	ntrid /= 3;
	BIO::SaveIntPrV(TridCnt, Root+TStr::Fmt("TridCnt_%d.dat", WK));
}

/**
 * degree distribution
 */
void dist_degree(const int WK){
	PNEGraph G = PNEGraph::TObj::New();
	TStr Root = "/media/WinE/workspace/triadic_cardinality/enron/";
	TSsParser Ss(Root+"enron_add_week_no.gz");
	while(Ss.Next()){
		const int week = Ss.GetInt(0), fid = Ss.GetInt(1), tid = Ss.GetInt(2);
		if(week==WK) {
			if (!G->IsNode(fid)) G->AddNode(fid);
			if (!G->IsNode(tid)) G->AddNode(tid);
			G->AddEdge(fid, tid);
		}else if(week>WK) break;
	}
	TIntPrV DegV;
	TSnap::GetDegCnt<PNEGraph>(G, DegV);
	BIO::SaveIntPrV(DegV, Root+TStr::Fmt("DegCnt_%d.dat", WK));
}

void stat_trids(const TStr& GFNm){
	PNEGraph G = TSnap::LoadEdgeList<PNEGraph>(GFNm);
	TIntPrV TridCnt;
	for(PNEGraph::TObj::TNodeI NI=G->BegNI(); NI<G->EndNI(); NI++){
		int nid = NI.GetId();
		int ntrids = TSnap::GetNodeTriadsAll(G, nid);
		TridCnt.Add(TIntPr(nid, ntrids));
	}
	BIO::SaveIntPrV(TridCnt, "NodeNTrids.dat");
}

int count_trids(PNEGraph& G){
	int ntrids = 0;
	for(PNEGraph::TObj::TNodeI NI=G->BegNI(); NI<G->EndNI(); NI++){
		ntrids += TSnap::GetNodeTriadsAll(G, NI.GetId());
	}
	return ntrids/3;
}

void gen_groundtruth(const TStr& GFNm){
	PNEGraph G = TSnap::LoadEdgeList<PNEGraph>(GFNm);
	double N = G->GetNodes();
	TIntPrV tridCnt;
	TExeTm2 tm;
	TSnap::GetTriadParticipAll(G, tridCnt);
	printf("time costs: %s\n", tm.GetStr());
	TStr ofnm = GFNm.GetFPath()+"groundtruth.dat";
	FILE* fw=fopen(ofnm.CStr(), "w");
	fprintf(fw, "# Time cost: %.2f seconds\n", tm.GetSecs());
	fprintf(fw, "# Nodes: %.0f\n", N);
	for (int i=0; i<tridCnt.Len(); i++) {
		int card = tridCnt[i].Val1;
		int freq = tridCnt[i].Val2;
		double prob = freq/N;
		fprintf(fw, "%d\t%d\t%.6e\n", card, freq, prob);
	}
	fclose(fw);
}

void eval_efficiency(ExamMgr& ExM){
	TIntPrV tridCnt;
	TExeTm2 tm;
//	TSnap::GetTriadParticipAll(ExM.GFull, tridCnt);
//	double secs = tm.GetSecs();
//	printf("Full graph: %.2f secs.\n", secs);

	PNEGraph G = PNEGraph::TObj::New();

	double ps[] = {.1, .15, .2, .25, .3};
	for (int i=0; i<5; i++){
		ExM.GetSampledGraph(G, ps[i]);
		tm.Tick();
		TSnap::GetTriadParticipAll(G, tridCnt);
		printf("PEdge = %g: %.2f secs.\n", ps[i], tm.GetSecs());
	}
}

void count_trids_after_sampling(ExamMgr& ExM){
	PNEGraph G = PNEGraph::TObj::New();
	for (int i=1; i<=10; i++){
		double pe=0.01*i;
		ExM.GetSampledGraph(G, pe);
		int nt = count_trids(G);
		printf("%g:  %g\n", pe, nt/pow(pe,3));
	}
}

void count_trids_per_node(const TStr& GFNm){
	PNEGraph G = TSnap::LoadEdgeList<PNEGraph>(GFNm);
	int nid, ntrids;
	TIntPrV NIdTrids;
	for(PNEGraph::TObj::TNodeI NI=G->BegNI(); NI<G->EndNI(); NI++){
		nid = NI.GetId();
		ntrids = TSnap::GetNodeTriadsAll(G, nid);
		NIdTrids.Add(TIntPr(nid, ntrids));
	}
	BIO::SaveIntPrV(NIdTrids, GFNm.GetFPath()+"NodeTrids.dat", "NId, NTrids");
}

int main(int argc, char* argv[]){
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	const TStr GFNm = Env.GetIfArgPrefixStr("-i:", "test.graph", "Input graph");
	const int W = Env.GetIfArgPrefixInt("-w:", 1000, "W. Default 1000");
	const int CPU = Env.GetIfArgPrefixInt("-c:", 8, "Cores. Default 8");
	const int Rpt = Env.GetIfArgPrefixInt("-r:", 12, "Repeat. Default 12");
	const double Pe = Env.GetIfArgPrefixFlt("-p:", 0.1, "Edge sampling rate. Default 0.1");
	const TStr Fmts = Env.GetIfArgPrefixStr("-c:", "", "What to compute:"
				"\n\tc: count trids per node"
				"\n\tg: get groundtruth"
				"\n\te: compare efficiency");
	if (Env.IsEndOfRun()) return 0;
	TExeTm2 tm;
	if (Fmts.SearchCh('c') != -1) count_trids_per_node(GFNm);
	if (Fmts.SearchCh('g') != -1) gen_groundtruth(GFNm);
	if (Fmts.SearchCh('e') != -1){
		ExamMgr ExM(GFNm, W, Pe, CPU, Rpt);
		eval_efficiency(ExM);
	}
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

