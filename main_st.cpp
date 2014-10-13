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

void gen_ground_truth(const TStr& GFNm, const int W){
	PNEGraph G = TSnap::LoadEdgeList<PNEGraph>(GFNm);
	double N = G->GetNodes();
	TIntPrV tridCnt;
	TExeTm2 tm;
	TSnap::GetTriadParticipAll(G, tridCnt);
	printf("time costs: %s\n", tm.GetStr());
	TStr ofnm = TStr::Fmt("groundtruth_W%dK.dat", W/1000);
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

void eval_efficiency(const TStr& GFNm){
	ExamMgr ExM(GFNm);
	TIntPrV tridCnt;
	TExeTm2 tm;
	TSnap::GetTriadParticipAll(ExM.GFull, tridCnt);
	double secs = tm.GetSecs();
	printf("Full graph: %.2f secs.\n", secs);

	PNEGraph G = PNEGraph::TObj::New();

	double ps[] = {.1, .2, .3};
	for (int i=0; i<3; i++){
		ExM.GetSampledGraph(G, ps[i]);
		tm.Tick();
		TSnap::GetTriadParticipAll(G, tridCnt);
		printf("PEdge = %g: %.2f secs.\n", ps[i], tm.GetSecs());
	}
}

void verify(const TStr& GFNm, const double PEdge){
	ExamMgr ExM(GFNm);
	PNEGraph G = PNEGraph::TObj::New();
	TIntPrV TridCnt;
	double g=0;
	for (int i=0; i<20; i++){
		ExM.GetSampledGraph(G, PEdge);
		TridCnt.Clr();
		TSnap::GetTriadParticipAll(G, TridCnt);
		for (int j=0; j<TridCnt.Len(); j++)
			if (TridCnt[j].Val1!=0) g+=TridCnt[j].Val2;
	}
	printf("g = %.2f\n", g/20);
}

int main(int argc, char* argv[]){
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs();
	const double p = Env.GetIfArgPrefixFlt("-p:", 0.1, "Edge sampling rate. Default 0.1.");
	const TStr GFNm = Env.GetIfArgPrefixStr("-i:", "graph.txt", "Input graph");
	TExeTm2 tm;
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}
