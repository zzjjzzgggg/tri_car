/*
 * main.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */
#include "stdafx.h"
#include "ExamMgr.h"

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

void verify(ExamMgr& ExM){
	PNEGraph G = PNEGraph::TObj::New();
	ExM.GetSampledGraph(G);
	TIntPrV TridCnt;
	TSnap::GetTriadParticipAll(G, TridCnt);
	int g = 0;
	for (int i=0; i<TridCnt.Len(); i++) {
		if (TridCnt[i].Val1 > 0) g += TridCnt[i].Val2;
	}
	printf("g: %d\n", g);

	TIntPrV CarFreq;
	BIO::LoadIntPrVec(ExM.GFNm.GetFPath()+"groundtruth.dat", CarFreq);
	TIntFltPrV Th;
	double norm = ExM.N - CarFreq[0].Val2;
	for (int i=1; i<CarFreq.Len(); i++){
		const int card = CarFreq[i].Val1, freq = CarFreq[i].Val2;
		Th.Add(TIntFltPr(card, freq/norm));
	}

	TIntFltPrV Th_hat;
	BIO::LoadIntFltPrVec(ExM.GetNTHFNm(), Th_hat);

	double qth=0;
	for (int i=0; i<Th.Len(); i++){
		qth += Th[i].Val2*TSpecFunc::Binomial(0, Th[i].Val1, pow(ExM.PEdge,3));
	}
	printf("q_th: %.6f\n", qth);

	double qth_hat=0;
	for (int i=1; i<Th_hat.Len(); i++){
		qth_hat += Th_hat[i].Val2*TSpecFunc::Binomial(0, Th_hat[i].Val1, pow(ExM.PEdge,3));
	}
	printf("q_th_hat: %.6f\n", qth_hat);

	double nest = g / (1-qth);
	printf("N est: %.2f (%.0f)\n", nest/norm, norm);

	double nest_hat = g / (1-qth_hat);
	printf("N est: %.2f (%.0f)\n", nest_hat/norm, norm);

}

void add_edge(ExamMgr& ExM){
	TIntPrV NodeTriads;
	BIO::LoadIntPrVec(ExM.GFNm.GetFPath()+"nodentrids.dat", NodeTriads);
	TIntV NoNds;
	for (int i=0; i<NodeTriads.Len(); i++){
		if (NodeTriads[i].Val2==0) NoNds.Add(NodeTriads[i].Val1);
	}
	NoNds.Shuffle(TInt::Rnd);
	for (int i=0; i<2000; i++){
		int id1 = NoNds[i];
		PNEGraph::TObj::TNodeI ni = ExM.GFull->GetNI(id1);
		int com = ni.GetNbhNId(0);
		ni = ExM.GFull->GetNI(com);
		for (int d=0; d<ni.GetDeg(); d++){
			int id2=ni.GetNbhNId(d);
			if (ExM.GFull->GetNI(id2).GetDeg()==1) {
				ExM.GFull->AddEdge(id1, id2);
				break;
			}
		}
	}
	TSnap::SaveEdgeList<PNEGraph>(ExM.GFull, ExM.GFNm.GetFPath()+"newgraph.gz");
}

int main(int argc, char* argv[]){
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	const TStr GFNm = Env.GetIfArgPrefixStr("-i:", "test.graph", "Input graph");
	const int W = Env.GetIfArgPrefixInt("-w:", 10000, "W");
	const int CPU = Env.GetIfArgPrefixInt("-n:", 8, "Cores to use, max=8");
	const int Rpt = Env.GetIfArgPrefixInt("-r:", 12, "Repeat");
	const double Pe = Env.GetIfArgPrefixFlt("-p:", 0.1, "Edge sampling rate");
	if (Env.IsEndOfRun()) return 0;

	TExeTm2 tm;
	ExamMgr ExM(GFNm, CPU, W, Pe, Rpt);
	add_edge(ExM);
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

