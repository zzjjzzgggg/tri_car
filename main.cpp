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

using namespace std;

const TStr DG_HEPTH = "../cit-HepTh.gz";
const TStr DG_TEST = "../test_dir_graph.txt";
const TStr DG_TRIADS = "../non_zero_triads.digraph.gz";

void test(){
	printf("%s\n", DG_TRIADS.GetFBase().CStr());
	printf("%s\n", DG_TRIADS.GetFMid().CStr());
	printf("%s\n", DG_TRIADS.GetFPath().CStr());
	printf("%s\n", DG_TRIADS.GetFNmStr("xx.dd").CStr());
}

void gen_test_data(){
	TIntV zeroNodes;
	PNEGraph G = TSnap::LoadEdgeList<PNEGraph>(DG_HEPTH);
	for (PNEGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
		const int Triads = TSnap::GetNodeTriads(G, NI.GetId());
		if(Triads == 0) zeroNodes.Add(NI.GetId());
	}
	for(int i=0; i<zeroNodes.Len(); i++) G->DelNode(zeroNodes[i]);
	TSnap::SaveEdgeList<PNEGraph>(G, "../non_zero_triads.digraph");
}

void stat_tc_dist(PNEGraph G, TStr ofnm){
	TIntPrV tridCnt;
	TExeTm2 tm;
	TSnap::GetTriadParticipAll(G, tridCnt);
	printf("time costs: %s\n", tm.GetStr());
	FILE* fw=fopen(ofnm.CStr(), "w");
	fprintf(fw, "# time cost: %.2f seconds\n", tm.GetSecs());
	double N = G->GetNodes();
	for (int i=0; i<tridCnt.Len(); i++) {
		int card = tridCnt[i].Val1;
		int freq = tridCnt[i].Val2;
		double prob = freq/N;
		fprintf(fw, "%d\t%d\t%.6e\n", card, freq, prob);
	}
	fclose(fw);
}

void ground_truth(TStr gfnm, TStr ofnm){
	PNEGraph G = TSnap::LoadEdgeList<PNEGraph>(gfnm);
	stat_tc_dist(G, ofnm);
}

void sample_stream(PNEGraph& G, TStr gfnm, const double p){
	TRnd rnd;
	TSsParser Ss(gfnm);
	while(Ss.Next()){
		if(rnd.GetUniDev() > p) continue;
		int src = Ss.GetInt(0), dst = Ss.GetInt(1);
		if(!G->IsNode(src)) G->AddNode(src);
		if(!G->IsNode(dst)) G->AddNode(dst);
		G->AddEdge(src, dst);
	}
	printf("sampled stream: \n\t nodes: %d, edges: %d\n", G->GetNodes(), G->GetEdges());
}

void state_sampled_stream(){
	PNEGraph G = PNEGraph::TObj::New();
	sample_stream(G, DG_HEPTH, 0.2);
	stat_tc_dist(G, "../sampled_p0.2.dist");
}

void em_estimate(){
	double p_delta = pow(0.2, 3);
	int W = 1000, M=200, N=27770;
	TIntH freq;
	BIO::LoadIntH("../sampled_p0.2.dist", freq);
//	for(int i=0; i<freq.Len(); i++) {
//		if(freq.GetKey(i)>M) M = freq.GetKey(i);
//	}
	TCEM EM(W, M, N, p_delta, freq);
	bool suc = EM.Run();
	if (suc) {
		EM.Save("../estimates.dist");
		printf("Converged and results are saved.");
	}else
		printf("Not converging, quit without saving results.");
}

void em_sub(TFltV& AvgThV, int& NSuc, const TStr& GFnm, const int W, const int M, const int N, const double p, const int Rpt){
	const double P = pow(p, 3);
	NSuc = 0;
	PNEGraph G = PNEGraph::TObj::New();
	TIntPrV TridCnt;
	TIntH Freq;
	for (int rpt=0; rpt<Rpt; rpt++){
		printf("************ rpt = %d ************\n", rpt);
		G->Clr(); TridCnt.Clr(); Freq.Clr();
		sample_stream(G, GFnm, p);
		printf("Stating sampled data...\n");
		TSnap::GetTriadParticipAll(G, TridCnt);
		for (int i=0; i<TridCnt.Len(); i++){
			const int key = TridCnt[i].Val1;
			const int val = TridCnt[i].Val2;
			if (key > M) break;
			Freq(key) = val;
		}
		TCEM EM(W, M, N, P, Freq);
		printf("Em running...\n");
		if (EM.Run()) {
			for (int i=0; i<=W; i++) AvgThV[i] += EM.ThV[i];
			NSuc++;
		}
	}
	printf("Experiment repeats %d times, and %d succeeded.\n", Rpt, NSuc);
}

void em_multi(){
	int W=1000, M=200, N=27770, PerRpt=2;
	double p=0.2;
	TStr GFnm = DG_HEPTH;
	TFltV AvgThV1(W+1), AvgThV2(W+1), AvgThV3(W+1), AvgThV4(W+1), AvgThV5(W+1);
	int NSuc1=0, NSuc2=0, NSuc3=0, NSuc4=0, NSuc5=0;
	std::vector<std::function<void()>> vec {
		[&AvgThV1, &NSuc1, &GFnm, W, M, N, p, PerRpt] () { em_sub(AvgThV1, NSuc1, GFnm, W, M, N, p, PerRpt); },
		[&AvgThV2, &NSuc2, &GFnm, W, M, N, p, PerRpt] () { em_sub(AvgThV2, NSuc2, GFnm, W, M, N, p, PerRpt); },
//		[&AvgThV3, &NSuc3, &GFnm, W, M, N, p, Rpt] () { em_sub(AvgThV3, NSuc3, GFnm, W, M, N, p, Rpt/5); },
//		[&AvgThV4, &NSuc4, &GFnm, W, M, N, p, Rpt] () { em_sub(AvgThV4, NSuc4, GFnm, W, M, N, p, Rpt/5); },
//		[&AvgThV5, &NSuc5, &GFnm, W, M, N, p, Rpt] () { em_sub(AvgThV5, NSuc5, GFnm, W, M, N, p, Rpt/5); },
	};
	std::vector<std::thread> threads;
	for(const auto& f: vec) threads.emplace_back((std::function<void()>)f);
	for(std::thread& t: threads) t.join();

	printf("Saving...\n");
	for (int i=0; i<=W; i++) {
		AvgThV1[i] += (AvgThV2[i] + AvgThV3[i] + AvgThV4[i] + AvgThV5[i]);
		AvgThV1[i] /= (NSuc1 + NSuc2 + NSuc3 + NSuc4 + NSuc5);
	}
	BIO::SaveFltsWithIdx(AvgThV1, TStr::Fmt("../est_%s_p%g_r%d.dist", GFnm.GetFMid().CStr(), p, PerRpt*vec.size()));

}

int main(int argc, char* argv[]){
//	test();
//	gen_test_data();
//	ground_truth(DG_TRIADS, "../triads.dist");
//	state_sampled_stream();
//	em_estimate();
//	em_sub(DG_HEPTH, 1000, 200, 27770, 0.2, 100);
	em_multi();
	return 0;
}

