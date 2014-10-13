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

void em_sub_n_unknown(const double PEdge, TFltV& ThV, int& NSuc, ExamMgr& ExM){
	const double PDelta = pow(PEdge, 3);
	NSuc = 0;
	PNEGraph G = PNEGraph::TObj::New();
	TIntPrV TridCnt;
	for (int rpt=0; rpt<PerRpt; rpt++){
		printf("************ rpt = %d ************\n", rpt);
		printf("Stating sampled data...\n");
		ExM.GetSampledGraph(G, PEdge);
		TridCnt.Clr();
		TSnap::GetTriadParticipAll(G, TridCnt);
		TCEMGeneral EM(W, PDelta, TridCnt);
		printf("Em running...\n");
		if (EM.Run()) {
			for (int i=0; i<=W; i++) ThV[i] += EM.ThV[i];
			NSuc++;
		}
	}
	printf("Experiment repeats %d times, and %d succeeded.\n", PerRpt, NSuc);
}

void em_multi_n_unknown(const double PEdge){
	ExamMgr ExM(GFNm);
	TFltV ThV1(W+1), ThV2(W+1), ThV3(W+1), ThV4(W+1), ThV5(W+1);
	int NSuc1=0, NSuc2=0, NSuc3=0, NSuc4=0, NSuc5=0;
	std::vector<std::function<void()>> vec = {
		[&PEdge, &ThV1, &NSuc1, &ExM] { em_sub_n_unknown(PEdge, ThV1, NSuc1, ExM); },
		[&PEdge, &ThV2, &NSuc2, &ExM] { em_sub_n_unknown(PEdge, ThV2, NSuc2, ExM); },
		[&PEdge, &ThV3, &NSuc3, &ExM] { em_sub_n_unknown(PEdge, ThV3, NSuc3, ExM); },
		[&PEdge, &ThV4, &NSuc4, &ExM] { em_sub_n_unknown(PEdge, ThV4, NSuc4, ExM); },
		[&PEdge, &ThV5, &NSuc5, &ExM] { em_sub_n_unknown(PEdge, ThV5, NSuc5, ExM); },
	};
	std::vector<std::thread> threads;
	for(const auto& f: vec) threads.emplace_back((std::function<void()>)f);
	for(std::thread& t: threads) t.join();

	printf("Saving...\n");
	for (int i=0; i<=W; i++) {
		ThV1[i] += (ThV2[i] + ThV3[i] + ThV4[i] + ThV5[i]);
		ThV1[i] /= (NSuc1 + NSuc2 + NSuc3 + NSuc4 + NSuc5);
	}
	TStr OFnm = TStr::Fmt("n_th_%s_p%g_r%d.dist", GFNm.GetFMid().CStr(),
		PEdge, PerRpt*vec.size());
	BIO::SaveFltsWithIdx(ThV1, OFnm, "# The first line denotes the estimated N");
	printf("Saved to %s\n", OFnm.CStr());
}

int main(int argc, char* argv[]){
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs();
	const double p = Env.GetIfArgPrefixFlt("-p:", 0.1, "Edge sampling rate. Default 0.1.");
	const TStr GFNm = Env.GetIfArgPrefixStr("-i:", "graph.txt", "Input graph");
	TExeTm2 tm;
//	verify();
//	em_single();
//	stat_trids();
//	gen_ground_truth();
	em_multi_n_known(p);
//	em_multi_n_unknown();
//	eval_efficiency();
//	count_triangles();
//	dist_triangles(50);
//	dist_degree(20);
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

