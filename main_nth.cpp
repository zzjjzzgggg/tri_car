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
#include "TCEMGeneral.h"
#include "ExamMgr.h"

using namespace std;

void em_sub(ExamMgr& ExM, TFltV& ThV){
	int NSuc = 0;
	PNEGraph G = PNEGraph::TObj::New();
	TIntPrV TridCnt;
	for (int rpt=0; rpt<ExM.Rpt; rpt++){
		printf("rpt = %d\n", rpt);
		ExM.GetSampledGraph(G);
		TridCnt.Clr();
		TSnap::GetTriadParticipAll(G, TridCnt);
		TCEMGeneral EM(ExM.W, pow(ExM.PEdge, 3), TridCnt);
		printf("Sampled: nodes: %d, edges: %d, M: %d\n", G->GetNodes(), G->GetEdges(), EM.M);
		if (EM.Run()) {
			for (int i=0; i<=ExM.W; i++) ThV[i] += EM.ThV[i];
			NSuc++;
		}
	}
	for (int i=0; i<ThV.Len(); i++) ThV[i] /= NSuc;
	printf("Experiment repeats %d times, and %d succeeded.\n", ExM.Rpt, NSuc);
}

void em_multi(ExamMgr& ExM){
	TFltV ThVs[ExM.CPU];
	for (int i=0; i<ExM.CPU; i++) ThVs[i] = TFltV(ExM.W+1);

	std::vector<std::thread> threads;
	for (int i=0; i<ExM.CPU; i++) threads.emplace_back([i, &ExM, &ThVs] { em_sub(ExM, ThVs[i]); });
	for(std::thread& t: threads) t.join();

	for (int i=0; i<=ExM.W; i++){
		for (int n=1; n<ExM.CPU; n++) ThVs[0][i] += ThVs[n][i];
		ThVs[0][i] /= ExM.CPU;
	}

	printf("Saving...\n");
	const TStr OFnm = ExM.GetNTHFNm();
	BIO::SaveFltsWithIdx(ThVs[0], OFnm, TStr::Fmt("First line is N est. Repeated: %d", ExM.Rpt*ExM.CPU));
	printf("Saved to %s\n", OFnm.CStr());
}

int main(int argc, char* argv[]){
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	const TStr GFNm = Env.GetIfArgPrefixStr("-i:", "test.graph", "Input graph");
	const int W = Env.GetIfArgPrefixInt("-w:", 10000, "W");
	const int CPU = Env.GetIfArgPrefixInt("-n:", 8, "Cores to use, max=8");
	const int Rpt = Env.GetIfArgPrefixInt("-r:", 12, "Repeat times");
	const double Pe = Env.GetIfArgPrefixFlt("-p:", 0.1, "Edge sampling rate");
	if (Env.IsEndOfRun()) return 0;

	TExeTm2 tm;
	ExamMgr ExM(GFNm, W, Pe, CPU, Rpt);
	em_multi(ExM);
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

