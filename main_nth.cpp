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
		printf("Em running...\n");
		if (EM.Run()) {
			for (int i=0; i<=ExM.W; i++) ThV[i] += EM.ThV[i];
			NSuc++;
		}
	}
	for (int i=0; i<ThV.Len(); i++) ThV[i] /= NSuc;
	printf("Experiment repeats %d times, and %d succeeded.\n", ExM.Rpt, NSuc);
}

void em_multi(ExamMgr& ExM){
	TFltV ThV1(ExM.W+1), ThV2(ExM.W+1), ThV3(ExM.W+1), ThV4(ExM.W+1),
			ThV5(ExM.W+1), ThV6(ExM.W+1), ThV7(ExM.W+1), ThV8(ExM.W+1);
	std::vector<std::function<void()>> Tasks = {
		[&ExM, &ThV1] { em_sub(ExM, ThV1); },
		[&ExM, &ThV2] { em_sub(ExM, ThV2); },
		[&ExM, &ThV3] { em_sub(ExM, ThV3); },
		[&ExM, &ThV4] { em_sub(ExM, ThV4); },
		[&ExM, &ThV5] { em_sub(ExM, ThV5); },
		[&ExM, &ThV6] { em_sub(ExM, ThV6); },
		[&ExM, &ThV7] { em_sub(ExM, ThV7); },
		[&ExM, &ThV8] { em_sub(ExM, ThV8); },
	};
	std::vector<std::thread> threads;
	for (int i=0; i<ExM.CPU; i++) threads.emplace_back(Tasks[i]);
	for(std::thread& t: threads) t.join();

	printf("Saving...\n");
	for (int i=0; i<ThV1.Len(); i++)
		ThV1[i] = (ThV1[i] + ThV2[i] + ThV3[i] + ThV4[i] + ThV5[i] + ThV6[i] + ThV7[i] + ThV8[i]) / ExM.CPU;
	const TStr OFnm = ExM.GetNTHFNm();
	BIO::SaveFltsWithIdx(ThV1, OFnm, TStr::Fmt("First line is \hat{N}. Repeated: %d", ExM.Rpt*ExM.CPU));
	printf("Saved to %s\n", OFnm.CStr());
}

int main(int argc, char* argv[]){
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	const TStr GFNm = Env.GetIfArgPrefixStr("-i:", "test.graph", "Input graph");
	const int W = Env.GetIfArgPrefixInt("-w:", 1000, "W. Default 1000");
	const int CPU = Env.GetIfArgPrefixInt("-c:", 8, "Cores. Default 8");
	const int Rpt = Env.GetIfArgPrefixInt("-r:", 12, "Repeat. Default 12");
	const double Pe = Env.GetIfArgPrefixFlt("-p:", 0.1, "Edge sampling rate. Default 0.1");
	if (Env.IsEndOfRun()) return 0;

	TExeTm2 tm;
	ExamMgr ExM(GFNm, W, Pe, CPU, Rpt);
	em_multi(ExM);
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

