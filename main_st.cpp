/*
 * main.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#include <queue>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <thread>

using namespace std;

#include "stdafx.h"
#include "Queue.h"
#include "ExamMgr.h"


/**
 * count the number of triangles in int-pair NIdCnt.
 */
void sub_gt(const int id, Queue<int>& Que, ExamMgr& ExM, TIntPrV& NIdCnt) {
	TExeTm2 tm;
	int nid, ntrids;
	printf("Thread [%d] is ready.\n", id);
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
	printf("Saving...\n");
	FILE* fw=fopen(ExM.GetNTFNm().CStr(), "w");
	fprintf(fw, "# Nodes: %d\n", ExM.N);
	TIntH TriadCntH;
	for (int n=0; n<ExM.CPU; n++) {
		for (int i=0; i<Vs[n].Len(); i++) {
			TriadCntH(Vs[n][i].Val2) ++;
			fprintf(fw, "%d\t%d\n", Vs[n][i].Val1.Val, Vs[n][i].Val2.Val);
		}
	}
	fclose(fw);

	TIntPrV TriadCntV;
	TriadCntH.GetKeyDatPrV(TriadCntV);
	TriadCntV.Sort();
	fw=fopen(ExM.GetGTFNm().CStr(), "w");
	fprintf(fw, "# Time cost: %.2f seconds\n", tm.GetSecs());
	fprintf(fw, "# Nodes: %d\n", ExM.N);
	for (int i=0; i<TriadCntV.Len(); i++) {
		int card = TriadCntV[i].Val1;
		int freq = TriadCntV[i].Val2;
		double prob = freq/(double)ExM.N;
		fprintf(fw, "%d\t%d\t%.6e\n", card, freq, prob);
	}
	fclose(fw);
}

void eval_efficiency(ExamMgr& ExM){
	TIntPrV tridCnt;
	TExeTm2 tm;
	PNEGraph G = PNEGraph::TObj::New();
	double ps[] = {.1, .15, .2, .25, .3};
	for (int i=0; i<5; i++){
		ExM.GetSampledGraph(G, ps[i]);
		tm.Tick();
		TSnap::GetTriadParticipAll(G, tridCnt);
		printf("PEdge = %g: %.2f secs.\n", ps[i], tm.GetSecs());
	}
}


int main(int argc, char* argv[]){
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	const TStr GFNm = Env.GetIfArgPrefixStr("-i:", "test.graph", "Input graph");
	const int CPU = Env.GetIfArgPrefixInt("-n:", 8, "Cores to use, max=8");
	const TStr Fmts = Env.GetIfArgPrefixStr("-c:", "", "What to compute:"
				"\n\tg: get groundtruth"
				"\n\te: compare efficiency");
	if (Env.IsEndOfRun()) return 0;
	TExeTm2 tm;
	ExamMgr ExM(GFNm, CPU);
	if (Fmts.SearchCh('g') != -1) multi_groundtruth(ExM);
	if (Fmts.SearchCh('e') != -1) eval_efficiency(ExM);
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

