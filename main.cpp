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
#include "ExamMgr.h"

using namespace std;

const static TStr GFNm = DG_HEPTH;
const static int W = 1000;
const static int M = 200;
const static int N = 27770;
const static double PEdge = 0.2;

void em_sub(TFltV& AvgThV, int& NSuc, ExamMgr& ExM, const int Rpt){
	const double PDelta = pow(PEdge, 3);
	NSuc = 0;
	PNEGraph G = PNEGraph::TObj::New();
	TIntPrV TridCnt;
	TIntH Freq;
	for (int rpt=0; rpt<Rpt; rpt++){
		printf("************ rpt = %d ************\n", rpt);
		G->Clr(); TridCnt.Clr(); Freq.Clr();
		printf("Stating sampled data...\n");
		ExM.GetSampledGraph(G);
		TSnap::GetTriadParticipAll(G, TridCnt);
		for (int i=0; i<TridCnt.Len(); i++){
			const int key = TridCnt[i].Val1;
			const int val = TridCnt[i].Val2;
			if (key > M) break;
			Freq(key) = val;
		}
		TCEM EM(W, M, N, PDelta, Freq);
		printf("Em running...\n");
		if (EM.Run()) {
			for (int i=0; i<=W; i++) AvgThV[i] += EM.ThV[i];
			NSuc++;
		}
	}
	printf("Experiment repeats %d times, and %d succeeded.\n", Rpt, NSuc);
}

void em_multi(){
	ExamMgr ExM(GFNm, PEdge);
	int PerRpt=1;
	TFltV AvgThV1(W+1), AvgThV2(W+1), AvgThV3(W+1), AvgThV4(W+1), AvgThV5(W+1);
	int NSuc1=0, NSuc2=0, NSuc3=0, NSuc4=0, NSuc5=0;
	std::vector<std::function<void()>> vec {
		[&AvgThV1, &NSuc1, &ExM, PerRpt] () { em_sub(AvgThV1, NSuc1, ExM, PerRpt); },
		[&AvgThV2, &NSuc2, &ExM, PerRpt] () { em_sub(AvgThV2, NSuc2, ExM, PerRpt); },
		[&AvgThV3, &NSuc3, &ExM, PerRpt] () { em_sub(AvgThV3, NSuc3, ExM, PerRpt); },
		[&AvgThV4, &NSuc4, &ExM, PerRpt] () { em_sub(AvgThV4, NSuc4, ExM, PerRpt); },
		[&AvgThV5, &NSuc5, &ExM, PerRpt] () { em_sub(AvgThV5, NSuc5, ExM, PerRpt); },
	};
	std::vector<std::thread> threads;
	for(const auto& f: vec) threads.emplace_back((std::function<void()>)f);
	for(std::thread& t: threads) t.join();

	printf("Saving...\n");
	for (int i=0; i<=W; i++) {
		AvgThV1[i] += (AvgThV2[i] + AvgThV3[i] + AvgThV4[i] + AvgThV5[i]);
		AvgThV1[i] /= (NSuc1 + NSuc2 + NSuc3 + NSuc4 + NSuc5);
	}
	TStr OFnm = TStr::Fmt("est_%s_p%g_r%d.dist", GFNm.GetFMid().CStr(), PEdge, PerRpt*vec.size());
	BIO::SaveFltsWithIdx(AvgThV1, OFnm);
	printf("Saved to %s\n", OFnm.CStr());
}

int main(int argc, char* argv[]){
	em_multi();
	return 0;
}

