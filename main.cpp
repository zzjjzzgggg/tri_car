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
const static int PerRpt=10;

void em_sub(TFltV& ThV, int& NSuc, ExamMgr& ExM){
	const double PDelta = pow(PEdge, 3);
	NSuc = 0;
	PNEGraph G = PNEGraph::TObj::New();
	TIntPrV TridCnt;
	TIntH Freq;
	for (int rpt=0; rpt<PerRpt; rpt++){
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
			for (int i=0; i<=W; i++) ThV[i] += EM.ThV[i];
			NSuc++;
		}
	}
	printf("Experiment repeats %d times, and %d succeeded.\n", PerRpt, NSuc);
}

void em_multi(){
	ExamMgr ExM(GFNm, PEdge);
	TFltV ThV1(W+1), ThV2(W+1), ThV3(W+1), ThV4(W+1), ThV5(W+1);
	int NSuc1=0, NSuc2=0, NSuc3=0, NSuc4=0, NSuc5=0;
	std::vector<std::function<void()>> vec {
		[&ThV1, &NSuc1, &ExM] () { em_sub(ThV1, NSuc1, ExM); },
		[&ThV2, &NSuc2, &ExM] () { em_sub(ThV2, NSuc2, ExM); },
		[&ThV3, &NSuc3, &ExM] () { em_sub(ThV3, NSuc3, ExM); },
		[&ThV4, &NSuc4, &ExM] () { em_sub(ThV4, NSuc4, ExM); },
		[&ThV5, &NSuc5, &ExM] () { em_sub(ThV5, NSuc5, ExM); },
	};
	std::vector<std::thread> threads;
	for(const auto& f: vec) threads.emplace_back((std::function<void()>)f);
	for(std::thread& t: threads) t.join();

	printf("Saving...\n");
	for (int i=0; i<=W; i++) {
		ThV1[i] += (ThV2[i] + ThV3[i] + ThV4[i] + ThV5[i]);
		ThV1[i] /= (NSuc1 + NSuc2 + NSuc3 + NSuc4 + NSuc5);
	}
	TStr OFnm = TStr::Fmt("est_%s_p%g_r%d.dist", GFNm.GetFMid().CStr(),
		PEdge, PerRpt*vec.size());
	BIO::SaveFltsWithIdx(ThV1, OFnm);
	printf("Saved to %s\n", OFnm.CStr());
}

int main(int argc, char* argv[]){
	em_multi();
	return 0;
}

