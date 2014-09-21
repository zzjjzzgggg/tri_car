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

const static TStr GFNm = DG_HEPTH;
const static int W = 1000;
const static int N = 27770;
const static double PEdge = 0.3;
const static int PerRpt=20;


void em_sub_n_known(TFltV& ThV, int& NSuc, ExamMgr& ExM){
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
		TCEM EM(W, N, PDelta, TridCnt);
		printf("Em running...\n");
		if (EM.Run()) {
			for (int i=0; i<=W; i++) ThV[i] += EM.ThV[i];
			NSuc++;
		}
	}
	printf("Experiment repeats %d times, and %d succeeded.\n", PerRpt, NSuc);
}

void em_multi_n_known(){
	ExamMgr ExM(GFNm);
	TFltV ThV1(W+1), ThV2(W+1), ThV3(W+1), ThV4(W+1), ThV5(W+1);
	int NSuc1=0, NSuc2=0, NSuc3=0, NSuc4=0, NSuc5=0;
	std::vector<std::function<void()>> vec = {
		[&ThV1, &NSuc1, &ExM] { em_sub_n_known(ThV1, NSuc1, ExM); },
		[&ThV2, &NSuc2, &ExM] { em_sub_n_known(ThV2, NSuc2, ExM); },
		[&ThV3, &NSuc3, &ExM] { em_sub_n_known(ThV3, NSuc3, ExM); },
		[&ThV4, &NSuc4, &ExM] { em_sub_n_known(ThV4, NSuc4, ExM); },
		[&ThV5, &NSuc5, &ExM] { em_sub_n_known(ThV5, NSuc5, ExM); },
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

void em_sub_n_unknown(TFltV& ThV, int& NSuc, ExamMgr& ExM){
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

void em_single(){
	ExamMgr ExM(GFNm);
	TFltV ThV(W+1);
	int NSuc=0;
	em_sub_n_unknown(ThV, NSuc, ExM);
}

void em_multi_n_unknown(){
	ExamMgr ExM(GFNm);
	TFltV ThV1(W+1), ThV2(W+1), ThV3(W+1), ThV4(W+1), ThV5(W+1);
	int NSuc1=0, NSuc2=0, NSuc3=0, NSuc4=0, NSuc5=0;
	std::vector<std::function<void()>> vec = {
		[&ThV1, &NSuc1, &ExM] { em_sub_n_unknown(ThV1, NSuc1, ExM); },
		[&ThV2, &NSuc2, &ExM] { em_sub_n_unknown(ThV2, NSuc2, ExM); },
		[&ThV3, &NSuc3, &ExM] { em_sub_n_unknown(ThV3, NSuc3, ExM); },
		[&ThV4, &NSuc4, &ExM] { em_sub_n_unknown(ThV4, NSuc4, ExM); },
		[&ThV5, &NSuc5, &ExM] { em_sub_n_unknown(ThV5, NSuc5, ExM); },
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

void eval_efficiency(){
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

void verify(){
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
	TExeTm2 tm;
//	verify();
//	em_single();
//	em_multi();
	em_multi_n_unknown();
//	eval_efficiency();
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

