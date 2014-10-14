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


void em_sub(TFltV& ThV, int& NSuc, ExamMgr& ExM){
	NSuc = 0;
	PNEGraph G = PNEGraph::TObj::New();
	TIntPrV TridCnt;
	for (int rpt=0; rpt<PerRpt; rpt++){
		printf("rpt = %d\n", rpt);
		ExM.GetSampledGraph(G);
		TridCnt.Clr();
		TSnap::GetTriadParticipAll(G, TridCnt);
		TCEM EM(ExM.W, ExM.N, pow(ExM.PEdge, 3), TridCnt);
		printf("Sampled: nodes: %d, edges: %d, M: %d\n", G->GetNodes(), G->GetEdges(), EM.M);
		if (EM.Run()) {
			for (int i=0; i<=ExM.W; i++) ThV[i] += EM.ThV[i];
			NSuc++;
		}
	}
	printf("Experiment repeats %d times, and %d succeeded.\n", PerRpt, NSuc);
}


void em_multi(const TStr& GFNm, const int W, const double Pe){
	ExamMgr ExM(GFNm, W, Pe);
	TFltV ThV1(W+1), ThV2(W+1), ThV3(W+1), ThV4(W+1), ThV5(W+1), ThV6(W+1), ThV7(W+1), ThV8(W+1);
	int NSuc1=0, NSuc2=0, NSuc3=0, NSuc4=0, NSuc5=0, NSuc6=0, NSuc7=0, NSuc8=0;
	std::vector<std::function<void()>> vec = {
		[&ThV1, &NSuc1, &ExM] { em_sub(ThV1, NSuc1, ExM); },
		[&ThV2, &NSuc2, &ExM] { em_sub(ThV2, NSuc2, ExM); },
		[&ThV3, &NSuc3, &ExM] { em_sub(ThV3, NSuc3, ExM); },
		[&ThV4, &NSuc4, &ExM] { em_sub(ThV4, NSuc4, ExM); },
		[&ThV5, &NSuc5, &ExM] { em_sub(ThV5, NSuc5, ExM); },
		[&ThV6, &NSuc6, &ExM] { em_sub(ThV6, NSuc6, ExM); },
		[&ThV7, &NSuc7, &ExM] { em_sub(ThV7, NSuc7, ExM); },
		[&ThV8, &NSuc8, &ExM] { em_sub(ThV8, NSuc8, ExM); },
	};
	std::vector<std::thread> threads;
	for(const auto& f: vec) threads.emplace_back((std::function<void()>)f);
	for(std::thread& t: threads) t.join();

	printf("Saving...\n");
	int NSuc = NSuc1 + NSuc2 + NSuc3 + NSuc4 + NSuc5 + NSuc6 + NSuc7 + NSuc8;
	for (int i=0; i<=W; i++) {
		ThV1[i] += (ThV2[i] + ThV3[i] + ThV4[i] + ThV5[i] + ThV6[i] + ThV7[i] + ThV8[i]);
		ThV1[i] /= NSuc;
	}
	TStr OFnm = TStr::Fmt("%sth_%s_W%dK_p%g.dist", GFNm.GetFPath().CStr(), GFNm.GetFMid().CStr(), W/1000, Pe);
	BIO::SaveFltsWithIdx(ThV1, OFnm, TStr::Fmt("Repeated: %d", PerRpt*vec.size()));
	printf("Saved to %s\n", OFnm.CStr());
}


int main(int argc, char* argv[]){
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs();
	const TStr GFNm = Env.GetIfArgPrefixStr("-i:", "../hepth/cit-HepTh.gz", "Input graph");
	const int W = Env.GetIfArgPrefixInt("-w:", 1000, "W. Default 1000");
	const double Pe = Env.GetIfArgPrefixFlt("-p:", 0.1, "Edge sampling rate. Default 0.1");
	TExeTm2 tm;
	em_multi(GFNm, W, Pe);
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

