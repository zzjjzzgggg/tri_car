/*
 * main.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */
#include "stdafx.h"
#include "ExamMgr.h"


void stat_trids(ExamMgr& ExM){
	PNEGraph G = ExM.GFull;
	TIntPrV TridCnt;
	for(PNEGraph::TObj::TNodeI NI=G->BegNI(); NI<G->EndNI(); NI++){
		int nid = NI.GetId();
		int ntrids = TSnap::GetNodeTriadsAll(G, nid);
		TridCnt.Add(TIntPr(nid, ntrids));
	}
	BIO::SaveIntPrV(TridCnt, "NodeNTrids.dat");
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
	BIO::LoadIntPrV(ExM.GetGTFNm(), CarFreq);
	TIntFltPrV Th;
	double norm = ExM.N - CarFreq[0].Val2;
	for (int i=1; i<CarFreq.Len(); i++){
		const int card = CarFreq[i].Val1, freq = CarFreq[i].Val2;
		Th.Add(TIntFltPr(card, freq/norm));
	}

	TFltV Th_hat;
	BIO::LoadFltV(ExM.GetNTHFNm(), Th_hat, 1);

	double qth=0;
	for (int i=0; i<Th.Len(); i++){
		qth += Th[i].Val2*TSpecFunc::Binomial(0, Th[i].Val1, pow(ExM.PEdge,3));
	}
	printf("q_th: %.6f\n", qth);

	double qth_hat=0;
	for (int i=1; i<Th_hat.Len(); i++){
		qth_hat += Th_hat[i]*TSpecFunc::Binomial(0, i, pow(ExM.PEdge,3));
	}
	printf("q_th_hat: %.6f\n", qth_hat);

	double nest = g / (1-qth);
	printf("N est: %.6f (%.0f)\n", nest/norm, norm);

	double nest_hat = g / (1-qth_hat);
	printf("[1] N est: %.6f\n", nest_hat/norm);

	double nest_hat2=0;
	for (int i=1; i<Th_hat.Len(); i++){
		nest_hat2 += Th_hat[i]/(1-TSpecFunc::Binomial(0, i, pow(ExM.PEdge,3)));
	}
	nest_hat2 *= g;
	printf("[2] N est: %.6f\n", nest_hat2/norm);

	Th_hat.Clr();
	BIO::LoadFltV(ExM.GFNm.GetFPath()+"nthp_mention-higgs_W1K_p0.1.dist", Th_hat, 1);
	double nest_hat3=0;
	for (int i=1; i<Th_hat.Len(); i++){
		nest_hat3 += Th_hat[i];
	}
	printf("[3] N est: %.6f\n", nest_hat3);
	nest_hat3 *= g;
	printf("[3] N est: %.6f\n", nest_hat3/norm);
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

	TExeTm tm;
	ExamMgr ExM(GFNm, CPU, W, Pe, Rpt);
//	verify(ExM);
	stat_trids(ExM);
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

