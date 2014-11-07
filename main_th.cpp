/*
 * main.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#include "stdafx.h"
#include "TCEM.h"
#include "ExamMgr.h"


void em_sub(const int id, ExamMgr& ExM, TFltV& ThV){
	int NSuc = 0; TIntPrV gV;
	for (int rpt=0; rpt<ExM.Rpt; rpt++){
		printf("rpt = %d\n", rpt);
		ExM.Sample(gV);
		TCEM EM(ExM.W, ExM.N, ExM.GetPdUU(), gV);
		printf("[%d] M: %d\n", id, EM.M);
		if (EM.Run()) {
			for (int i=0; i<=ExM.W; i++) ThV[i] += EM.ThV[i];
			NSuc++;
		}
//		printf("[%d] EM: %d\n", id, EM.Iters);
	}
	for (int i=0; i<ThV.Len(); i++) ThV[i] /= NSuc;
	printf("[%d] Experiment repeats %d times, and %d succeeded.\n", id, ExM.Rpt, NSuc);
}

void trim_tail(ExamMgr& ExM, TFltV& ThV){
	double minval=1, Pd = pow(ExM.PEdge,3);
	int Wp=1;
	for (int i=0; i<=ExM.W; i++) {
		if (ThV[i] < minval){
			minval = ThV[i];
			Wp = i;
		}
	}
	double qth = 0, rem = 0; // q_theta
	for (int i=Wp+1; i<=ExM.W; i++) {
		rem += ThV[i];
		ThV[i] = 0;
	}
	for (int i=0; i<=Wp; i++){
		ThV[i] /= (1-rem);
		qth += ThV[i]*TSpecFunc::Binomial(0, i, Pd);
	}
	printf("min val = %.2e   Wp=%d  rem = %.2e\n", minval, Wp, rem);
}

void em_multi(ExamMgr& ExM){
	TExeTm tm;
	TFltV ThVs[ExM.CPU];
	for (int i=0; i<ExM.CPU; i++) ThVs[i] = TFltV(ExM.W+1);
	std::vector<std::thread> threads;
	for (int i=0; i<ExM.CPU; i++) threads.emplace_back([i, &ExM, &ThVs] { em_sub(i, ExM, ThVs[i]); });
	for(std::thread& t: threads) t.join();
	for (int i=0; i<=ExM.W; i++){
		for (int n=1; n<ExM.CPU; n++) ThVs[0][i] += ThVs[n][i];
		ThVs[0][i] /= ExM.CPU;
	}
	if (ExM.TrimTail) trim_tail(ExM, ThVs[0]);
	const TStr OFnm = ExM.GetTHFNm();
	BIO::SaveFltVWithIdx(ThVs[0], OFnm, TStr::Fmt("# Nodes: %d\n# Repeated: %d. \n# Avg time cost: %.2f secs.", ExM.N, ExM.GetRpt(), tm.GetSecs()/ExM.GetRpt()));
	printf("Saved to %s\n", OFnm.CStr());
}


int main(int argc, char* argv[]){
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	const TStr GFNm = Env.GetIfArgPrefixStr("-i:", "", "Input graph");
	const int W = Env.GetIfArgPrefixInt("-w:", 10000, "W");
	const int Rpt = Env.GetIfArgPrefixInt("-r:", 10, "Repeat times");
	const double Pe = Env.GetIfArgPrefixFlt("-p:", 0.1, "Edge sampling rate");
	const bool TrimTail = Env.GetIfArgPrefixBool("-t:", false, "Trim tail");
	if (Env.IsEndOfRun()) return 0;

	TExeTm2 tm;
	ExamMgr ExM(GFNm, W, Pe, Rpt, TrimTail);
	em_multi(ExM);
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

