/*
 * main.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */

#include "stdafx.h"
#include "TCEMBetaBinomGeneral.h"
#include "ExamMgr.h"

void em_sub(const int id, ExamMgr& ExM, TFlt& alpha, TFltV& ThV){
	int NSuc = 0; TIntPrV gV;
	for (int rpt=0; rpt<ExM.Rpt; rpt++){
		printf("[%d] rpt = %d\n", id, rpt);
		ExM.SampleUC(gV);
		TCEMBetaBinomGeneral EM(ExM.W, ExM.BoundW, ExM.GetPdUC(), ExM.alpha, gV);
		if (EM.Run()) {
			alpha += EM.alpha;
			for (int i=0; i<=ExM.W; i++) ThV[i] += EM.ThV[i];
			NSuc++;
		}
	}
	alpha /= NSuc;
	for (int i=0; i<ThV.Len(); i++) ThV[i] /= NSuc;
	printf("[%d] Experiment repeats %d times, and %d succeeded.\n", id, ExM.Rpt, NSuc);
}

void trim_tail(ExamMgr& ExM, TFltV& ThV, const double alpha){
	double minval=1, Pd = ExM.GetPdUC(), qth=0, qth_new = 0, rem = 0; // q_theta
	int Wp=1;
	for (int i=1; i<=ExM.W; i++) {
		qth += ThV[i]*TSpecFunc::BetaBinomial(0, i, Pd/alpha, (1-Pd)/alpha);
		if (ThV[i] < minval){
			minval = ThV[i];
			Wp = i;
		}
	}
	for (int i=Wp+1; i<=ExM.W; i++) {
		rem += ThV[i];
		ThV[i] = 0;
	}
	for (int i=1; i<=Wp; i++){
		ThV[i] /= 1-rem;
		qth_new += ThV[i]*TSpecFunc::BetaBinomial(0, i, Pd/alpha, (1-Pd)/alpha);
	}
	ThV[0] *= (1-qth)/(1-qth_new);
	printf("min val = %.2e   Wp=%d  rem = %.2e\n", minval, Wp, rem);
}

void em_multi(ExamMgr& ExM) {
	TExeTm tm;
	TFltV Alphas(ExM.CPU), ThVs[ExM.CPU];
	for (int i=0; i<ExM.CPU; i++) ThVs[i] = TFltV(ExM.W+1);
	std::vector<std::thread> threads;
	for (int i=0; i<ExM.CPU; i++) threads.emplace_back([i, &ExM, &Alphas, &ThVs] { em_sub(i, ExM, Alphas[i], ThVs[i]); });
	for(std::thread& t: threads) t.join();
	for (int n=1; n<ExM.CPU; n++) Alphas[0] += Alphas[n];
	Alphas[0] /= ExM.CPU;
	for (int i=0; i<=ExM.W; i++) {
		for (int n=1; n<ExM.CPU; n++) ThVs[0][i] += ThVs[n][i];
		ThVs[0][i] /= ExM.CPU;
	}
	if (ExM.TrimTail) trim_tail(ExM, ThVs[0], Alphas[0]);
	const TStr OFnm = ExM.GetBNTHFNm();
	BIO::SaveFltVWithIdx(ThVs[0], OFnm, TStr::Fmt("# Nodes: %d\n# Repeated: %d\n# Avg time cost: %.2f secs.\n# Alpha: %.6e",
			ExM.N, ExM.GetRpt(), tm.GetSecs()/ExM.GetRpt(), Alphas[0].Val));
	printf("Saved to %s\n", OFnm.CStr());
}


int main(int argc, char* argv[]){
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	const TStr GFNm = Env.GetIfArgPrefixStr("-i:", "uc.graph", "User-cont interaction graph");
	const TStr FGFNm = Env.GetIfArgPrefixStr("-f:", "fg.graph", "Follower graph");
	const int W = Env.GetIfArgPrefixInt("-w:", 10000, "W");
	const int CPU = Env.GetIfArgPrefixInt("-cpu:", std::thread::hardware_concurrency(), "# of CPUs");
	const int Rpt = Env.GetIfArgPrefixInt("-r:", 100/CPU, "Repeat");
	const int BW = Env.GetIfArgPrefixInt("-bw:", 100, "W upper bound");
	const double Pe = Env.GetIfArgPrefixFlt("-pe:", 0.1, "Edge sampling rate");
	const double Pr = Env.GetIfArgPrefixFlt("-pr:", 0.1, "Relation sampling rate");
	const double alpha = Env.GetIfArgPrefixFlt("-alpha:", 0.0001, "alpha");
	const bool TrimTail = Env.GetIfArgPrefixBool("-t:", false, "Trim tail");
	if (Env.IsEndOfRun())  return 0;

	ExamMgr ExM;
	TExeTm2 tm;
	ExM.SetSocialGraph(FGFNm).SetActionGraph(GFNm).SetW(W).SetPEdge(Pe).SetPSocial(Pr).SetRepeat(Rpt).SetCPU(CPU).IsTrimTail(TrimTail);
	em_multi(ExM);
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

