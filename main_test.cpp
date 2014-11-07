/*
 * main.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: jzzhao
 */
#include "stdafx.h"
#include "ExamMgr.h"

const TStr Root = "/home/jzzhao/workspace/trc/enron/cleaned/";

/**
 * count triangles per week
 */
void count_triangles(){
	PNEGraph G = PNEGraph::TObj::New();
	int pre_week = 0, ntrid = 0;
	TIntPrV TridCnt; TIntV NTridsPerWeek;
	TSsParser Ss(Root+"enron_add_week_no.gz");
	while(Ss.Next()){
		const int fid = Ss.GetInt(0), tid = Ss.GetInt(1), week = Ss.GetInt(2);
		if(week!=pre_week && G->GetEdges()!=0) {
			TSnap::GetTriadParticipAll(G, TridCnt);
			for (int i=0; i<TridCnt.Len(); i++) ntrid += TridCnt[i].Val1*TridCnt[i].Val2;
			NTridsPerWeek.Add(ntrid/3);
			G->Clr(); TridCnt.Clr(); pre_week = week; ntrid = 0;
		}
		if (!G->IsNode(fid)) G->AddNode(fid);
		if (!G->IsNode(tid)) G->AddNode(tid);
		G->AddEdge(fid, tid);
	}
	BIO::SaveIntVWithIdx(NTridsPerWeek, Root+"ntrids_week.dat");
}

/**
 * triadic cardinality distribution
 */
void dist_triangles(const int WK){
	PNEGraph G = PNEGraph::TObj::New();
	TSsParser Ss(Root+"enron_cleaned.gz");
	while(Ss.Next()){
		const int fid = Ss.GetInt(0), tid = Ss.GetInt(1), week = Ss.GetInt(2);
		if(week==WK) {
			if (!G->IsNode(fid)) G->AddNode(fid);
			if (!G->IsNode(tid)) G->AddNode(tid);
			G->AddEdge(fid, tid);
		}else if(week>WK) break;
	}
	TSnap::SaveEdgeList(G, Root+TStr::Fmt("enron_%d.gz", WK));
	return;
	TIntPrV TriadCntV;
	TSnap::GetTriadParticipAll(G, TriadCntV);
	int N = G->GetNodes();
	FILE* fw=fopen((Root+TStr::Fmt("tcd_%d.dist", WK)).CStr(), "w");
	fprintf(fw, "# Nodes: %d\n", N);
	for (int i=0; i<TriadCntV.Len(); i++) {
		int card = TriadCntV[i].Val1;
		int freq = TriadCntV[i].Val2;
		double prob = freq/(double)N;
		fprintf(fw, "%d\t%d\t%.6e\n", card, freq, prob);
	}
	fclose(fw);
}

/**
 * degree distribution
 */
void dist_degree(const int WK){
	PNEGraph G = PNEGraph::TObj::New();
	TStr Root = "/media/WinE/workspace/triadic_cardinality/enron/";
	TSsParser Ss(Root+"enron_add_week_no.gz");
	while(Ss.Next()){
		const int week = Ss.GetInt(0), fid = Ss.GetInt(1), tid = Ss.GetInt(2);
		if(week==WK) {
			if (!G->IsNode(fid)) G->AddNode(fid);
			if (!G->IsNode(tid)) G->AddNode(tid);
			G->AddEdge(fid, tid);
		}else if(week>WK) break;
	}
	TIntPrV DegV;
	TSnap::GetDegCnt<PNEGraph>(G, DegV);
	BIO::SaveIntPrV(DegV, Root+TStr::Fmt("DegCnt_%d.dat", WK));
}

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

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
const int SID = 100000, N = 22477;
////////////////////////////////////////////////////////

int add_spamming_rnd(const TIntV& Nodes, const PNEGraph& OriG, const int NSpam, TIntPrV& TriadCntV){
	PNEGraph G = PNEGraph::New();
	for (TNEGraph::TEdgeI EI = OriG->BegEI(); EI < OriG->EndEI(); EI++) {
		if (!G->IsNode(EI.GetSrcNId())) G->AddNode(EI.GetSrcNId());
		if (!G->IsNode(EI.GetDstNId())) G->AddNode(EI.GetDstNId());
		G->AddEdge(EI);
	}
	G->AddNode(SID);

	TIntV targets;
	TRandom::ChooseWithReplacement(Nodes, targets, NSpam);
	for (int i=0; i<NSpam; i++){
		if (!G->IsNode(targets[i])) G->AddNode(targets[i]);
		G->AddEdge(SID, targets[i]);
	}
	printf("nspam: %d, nodes: %d, edges: %d\n", NSpam, G->GetNodes(), G->GetEdges());

	TriadCntV.Clr();
	TSnap::GetTriadParticipAll(G, TriadCntV);
	/*
	FILE* fw=fopen((Root+TStr::Fmt("enron_22_SP%dK.dist", NSpam/1000)).CStr(), "w");
	fprintf(fw, "# Nodes: %d\n", N);
	for (int i=0; i<TriadCntV.Len(); i++) {
		int card = TriadCntV[i].Val1;
		int freq = TriadCntV[i].Val2;
		if (card == 0) freq += (N - G->GetNodes());
		double prob = freq/(double)N;
		fprintf(fw, "%d\t%d\t%.6e\n", card, freq, prob);
	}
	fclose(fw);
	*/
	return N - G->GetNodes();
}

double KL(const TIntPrV& PV, const TIntH& QH){
	int k;
	double p, q, rst=0;
	for (int i=0; i<PV.Len(); i++) {
		k = PV[i].Val1;
		p = PV[i].Val2/(double)N;
		if (!QH.IsKey(k)) q = 1E-10;
		else q = QH(k)/(double)N;
		rst += p*log(p/q);
	}
	return rst;
}

void spam_rnd(const int R){
	TIntPrV PV, QV;
	TIntH QH;
	TIntFltKdV klV;
	TIntV Nodes;
	BIO::LoadIntV(Root+"nodes.gz", Nodes);
	BIO::LoadIntPrV(Root+"tcd_22.dist", PV);
	PNEGraph OriG = TSnap::LoadEdgeList<PNEGraph>(Root+"enron_22.gz");
	for (int i=1; i<=10; i++){
		const int NSpam = i*1000;
		double kl = 0;
		for (int r=0; r<R; r++){
			const int delta = add_spamming_rnd(Nodes, OriG, NSpam, QV);
			QH.Clr();
			for (int i=0; i<QV.Len(); i++) QH(QV[i].Val1) = QV[i].Val2;
			QH[0] += delta;
			kl += KL(PV, QH);
		}
		klV.Add(TIntFltKd(i, kl/R));
	}
	BIO::SaveIntFltKdV(klV, Root+"KL_rnd.dat");
}

int add_spamming_friend(const PUNGraph& FG, const PNEGraph& OriG, const int NSpam, TIntPrV& TriadCntV){
	PNEGraph G = PNEGraph::New();
	for (TNEGraph::TEdgeI EI = OriG->BegEI(); EI < OriG->EndEI(); EI++) {
		if (!G->IsNode(EI.GetSrcNId())) G->AddNode(EI.GetSrcNId());
		if (!G->IsNode(EI.GetDstNId())) G->AddNode(EI.GetDstNId());
		G->AddEdge(EI);
	}
	G->AddNode(SID);
	for (int i=0; i<NSpam/2; i++){
		const int tar1 = FG->GetRndNId();
		TUNGraph::TNodeI NI = FG->GetNI(tar1);
		const int tar2 = NI.GetDeg()==0 ? tar1 : NI.GetRndNbhNId();
		if (!G->IsNode(tar1)) G->AddNode(tar1);
		if (!G->IsNode(tar2)) G->AddNode(tar2);
		G->AddEdge(SID, tar1);
		G->AddEdge(SID, tar2);
	}
	printf("nspam: %d, nodes: %d, edges: %d\n", NSpam, G->GetNodes(), G->GetEdges());

	TriadCntV.Clr();
	TSnap::GetTriadParticipAll(G, TriadCntV);
/*
	FILE* fw=fopen((Root+TStr::Fmt("enron_22_SP%dK_fnd.dist", NSpam/1000)).CStr(), "w");
	fprintf(fw, "# Nodes: %d\n", N);
	for (int i=0; i<TriadCntV.Len(); i++) {
		int card = TriadCntV[i].Val1;
		int freq = TriadCntV[i].Val2;
		if (card == 0) freq += (N - G->GetNodes());
		double prob = freq/(double)N;
		fprintf(fw, "%d\t%d\t%.6e\n", card, freq, prob);
	}
	fclose(fw);
*/
	return N - G->GetNodes();
}

void spam_friend(const int R){
	TIntPrV PV, QV;
	TIntH QH;
	TIntFltKdV klV;
	BIO::LoadIntPrV(Root+"tcd_22.dist", PV);
	PUNGraph FriendG = TSnap::LoadEdgeList<PUNGraph>(Root+"enron_cleaned.gz");
	PNEGraph OriG = TSnap::LoadEdgeList<PNEGraph>(Root+"enron_22.gz");
	for (int i=1; i<=10; i++){
		const int NSpam = i*1000;
		double kl = 0;
		for (int r=0; r<R; r++){
			const int delta = add_spamming_friend(FriendG, OriG, NSpam, QV);
			QH.Clr();
			for (int i=0; i<QV.Len(); i++) QH(QV[i].Val1) = QV[i].Val2;
			QH[0] += delta;
			kl += KL(PV, QH);
		}
		klV.Add(TIntFltKd(i, kl/R));
	}
	BIO::SaveIntFltKdV(klV, Root+"KL_friend.dat");
}


int main(int argc, char* argv[]){
//	dist_triangles(22);
//	spam_rnd(10);
//	spam_friend(10);
//	return 0;

	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	const TStr GFNm = Env.GetIfArgPrefixStr("-i:", "test.graph", "Input graph");
	const TStr FGFNm = Env.GetIfArgPrefixStr("-f:", "test.graph", "Follower graph");
	const int W = Env.GetIfArgPrefixInt("-w:", 10000, "W");
	const int CPU = Env.GetIfArgPrefixInt("-n:", 8, "Cores to use, max=8");
	const int Rpt = Env.GetIfArgPrefixInt("-r:", 12, "Repeat");
	const double Pe = Env.GetIfArgPrefixFlt("-p:", 0.1, "Edge sampling rate");
	const double Pr = Env.GetIfArgPrefixFlt("-q:", 0.1, "Relation sampling rate");
	const bool TrimTail = Env.GetIfArgPrefixBool("-t:", false, "Trim tail");
	if (Env.IsEndOfRun()) return 0;

	TExeTm tm;
//	ExamMgr ExM(GFNm, CPU, W, Pe, Rpt);
	ExamMgr ExM(GFNm, FGFNm, CPU, W, Pe, Pr, Rpt, TrimTail);
	TIntPrV gV;
	ExM.SampleUC(gV);
	printf("Cost time: %s.\n", tm.GetStr());
	return 0;
}

