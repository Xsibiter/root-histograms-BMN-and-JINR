#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTimer.h"
#include <iostream>
#include <cmath>
#include <unordered_map>
#include <map>
#include <iostream>
#include "Math/Vector4D.h"

#define VER_STR "v45"
/*
	Версия 37:
	1. Добавил EventInfo - сведения о числе хитов и треков в dp парах
	2. Добавил гистограммы соответствующие п.1
	3. В миксинге исправил подсчет входов для корректного учета веса в гистограммах

	Версия 38:
	1. Добавил проверку совпадения числа треков в первичной вершине для каждой частицы пары в миксинге
	
	Версия 39:
	1. Праедварительный подсчет числа pd пар доступных для миксинга
	Версия не закончена. Нужно выкинуть из Ak пары, которые по числу входов в Bk ниже порога на ratio!!!

	Версия 44:
	1. Много изменений. Переделан алгоритм:
		(1) найти пару,
		(2) подобрать протоны для миксинга с дейтроном из этой пары,
		(3) заполнить B(k) с весом по числу найденных протонов
		(4) заполнить A(k)
	2. Добавил картинки для разных мишеней
	3. Картинки для изучения особенности на dPLab_vs_k
	
	Версия 45:
	1. Изменил диапазоны Po было 2.5 и 5.0 GeV, стало 3 и 6 GeV
	
	ИНТЕРВАЛЫ ПО NTPV 0-11 12-14 15-50
	
*/
/// Индексы в массиве ядер мишени соответствуют нумерации во входном файле
const char * gTargets[] = {"empty","Al","C","Cu","Pb","Sn"};

#define mass2proton 0.880354
#define mass2deutron 3.51938

#define SigmaCriteria 2.5

#define PTMIN 0.0
#define PTMAX 1.4
#define YLMIN 0.0
#define YLMAX 3.0

UInt_t gNPVTMin = 0;
UInt_t gNPVTMax = 11;

Int_t gPoRangeId = 0;
const char * gCutTitle[3][2] = {
	{"",""},
	{"_PltPo",", P_{p}<3, P_{d}<6"},
	{"_PgtPo",", P_{p}#geq3, P_{d}#geq6"},
};

Bool_t CheckPoRange(Double_t Pp, Double_t Pd) {
	const Double_t Pp_lim = 3.0; // для протонов
	const Double_t Pd_lim = 6.0; // для дейтронов
	if (gPoRangeId == 0) {
		return kTRUE;
	} else if (gPoRangeId == 1) {
		if (Pp < Pp_lim && Pd < Pd_lim) {
			return kTRUE;
		}
	} else if (gPoRangeId == 2) { 
		if (Pp >= Pp_lim && Pd >= Pd_lim) {
			return kTRUE;
		}
	}
	return kFALSE;
}

const char * GetName(const char * name, const char * suffix = "") {
	static TString str;
	str = TString::Format("%s%s",name,suffix);
	return str.Data();
}

const char * GetTitle(const char * format, const char * name, Int_t suffix_index = 1) {
	static TString str;
	str = TString::Format(format,name,gCutTitle[gPoRangeId][suffix_index]);
	return str.Data();
}

/*
	2D k vs p, 2D k/k_mix vs p/p_mix, средний импульс d и p при k < 0.3
*/
struct Index {
	UInt_t run;
	UInt_t event;
	friend bool operator<(Index a, Index b) {
		return a.run < b.run || (a.run == b.run && a.event < b.event);
	}
	friend bool operator==(Index a, Index b) {
		return a.run == b.run && a.event == b.event;
	}
};

struct Event {
	Double_t pxfr;
	Double_t pyfr;
	Double_t pzfr;
	Double_t m2fr;
	Double_t Pt;
	Double_t Y;
	Double_t GetP() {return TMath::Sqrt(pxfr*pxfr+pyfr*pyfr+pzfr*pzfr);}
};

struct EventInfo {
	UInt_t target;
	UInt_t trigger;
	UInt_t NhitsSI;
	UInt_t NhitsGem;
	UInt_t PVNtracks;
	UInt_t Npd;
	UInt_t ratio;
};


TMap fTHMap;
void AddTH(TH1 * p) {
	fTHMap.Add(p, p);
}
void FillTH(const char * name, Double_t x, Double_t weight) {
	TH1D * pTH1D = (TH1D *) (fTHMap.GetValue(name));
	if (pTH1D == 0) {
		Fatal(__func__,"Use AddTH(\"%s\",...) first!\n",name);
	}
	pTH1D->Fill(x, weight);
}
void FillTH(const char * name, Double_t x, Double_t y, Double_t weight) {
	TH2D * pTH2D = (TH2D *) (fTHMap.GetValue(name));
	if (pTH2D == 0) {
		Fatal(__func__,"Use AddTH(\"%s\",...) first!\n",name);
	}
	pTH2D->Fill(x, y, weight);
}
void WriteTH() {
	TIterator * it = fTHMap.MakeIterator();
	TObject * p = 0;
	while ((p = it->Next())) {
		((TH1 *) p)->Write();
	}
}

struct THD {
	TH1D * k_All;
	TH1D * k_Cut;
	TH1D * dPLab_All;
	TH1D * dPLab_Cut;	
	TH2D * dPLab_vs_k;
	void Init(const char * name, const char * title) {
		k_All = new TH1D(GetName(name,"_All"),GetTitle("%s%s;GeV/c",title), 150,0,1.5);
		k_Cut = new TH1D(GetName(name,"_Cut"),GetTitle("%s, P_{lab}>0.4%s;GeV/c",title), 150,0,1.5);
		dPLab_All = new TH1D(GetName(name,"_dPLab_All"),GetTitle("%s dP_{lab}, no cut%s;dP_{lab}, Gev/c",title),150,0,10);
		dPLab_Cut = new TH1D(GetName(name,"_dPLab_Cut"),GetTitle("%s dP_{lab}, k>0.2%s;dP_{lab}, GeV/c",title),150,0,10);	
		dPLab_vs_k = new TH2D(GetName(name,"_dPLab_vs_k"),GetTitle("%s%s;k, Gev/c;dP_{lab}, GeV/c",title),150,0,5,150,0,15);
	}
	void Fill(Double_t k, Double_t dPLab, Double_t weight = 1) {
		k_All->Fill(k,weight);
		if (dPLab > 0.4) k_Cut->Fill(k,weight);
		dPLab_All->Fill(dPLab,weight);
		if (2.*k > 0.4) dPLab_Cut->Fill(dPLab,weight); // k is dPSCM/2
		dPLab_vs_k->Fill(k,dPLab);
	}
	void Write() {
		k_All->Write();
		k_Cut->Write();
		dPLab_All->Write();
		dPLab_Cut->Write();	
		dPLab_vs_k->Write();
	}
	THD * Clone(const char * name, const char * title) {
		THD * clone = new THD();
		clone->k_All = (TH1D*)k_All->Clone(GetName(name,"_All"));
		clone->k_All->SetTitle(GetTitle("%s%s;GeV/c",title));
		clone->k_Cut = (TH1D*)k_Cut->Clone(GetName(name,"_Cut"));
		clone->k_Cut->SetTitle(GetTitle("%s, P_{lab}>0.4%s;GeV/c",title));
		clone->dPLab_All = (TH1D*)dPLab_All->Clone(GetName(name,"_dPLab_All"));
		clone->dPLab_All->SetTitle(GetTitle("%s dP_{lab}, no cut%s;dP_{lab}, GeV/c",title));
		clone->dPLab_Cut = (TH1D*)dPLab_Cut->Clone(GetName(name,"_dPLab_Cut"));
		clone->dPLab_Cut->SetTitle(GetTitle("%s dP_{lab}, k>0.2%s;dP_{lab}, GeV/c",title));
		clone->dPLab_vs_k = (TH2D*)dPLab_vs_k->Clone(GetName(name,"_dPLab_vs_k"));
		clone->dPLab_vs_k->SetTitle(GetTitle("%s%s;k, GeV/c;dP_{lab}, GeV/c",title));
		return clone;
	}
	void Divide(THD * other) {
		k_All->Divide(other->k_All);
		k_Cut->Divide(other->k_Cut);
		dPLab_All->Divide(other->dPLab_All);
		dPLab_Cut->Divide(other->dPLab_Cut);
		dPLab_vs_k->Divide(other->dPLab_vs_k);
	}
};

bool in_window(double variable, double min, double max) {
	bool pass = false;
	if (variable > min && variable <= max) pass = true;
	return pass;
}

multimap<Index,  Event> gProtons;
multimap<Index,  Event> gDeutrons;
map<Index, EventInfo> gEventInfo;

void ReadRoot() {
	TFile * file = TFile::Open("/home/pnaleks/femtoscopy/Fragments_Ar_full_statistic.root");
	TTree * tree = (TTree*) file->Get("tfr");

	Double_t m2fr=0., pfr=0.,Dch_dx=0., Dch_dy=0., Dch_SigmaX=0., Dch_SigmaY=0.,Tof700_dx=0., Tof700_dy=0., Tof700_SigmaX=0., Tof700_SigmaY=0., pvx=0., pvy=0., pvz=0.,dyfr=0.,ptfr=0., prtnyfr=0., yfr=0, x0pv=0, y0pv=0, z0pv=0, piyfr=0, pxfr=0, pyfr=0, pzfr =0, bfr=0;
	UInt_t  NhitsSI=0., NhitsGem=0.,nurun=0.,nuev=0.,PVNtracks=0., target=0., trigger=0.;
	
	tree->SetBranchAddress("m2fr", &m2fr);		// квадрат массы
	tree->SetBranchAddress("pfr", &pfr);		// импульс
	tree->SetBranchAddress("ptfr", &ptfr);		// поперечная компонента импульса
	tree->SetBranchAddress("dyfr", &dyfr);		// быстрота дейтона
	tree->SetBranchAddress("prtnyfr", &prtnyfr);// быстрота протона
	tree->SetBranchAddress("piyfr", &piyfr);	// быстрота пиона
	tree->SetBranchAddress("yfr", &yfr);		// быстрота
	tree->SetBranchAddress("Dch_dx", &Dch_dx);	// ресидуал хита от трека
	tree->SetBranchAddress("Dch_dy", &Dch_dy);
	tree->SetBranchAddress("Dch_SigmaX", &Dch_SigmaX);	// ширина распределения ресидуалов от быстроты?
	tree->SetBranchAddress("Dch_SigmaY", &Dch_SigmaY);
	tree->SetBranchAddress("Tof700_dx", &Tof700_dx);
	tree->SetBranchAddress("Tof700_dy", &Tof700_dy);       
	tree->SetBranchAddress("Tof700_SigmaX", &Tof700_SigmaX);
	tree->SetBranchAddress("Tof700_SigmaY", &Tof700_SigmaY);  
	tree->SetBranchAddress("NhitsSI", &NhitsSI);		// хитов в силиконах
	tree->SetBranchAddress("NhitsGem", &NhitsGem); 		// хиты в джемах
	tree->SetBranchAddress("nurun", &nurun);  			// номер рана
	tree->SetBranchAddress("nuev", &nuev); 				// номер события
	tree->SetBranchAddress("PVNtracks", &PVNtracks); 	// число треков в первичной вершине
	tree->SetBranchAddress("target", &target); 			// номер мишени
	tree->SetBranchAddress("trigger", &trigger); 		// номер триггера
	tree->SetBranchAddress("pvx", &pvx);				// восстановленная вершина
	tree->SetBranchAddress("pvy", &pvy);
	tree->SetBranchAddress("pvz", &pvz);
	tree->SetBranchAddress("x0pv", &x0pv);				// экстраполяция трека в z вершины
	tree->SetBranchAddress("y0pv", &y0pv);
	tree->SetBranchAddress("z0pv", &z0pv);

	tree->SetBranchAddress("pxfr", &pxfr);				// компоненты импульса
	tree->SetBranchAddress("pyfr", &pyfr);
	tree->SetBranchAddress("pzfr", &pzfr);
	tree->SetBranchAddress("bfr", &bfr);				// бета

	AddTH( new TH1D("P_p_all","Protons momentum;P, GeV/c",150,0,15) );
	AddTH( new TH1D("P_d_all","Deutrons momentum;P, GeV/c",150,0,15) );

	// Pt vs rapidity Pt GeV/c
	AddTH( new TH2D("PtVsY_p_all","All p^{+};rapidity;P_{t}, GeV/c",200,-1,3.5,200,0,4) );
	AddTH( new TH2D("PtVsY_d_all","All d^{+};rapidity;P_{t}, GeV/c",200,-1,3.5,200,0,4) );

	AddTH( new TH2D("target_hit","Number of reconstructed particles;Target number; particle type (1 is p^{+}, 2 is d^{+})",8,0,8,4,0,4) );

	AddTH( new TH2D("NTPV_Vs_Pp","Num P.V.T. vs Pp;P, GeV/c;NTPV",200,0,10,50,0,50) );


	Double_t NumEv  = tree->GetEntries();
	Int_t print_cnt = 0;
	UInt_t ev_cnt = 0;
	for (int ev = 0;  ev < NumEv; ev++) {
		if (++ev_cnt % 100000 == 0) cout << "ReadRoot: " << ev_cnt << endl;

		tree->GetEntry(ev);
		if (PVNtracks < 2) continue;

		if (pvx < -2.0 || pvx > 3.0 || pvy < 1.0 || pvy > 5.0 || pvz < -1.0 || pvz > 2.5) continue;
		if (x0pv < -2.0 || x0pv > 3.0 || y0pv <  0.5 || y0pv > 4.5) continue;
		
		bool matched =
			Dch_dx < SigmaCriteria * Dch_SigmaX && Dch_dx > -SigmaCriteria * Dch_SigmaX && 
			Dch_dy < SigmaCriteria * Dch_SigmaY && Dch_dy > -SigmaCriteria * Dch_SigmaY && 
			Tof700_dx < SigmaCriteria * Tof700_SigmaX && Tof700_dx > -SigmaCriteria * Tof700_SigmaX && 
			Tof700_dy < SigmaCriteria * Tof700_SigmaY && Tof700_dy > -SigmaCriteria * Tof700_SigmaY && 
			NhitsSI > 0 && 
			NhitsGem > 3;
		
		if (!matched) continue;
		
		if ( in_window(m2fr, 2.4 , 6) ) {
			Index index = {nurun, nuev};
			Event event = {pxfr, pyfr, pzfr,m2fr,ptfr,yfr};
			gDeutrons.insert(make_pair(index, event));
			FillTH("P_d_all", event.GetP(), 1);
			FillTH("PtVsY_d_all", event.Y, event.Pt, 1);
			FillTH("target_hit",target,2,1);
			EventInfo info = {target,trigger,NhitsSI,NhitsGem,PVNtracks,0,0};
			gEventInfo[index] = info;
		}

		if ( in_window(m2fr, 0.5 , 1.5) ) {
			Index index = {nurun, nuev};
			Event event = {pxfr, pyfr, pzfr,m2fr,ptfr,yfr};
			gProtons.insert(make_pair(index, event));
			FillTH("P_p_all", event.GetP(), 1);
			FillTH("PtVsY_p_all", event.Y, event.Pt, 1);
			FillTH("target_hit",target,1,1);
			EventInfo info = {target,trigger,NhitsSI,NhitsGem,PVNtracks,0,0};
			gEventInfo[index] = info;
			FillTH("NTPV_Vs_Pp",event.GetP(),PVNtracks,1);
		}
		
	}
	//file->Close();
	cout << "ReadRoot done: " << ev_cnt << " entries" << endl;
}

THD Ak, Bk;
multimap<Index,  Event>::iterator gMixIterator;
Bool_t Mix(Event &ev_d, multimap<Index,  Event>::iterator &it_p) {
	const Int_t mix_max = 1000;
	const Int_t mix_min = 100;
	
	auto first_it = gMixIterator;
	multimap<Index,  Event> mix;
	do {
		if (gMixIterator != it_p) {
			Index index = gMixIterator->first;
			if (gDeutrons.count(index)>0) {
				Event ev_p = gMixIterator->second;
				if (CheckPoRange(ev_p.GetP(), ev_d.GetP())) {
					mix.insert(make_pair(index, gMixIterator->second));
				}					
			}				
		}
		gMixIterator++;
		if (gMixIterator == gProtons.end()) gMixIterator = gProtons.begin();
	} while(gMixIterator != first_it && mix.size() < mix_max);
	if (mix.size() < mix_min) return kFALSE;
	
	for (auto it = mix.begin(); it != mix.end(); it++) { 
		Event ev_p = it->second;		
		// Define the masses of the particles
		Double_t m_d = std::sqrt(ev_d.m2fr);
		Double_t m_p = std::sqrt(ev_p.m2fr);
		// Create 4-momentum vectors for each particle in the lab frame
		ROOT::Math::PxPyPzMVector p_d(ev_d.pxfr, ev_d.pyfr, ev_d.pzfr, m_d);
		ROOT::Math::PxPyPzMVector p_p(ev_p.pxfr, ev_p.pyfr, ev_p.pzfr, m_p);
		// Create the total 4-momentum in the lab frame
		ROOT::Math::PxPyPzMVector p_total = p_d + p_p;
		ROOT::Math::PxPyPzMVector dPLab = p_d - p_p;
		// Boost to the center-of-mass frame
		ROOT::Math::Boost to_cm(p_total.BoostToCM());
		p_d = to_cm(p_d);
		p_p = to_cm(p_p);
		
		Bk.Fill(p_d.P(),dPLab.P(),1./mix.size());
	}
	return kTRUE;
}

int femtoscopy(Int_t cut_index = 0) {
	cout << "INITIALIZE" << endl;
	gPoRangeId = cut_index;
	ReadRoot();

#define ANALIZE
#ifdef ANALIZE
	Ak.Init("Ak","A(k)");
	Bk.Init("Bk","B(k)");
	
	AddTH( new TH1D("Ak_NPVT_00_11", "Ak_NPVT_00_11", 150, 0, 1.5) );
	AddTH( new TH1D("Ak_NPVT_12_14", "Ak_NPVT_12_14", 150, 0, 1.5) );
	AddTH( new TH1D("Ak_NPVT_15_50", "Ak_NPVT_15_50", 150, 0, 1.5) );

	AddTH( new TH1D("Pp","Protons momentum (pd events);P_{p}, GeV/c",150,0,15) );
	AddTH( new TH1D("Pd","Deutrons momentum (pd events);P_{d}, GeV/c",150,0,15) );
	AddTH( new TH1D("Pp_dCut","Protons momentum (P_{d}>7.6GeV/c ~20%);P_{p}, GeV/c",150,0,15) );
	AddTH( new TH1D("Pd_pCut","Deutrons momentum (P_{p}>4.2GeV/c ~20%);P_{d}, GeV/c",150,0,15) );

	AddTH( new TH2D("PpPd",";P_{d}, GeV/c;P_{p}, GeV/c",160,0,14,130,0,8) );
	AddTH( new TH2D("PpPd_Lower",";P_{d}, GeV/c;P_{p}, GeV/c",160,0,14,130,0,8) );

	AddTH( new TH2D("NdNp","Number of d^{+} and p^{+} in the event;N_{p};N_{d}",7,0,7,5,0,5) );
	
	AddTH(new TH2D("PtVsY_d","Deutrons from pd events;rapidity;P_{t}, GeV/c",200,-1,3.5,200,0,4) );
	for (Int_t i = 1; i < sizeof(gTargets)/sizeof(const char *); i++) {
		TString name = TString::Format("PtVsY_d_Ar%s",gTargets[i]);
		TString title = TString::Format("Deutrons from pd events, %s target);rapidity;P_{t}, GeV/c",gTargets[i]);
		AddTH(new TH2D(name,title,200,-1,3.5,200,0,4));
	}
	
	AddTH(new TH2D("PtVsY_p","Protons from pd events;rapidity;P_{t}, GeV/c",200,-1,3.5,200,0,4) );
	for (Int_t i = 1; i < sizeof(gTargets)/sizeof(const char *); i++) {
		TString name = TString::Format("PtVsY_p_Ar%s",gTargets[i]);
		TString title = TString::Format("Protons from pd events, %s target);rapidity;P_{t}, GeV/c",gTargets[i]);
		AddTH(new TH2D(name,title,200,-1,3.5,200,0,4));
	}
	AddTH(new TH2D("PtVsY_p_Upper","Protons from pd events (upper);rapidity;P_{t}, GeV/c",200,-1,3.5,200,0,4) );
	AddTH(new TH2D("PtVsY_p_Lower","Protons from pd events (lower);rapidity;P_{t}, GeV/c",200,-1,3.5,200,0,4) );
	AddTH(new TH2D("PtVsY_d_Upper","Deutrons from pd events (upper);rapidity;P_{t}, GeV/c",200,-1,3.5,200,0,4) );
	AddTH(new TH2D("PtVsY_d_Lower","Deutrons from pd events (lower);rapidity;P_{t}, GeV/c",200,-1,3.5,200,0,4) );
	AddTH(new TH2D("Ak_dPLab_vs_k_Lower",GetTitle("%s (Lower)%s;k, Gev/c;dP_{lab}, GeV/c","A(k) "),150,0,5,150,0,15) );
	AddTH(new TH2D("YdYp_Upper","PtPt (upper);d rapidity;p rapidity",200,-1,3.5,200,-1,3.5) );
	AddTH(new TH2D("YdYp_Lower","PtPt (lower);d rapidity;p rapidity",200,-1,3.5,200,-1,3.5) );

	TH1D * hAk_Target = new TH1D("hAk_Target","hAk_Target",10,0,10);
	TH1D * hAk_NhitsSI = new TH1D("hAk_NhitsSI","hAk_NhitsSI",10,0,10);
	TH1D * hAk_NhitsGem = new TH1D("hAk_NhitsGem","hAk_NhitsGem",10,0,10);
	TH1D * hAk_PVNtracks = new TH1D("hAk_PVNtracks","hAk_PVNtracks",50,0,50);

	cout << "ANALIZE" << endl;
	UInt_t cnt_Ak = 0;
	UInt_t ev_cnt = 0;
	
	////////////////////////////////////////////////////////////////////////////
	// Заполняю NdNp для Nd == 0 && Np > 0
	for (auto it = gProtons.begin(); it != gProtons.end(); it = gProtons.upper_bound(it->first)) {
		Index index = it->first;
		Int_t p_hit = gProtons.count(index);
		Int_t d_hit = gDeutrons.count(index);
		if (d_hit == 0) {
			FillTH("NdNp",p_hit,d_hit,1);
		}
	}
	////////////////////////////////////////////////////////////////////////////
	// Основной цикл
	gMixIterator = gProtons.begin();
	for (auto it = gDeutrons.begin(); it != gDeutrons.end(); it = gDeutrons.upper_bound(it->first)) {
		if (++ev_cnt % 1000 == 0) cout << "ANALIZE: " << ev_cnt << endl;
		Index index = it->first;
		
		Int_t p_hit = gProtons.count(index);
		Int_t d_hit = gDeutrons.count(index);
		// Заполняю NdNp для Nd > 0 && Np >= 0
		FillTH("NdNp",p_hit,d_hit,1);
		
		if (!p_hit || !d_hit) continue;

		auto proton_range = gProtons.equal_range(index);
		auto deutron_range = gDeutrons.equal_range(index);
		
		gEventInfo[index].Npd = 0;
		auto event_info_copy = gEventInfo[index];
		
		for (auto p_it = proton_range.first; p_it != proton_range.second; p_it++) {
			for (auto d_it = deutron_range.first; d_it != deutron_range.second; d_it++) {

				Event proton = p_it->second;
				Event deutron = d_it->second;
				
				if (CheckPoRange(proton.GetP(), deutron.GetP())) {
					if (Mix(deutron,p_it)) {
						// Define the masses of the particles
						Double_t m_proton = std::sqrt(proton.m2fr);
						Double_t m_deutron = std::sqrt(deutron.m2fr);
						// Create 4-momentum vectors for each particle in the lab frame
						ROOT::Math::PxPyPzMVector p_proton(proton.pxfr, proton.pyfr, proton.pzfr, m_proton);
						ROOT::Math::PxPyPzMVector p_deutron(deutron.pxfr, deutron.pyfr, deutron.pzfr, m_deutron);
						FillTH("Pp",p_proton.P(),1);
						FillTH("Pd",p_deutron.P(),1);
						FillTH("PpPd",p_deutron.P(),p_proton.P(),1);
						if (p_proton.P()>4.2) FillTH("Pd_pCut",p_deutron.P(),1);
						if (p_deutron.P()>7.6) FillTH("Pp_dCut",p_proton.P(),1);
						// Create the total 4-momentum in the lab frame
						ROOT::Math::PxPyPzMVector p_total = p_proton + p_deutron;
						ROOT::Math::PxPyPzMVector dPLab = p_proton - p_deutron;
						// Boost to the center-of-mass frame
						ROOT::Math::Boost to_cm(p_total.BoostToCM());
						ROOT::Math::PxPyPzMVector k_proton = to_cm(p_proton);
						ROOT::Math::PxPyPzMVector k_deutron = to_cm(p_deutron);

						Ak.Fill(k_deutron.P(),dPLab.P(),1);
						
						// Double_t cut_val = 10./3.* k_deutron.P() - 1.;
						Double_t cut_val = 4.* k_deutron.P();
						if (dPLab.P() < cut_val) {
							FillTH("Ak_dPLab_vs_k_Lower",k_deutron.P(),dPLab.P());
							FillTH("PpPd_Lower",p_deutron.P(),p_proton.P(),1);
							FillTH("PtVsY_p_Lower",proton.Y,proton.Pt,1);
							FillTH("PtVsY_d_Lower",deutron.Y,deutron.Pt,1);
							FillTH("YdYp_Lower",deutron.Y,proton.Y,1);
						} else {
							FillTH("PtVsY_p_Upper",proton.Y,proton.Pt,1);
							FillTH("PtVsY_d_Upper",deutron.Y,deutron.Pt,1);
							FillTH("YdYp_Upper",deutron.Y,proton.Y,1);
						}

						cnt_Ak++;
						gEventInfo[index].Npd++;
						if (event_info_copy.PVNtracks < 12) {
							FillTH("Ak_NPVT_00_11",k_deutron.P(),1);
						} else if (event_info_copy.PVNtracks < 15) {
							FillTH("Ak_NPVT_12_14",k_deutron.P(),1);
						} else {
							FillTH("Ak_NPVT_15_50",k_deutron.P(),1);
						}
						FillTH("PtVsY_p",proton.Y,proton.Pt,1);
						FillTH("PtVsY_d",deutron.Y,deutron.Pt,1);
						
						hAk_Target->Fill(event_info_copy.target);
						hAk_NhitsSI->Fill(event_info_copy.NhitsSI);
						hAk_NhitsGem->Fill(event_info_copy.NhitsGem);
						hAk_PVNtracks->Fill(event_info_copy.PVNtracks);
						
						TString name1 = TString::Format("PtVsY_d_Ar%s",gTargets[event_info_copy.target]);
						TString name2 = TString::Format("PtVsY_p_Ar%s",gTargets[event_info_copy.target]);
						FillTH(name1,deutron.Y,deutron.Pt,1);
						FillTH(name2,proton.Y,proton.Pt,1);
					}
				}								
			}
		}
	}
	cout << "DONE" << endl;

#endif	
	TFile *fReco = new TFile(GetTitle("femtoscopy_%s%s.root",VER_STR,0), "RECREATE");
	cout << "Writing into \"" << fReco->GetName() << "\"\n";

	WriteTH();

#ifdef ANALIZE
	fReco->cd();
	THD * Ck = Ak.Clone("Ck","C(k)");
	Ck->Divide(&Bk);

	Ak.Write();
	Bk.Write();
	Ck->Write();
	
	hAk_Target->Write();
	hAk_NhitsSI->Write();
	hAk_NhitsGem->Write();
	hAk_PVNtracks->Write();
#endif
	fReco->Write();
	fReco->Close();

	return 0;
}
// Fit func exp([0]*x+[1])*[2]+[3]
// широкие кореляции 
// dP_ЛАБ vs dP_СЦМ
