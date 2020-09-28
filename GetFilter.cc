#include "GetFilter.h"

using namespace std;

bool GetFilter(int argc, char** argv)
{
	double Scale, Threshold, Charge_Thr;// rise time 1525: 0.55; 1305: 0.4
	int Nsize, WavNum, step, type, RangeA, RangeB;
	if (argc != 12) {
		cout << "Wrong number (" << argc << ") of input arguments" << endl;
		return false;
	}
	sscanf(argv[1], "%d", &step);
	sscanf(argv[2], "%d", &type);
	TString name(argv[3]);
	sscanf(argv[4], "%d", &RangeA);
	sscanf(argv[5], "%d", &RangeB);
	sscanf(argv[6], "%lf", &Scale);
	sscanf(argv[7], "%d", &Nsize);
	sscanf(argv[8], "%d", &WavNum);
	sscanf(argv[9], "%lf", &Threshold);
	sscanf(argv[10], "%lf", &Charge_Thr);
	TString basename(argv[11]);
	ofstream out(name + "_Filter.txt", std::ofstream::out);
	TString ro(".root");
	TString tx(".txt");
	ifstream in_bl(argv[11] + tx);

	const int NBaseline = 5000;
	double FilterRange = 40;
	double ori, nan, baseline;
	TH1D* Raw;
	TH1D* Raw_baseline;
	TH1* Freq;
	Raw = NULL;
	Freq = NULL;
	TH2D* Freq_Spec = new TH2D("Frequency", "Frequency", Nsize, 0, Nsize, 2000, 0, 2000);
	TH2D* Freq_Spec_bl = new TH2D("Frequency of baseline", "Frequency of baseline", Nsize, 0, Nsize, 2000, 0, 2000);
	TH1D* proj;
	TH1D* mean = new TH1D("Frequency Mean", "Frequency Mean", Nsize, 0, Nsize);
	TH1D* mean_bl = new TH1D("Baseline Frequency Mean", "Baseline Frequency Mean", Nsize, 0, Nsize);
	TH1D* Filter = new TH1D("Filter", "Filter", Nsize, 0, Nsize);
	TFile* file = new TFile(name + "_Filter" + ro, "recreate");
	ifstream in;
	in.open(name + tx);
	if (!in.good()) {
		cout << "Wrong" << endl;
		return false;
	}
	for (int i = 0; i < WavNum; i ++)
	{
		Raw = new TH1D("Waveform", "Waveform", Nsize, 0, Nsize);
		Raw_baseline = new TH1D("Noise", "Noise", Nsize, 0, Nsize);
		// Read signal file
		for (int j = 0; j < Nsize; j ++)
		{
			in >> ori >> nan;
			Raw->SetBinContent(j + 1, - ori * 1000);
		}
		// Read noise file
		if (i < NBaseline)
			for (int j = 0; j < Nsize; j ++)
			{
				in_bl >> ori >> nan;
				Raw_baseline->SetBinContent(j + 1, - ori * 1000);
			}
		// Get noise frequency domain spectrum
		if (i)
			delete TVirtualFFT::GetCurrentTransform();
		TVirtualFFT::SetTransform(0);
		Freq = Raw_baseline->FFT(Freq, "MAG");
		for (int j = 1; j <= Nsize; j ++)
			Freq_Spec_bl->Fill(Freq->GetBinCenter(j), Freq->GetBinContent(j));
		Freq->Delete();
		Raw_baseline->Delete();
		Raw_baseline = NULL;
		Freq = NULL;
		baseline = 0;
		// Get signal frequency domain spectrum
		for (int j = 0; j < 300; j ++)
			baseline += Raw->GetBinContent(j + 1);
		baseline /= 300;
		if (Raw->GetBinContent(Raw->GetMaximumBin()) - baseline < Threshold)
		{
			Raw->Delete();
			Raw = NULL;
			continue;
		}
		delete TVirtualFFT::GetCurrentTransform();
		TVirtualFFT::SetTransform(0);
		Freq = Raw->FFT(Freq, "MAG");
		for (int j = 1; j <= Nsize; j ++)	
			Freq_Spec->Fill(Freq->GetBinCenter(j), Freq->GetBinContent(j));
		Raw->Delete();
		Freq->Delete();
		Raw = NULL;
		Freq = NULL;
	}
	file->cd();
	for (int i = 1; i <= Nsize; i ++)
	{
		name = "ProjectionY_Baseline";
		proj = Freq_Spec_bl->ProjectionY(name, i, i ,"");
		proj->Smooth(3);
		mean_bl->SetBinContent(i, proj->GetMaximumBin());
		proj->Delete();
		proj = NULL;
	}
	for (int i = 1; i <= Nsize; i ++)
	{
		name = "ProjectionY";
		proj = Freq_Spec->ProjectionY(name, i, i, "");
		proj->Smooth(3);
		mean->SetBinContent(i, proj->GetMaximumBin());
		proj->Delete();
		proj = NULL;
	}
	double Sig, Noi = 0;
	// for (int i = 100; i < 200; i ++)
	// 	Noi += mean->GetBinContent(i);
	// Noi /= 100;
	out << "Signal\tNoise" << endl;
	for (int i = 1; i <= Nsize; i ++)
	{
		if (i < 701)
		{
			Sig = mean->GetBinContent(Nsize + 1 - i);
			Noi = mean_bl->GetBinContent(Nsize + 1 - i);
		}
		else
		{
			Sig = mean->GetBinContent(i);
			Noi = mean_bl->GetBinContent(i);
		}
		Filter->SetBinContent(i, pow(Sig, 2) < pow(Noi, 2) ? 0 : ((pow(Sig, 2) - pow(Noi, 2)) / pow(Sig, 2)));
		
		if (i > FilterRange && i < Nsize - FilterRange)
			Filter->SetBinContent(i, 0);
		out << Sig << "\t" << Noi << "\t" << endl;
	}
	TCanvas* can2=new TCanvas("can2","can2",800,600);
	mean->Smooth(5);
	mean_bl->Smooth(5);
	mean->Draw();
	mean_bl->Draw("same");

	Filter->Write();
	mean->Write();
	mean_bl->Write();
	Freq_Spec->Write();
	Freq_Spec_bl->Write();
	file->Close();
	in.close();
	return true;
}
