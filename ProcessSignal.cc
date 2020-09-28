#include "GetFilter.h"

using namespace std;

void format_h(TH1* h){
    h->SetLineWidth(3);
    // h->SetLineColor(linecolor);
    h->SetTitleOffset(0.5, "x");
    h->SetTitleOffset(0.8, "yz");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.85);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.65);
    }
double interpolate(double x1, double y1, double x2, double y2, double y)
{
	return (y - y1) * (x2 - x1) / (y2 - y1) + x1;
}
int main(int argc, char **argv)
{
	double Scale, Threshold, Charge_Thr;// rise time 1525: 0.55; 1305: 0.4
	int Nsize, WavNum, step, type, RangeA, RangeB;
	if (argc != 11 && argc != 12)
	{
		cout << "Input arguments: " << endl;
		cout << "[0]: './getParas'." << endl;
		cout << "[1]: Step (1 or 2)." << endl;
		cout << "[2]: Input type (1 or 2)." << endl;
		cout << "[3]: Input file name." << endl;
		cout << "[4]: Lower bound." << endl;
		cout << "[5]: Higher bound." << endl;
		cout << "[6]: Scale." << endl;
		cout << "[7]: The number of sampling points." << endl;
		cout << "[8]: The number of waveforms." << endl;
		cout << "[9]: Amplitude threshold." << endl;
		cout << "[10]: Charge threshold." << endl;
		cout << "[11]: Baseline input file name (just for step 1)." << endl;
		cout << "Wrong number (" << argc << ") of input arguments" << endl;
		return 0;
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

	if (step == 1 && argc == 12) {
		if (GetFilter(argc, argv))
			return 1;
		return 0;
	}
	else if (step != 2) {
		cout << "Wrong step number (" << step << ")." << endl;
		cout << "It must be 1 or 2." << endl;
		return 0;
	}

	if (RangeA < 0) {
		cout << "Error: lower bound (" << RangeA <<") is smaller than 0." << endl;
		return 0;
	}
	if (RangeB > Nsize) {
		cout << "Error: higher bound (" << RangeB << ") is bigger than the number of sample points (" << Nsize << ")." << endl;
		return 0;
	}
	if (RangeA > RangeB) {
		cout << "Error: lower bound (" << RangeA << ") is bigge than higher bound (" << RangeB << ")." << endl;
		return 0;
	}

	char nn[16];
	sprintf(nn, "_%lf", Threshold);
	int counter = 0;
	double ori, nan;
	TFile* file = new TFile(name + "_Filter.root", "read");
	ifstream in(name + ".txt");
	TH1D* Filter;
	TH1D* Raw;
	TH1D* tmp;
	TH1* Freq;
	TH1* InFFT_RE;
	TH2D* InFFT_Coll = new TH2D("Inverse_FFT_Collection", "Inverse_FFT_Collection;/ mV;#times 0.2ns", Nsize, 0, Nsize, 300, -10, 20);
	TH1D* Raw_samp = new TH1D("Raw_Sample", "Raw_Sample", Nsize, 0, Nsize);
	TH1D* InFFT_samp = new TH1D("Inverse_FFT_Sample", "Inverse_FFT_Sample", Nsize, 0, Nsize);
	TH1D* Amplitude = new TH1D("Amplitude Spectrum", "Amplitude Spectrum;/ mV;Counts", 200, 0, 20);
	TH1D* Charge = new TH1D("Charge Spectrum", "Charge Spectrum;/ pC;Counts", 500, -1, 4);// 80ns
	TH1D* Charge_100ns = new TH1D("Charge Spectrum with 100ns integral interval", "Charge Spectrum with 100ns integral interval;/ pC;Counts", 500, -1, 4);
	TH1D* Charge_120ns = new TH1D("Charge Spectrum with 120ns integral interval", "Charge Spectrum with 120ns integral interval;/ pC;Counts", 500, -1, 4);
	TH1D* Charge_140ns = new TH1D("Charge Spectrum with 140ns integral interval", "Charge Spectrum with 140ns integral interval;/ pC;Counts", 500, -1, 4);
	TH1D* Rise = new TH1D("Rise Time", "Rise Time;#times 0.2ns;Counts", 60, 0, 30);
	TH1D* Fall = new TH1D("Fall Time", "Fall Time;#times 0.2ns;Counts", 60, 0, 30);
	TH1D* FWHM = new TH1D("FWHM", "FWHM;#times 0.2ns;Counts", 60, 0, 30);
	TH1D* AveWav = new TH1D("Average waveform", "Average waveform", 1000, 0, 1000);
	TH2D* risevamp = new TH2D("Rise time vs. Amplitude", "Rise time vs. Amplitude;/ mV;/ ns", 200, 0, 20, 60, 0, 30);
	TH2D* fallvamp = new TH2D("Fall time vs. Amplitude", "Fall time vs. Amplitude;/ mV;/ ns", 200, 0, 20, 60, 0, 30);
	TH2D* risevcharge = new TH2D("Rise time vs. Charge", "Rise time vs. Charge;/ pC;/ ns", 500, -1, 4, 60, 0, 30);
	TH2D* fallvcharge = new TH2D("Fall time vs. Charge", "Fall time vs. Charge;/ pC;/ ns", 500, -1, 4, 60, 0, 30);
	TH2D* ampvcharge = new TH2D("Amplitude vs. Charge", "Amplitude vs. Charge;/ pC;/ mV", 500, -1, 4, 200, 0, 20);
	InFFT_RE = NULL;
	Filter = NULL;
	Raw = NULL;
	Freq = NULL;
	bool done = false;
	// ifstream in("./" + name + ".txt");
	ofstream out(name + "_wf.txt");
	ofstream ampcharge("./ampcharge.txt", std::ios::app);
	file->GetObject("Filter", Filter);
	if (!in.good()) {
		cout << "Open input txt error!" << endl;
		return 0;
	}
	if (!file) {
		cout << "Read in filter error!" << endl;
		return 0;
	}
	cout << "Processing..." << endl;
	for (int i = 0; i < WavNum; i ++)
	{
		Raw = new TH1D("Waveform", "Waveform", Nsize, 0, Nsize);
		tmp = new TH1D("tmp", "tmp", Nsize, 0, Nsize);
		for (int j = 1; j <= Nsize; j ++)
		{
			if (type == 2)
				in >> ori >> nan;
			else
				in >> ori;
			Raw->SetBinContent(j, - ori * 1000);
		}
		if (i)
			delete TVirtualFFT::GetCurrentTransform();
		TVirtualFFT::SetTransform(0);
		Freq = Raw->FFT(Freq, "MAG");
		std::vector<double> re_full_vec(Nsize);
		std::vector<double> im_full_vec(Nsize);
		double *re_full = &re_full_vec.front();
		double *im_full = &im_full_vec.front();
		for (int j = 0; j < Nsize; j ++)
		{
			re_full[j] = 0;
			im_full[j] = 0;
		}
		TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
		if (!fft)
			cout << "error" << endl;
		fft->GetPointsComplex(re_full, im_full);
		for (int j = 0; j < Nsize; j ++)
		{
			re_full[j] *= Filter->GetBinContent(j + 1);
			im_full[j] *= Filter->GetBinContent(j + 1);
		}
		int N = Nsize;
		TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2R M K");
		fft_back->SetPointsComplex(re_full, im_full);
		fft_back->Transform();
		InFFT_RE = TH1::TransformHisto(fft_back, InFFT_RE, "RE");
		InFFT_RE->Scale(1. / Nsize);
		for (int j = 1; j <= Nsize; j ++)
			InFFT_Coll->Fill(j, InFFT_RE->GetBinContent(j));	

		//Get signals' parameters. Amplitude, Rise time, Fall time, FWHM, Gain.
		double baseline = 0;
		double peakpos = InFFT_RE->GetMaximumBin();
		if (peakpos > RangeB || peakpos < RangeA)
		{
			Raw->Delete();
			Freq->Delete();
			InFFT_RE->Delete();
			tmp->Delete();
			Raw = NULL;
			Freq = NULL;
			InFFT_RE = NULL;
			continue;
		}
		double amp, charge = 0;
		for (int j = 1; j < 300; j ++)
			baseline += InFFT_RE->GetBinContent(j);
		baseline /= 300;
		// 80ns integral interval
		for (int j = 1; j <= Nsize; j ++)
		{
			double p = InFFT_RE->GetBinContent(j) - baseline;
			tmp->SetBinContent(j, p);
			if (j >= 500 && j < 900)
				charge += p;
		}
		// charge *= Scale / 50;
		amp = tmp->GetBinContent(peakpos);
		Charge->Fill(charge * Scale / 50);
		// 100ns integral interval
		for (int j = 450; j < 500; j ++)
		{
			double p = InFFT_RE->GetBinContent(j) - baseline;
			charge += p;
		}
		for (int j = 900; j < 950; j ++)
		{
			double p = InFFT_RE->GetBinContent(j) - baseline;
			charge += p;
		}
		Charge_100ns->Fill(charge * Scale / 50);
		for (int j = 400; j < 450; j ++)
		{
			double p = InFFT_RE->GetBinContent(j) - baseline;
			charge += p;
		}
		for (int j = 950; j < 1000; j ++)
		{
			double p = InFFT_RE->GetBinContent(j) - baseline;
			charge += p;
		}
		Charge_120ns->Fill(charge * Scale / 50);
		for (int j = 350; j < 400; j ++)
		{
			double p = InFFT_RE->GetBinContent(j) - baseline;
			charge += p;
		}
		for (int j = 1000; j < 1050; j ++)
		{
			double p = InFFT_RE->GetBinContent(j) - baseline;
			charge += p;
		}
		Charge_140ns->Fill(charge * Scale / 50);
		Amplitude->Fill(amp);
		double rise_s = 0, rise_e = 0, mid_s = 0;
		for (int j = peakpos - 200; j <= peakpos; j ++)
		{
			double p = tmp->GetBinContent(j);
			if (p < 0.1 * amp)
				rise_s = j;
			if (p < 0.5 * amp)
				mid_s = j;
			if (p < 0.9 * amp)
				rise_e = j;
		}
		double fall_s = 0, fall_e = 0, mid_e = 0;
		for (int j = peakpos + 200; j >= peakpos; j --)
		{
			double p = tmp->GetBinContent(j);
			if (p < 0.9 * amp)
				fall_s = j;
			if (p < 0.5 * amp)
				mid_e = j;
			if (p < 0.1 * amp)	
				fall_e = j;
		}
		// interpolation
		rise_s = interpolate(rise_s, tmp->GetBinContent(rise_s), rise_s + 1, tmp->GetBinContent(rise_s + 1), 0.1 * amp);
		rise_e = interpolate(rise_e, tmp->GetBinContent(rise_e), rise_e + 1, tmp->GetBinContent(rise_e + 1), 0.9 * amp);
		fall_s = interpolate(fall_s, tmp->GetBinContent(fall_s), fall_s + 1, tmp->GetBinContent(fall_s + 1), 0.9 * amp);
		fall_e = interpolate(fall_e, tmp->GetBinContent(fall_e), fall_e + 1, tmp->GetBinContent(fall_e + 1), 0.1 * amp);
		mid_s = interpolate(mid_s, tmp->GetBinContent(mid_s), mid_s + 1, tmp->GetBinContent(mid_s + 1), 0.5 * amp);
		mid_e = interpolate(mid_e, tmp->GetBinContent(mid_e), mid_e + 1, tmp->GetBinContent(mid_e + 1), 0.5 * amp);

		double rise_t = (rise_e - rise_s) * Scale, fall_t = (fall_e - fall_s) * Scale, fwhm = (mid_e - mid_s) * Scale;

		if (amp > 0.4 && charge > 0.3 && !done)
		{
			for (int j = 1; j < Nsize; j ++)
			{
				Raw_samp->SetBinContent(j, Raw->GetBinContent(j));
				InFFT_samp->SetBinContent(j, InFFT_RE->GetBinContent(j));
			}
			done = true;
		}

		// 10m_1305: (2 * amp - 0.98 * charge * Scale / 50 - 1) > 0
		// 1m_1305: amp > 0.4 && charge > 0.3
		// 5m_1305: amp > 0.4 && charge > 0.3
		// 1m_1525: amp > 0.7 && charge > 0.5
		// 5m_1525: amp > 0.7 && charge > 0.5
		// 10m_1525: (1.6 * amp - 0.6 * charge * Scale / 50 - 1.14) > 0
		// if (amp > 0.4 && charge > 0.3)
		{
			Fall->Fill(fall_t);
			Rise->Fill(rise_t);
			FWHM->Fill(fwhm);
		}
		risevamp->Fill(amp, rise_t);
		fallvamp->Fill(amp, fall_t);
		risevcharge->Fill(charge * Scale / 50, rise_t);
		fallvcharge->Fill(charge * Scale / 50, fall_t);
		ampvcharge->Fill(charge * Scale / 50, amp);
		// 10m_1305: (2 * amp - 0.98 * charge * Scale / 50 - 1) > 0
		// 1m_1305: amp > 0.4 && charge > 0.3
		// 5m_1305: amp > 0.4 && charge > 0.3
		// 1m_1525: amp > 0.7 && charge > 0.5
		// 5m_1525: amp > 0.7 && charge > 0.5
		// 10m_1525: (1.6 * amp - 0.6 * charge * Scale / 50 - 1.14) > 0
		if (amp > Threshold && charge > Charge_Thr)
		{
			int n = rise_s - 400;
			for (int j = 1; j <= 1000; j ++)
			{
				AveWav->SetBinContent(j, tmp->GetBinContent(n) + AveWav->GetBinContent(j));	
				n ++;
			}
			counter ++;
		}
		tmp->Delete();
		//
		Raw->Delete();
		Freq->Delete();
		InFFT_RE->Delete();
		Raw = NULL;
		Freq = NULL;
		InFFT_RE = NULL;
	}
	double charge = 0;
	for (int j = 1; j <= 1000; j ++)
	{
		AveWav->SetBinContent(j, AveWav->GetBinContent(j) / counter);
		charge += AveWav->GetBinContent(j);
		out << AveWav->GetBinContent(j) << endl;
	}
	charge *= Scale / 50;
	out << charge << endl;
	double peakpos = AveWav->GetMaximumBin();
	double amp = AveWav->GetBinContent(peakpos);
	double rise_s = 0, rise_e = 0, mid_s = 0;
	for (int j = peakpos - 200; j <= peakpos; j ++)
	{
		double p = AveWav->GetBinContent(j);
		if (p < 0.1 * amp)
			rise_s = j;
		if (p < 0.5 * amp)
			mid_s = j;
		if (p < 0.9 * amp)
			rise_e = j;
	}
	double fall_s = 0, fall_e = 0, mid_e = 0;
	for (int j = peakpos + 200; j >= peakpos; j --)
	{
		double p = AveWav->GetBinContent(j);
		if (p < 0.9 * amp)
			fall_s = j;
		if (p < 0.5 * amp)
			mid_e = j;
		if (p < 0.1 * amp)	
			fall_e = j;
	}
	double rise_t = (rise_e - rise_s) * Scale, fall_t = (fall_e - fall_s) * Scale, fwhm = (mid_e - mid_s) * Scale;
	out << amp << endl;
	out << rise_t << endl;
	out << fall_t << endl;
	out << fwhm << endl;
	// ampcharge << charge << "\t" << amp << "\t" << rise_t << "\t" << fall_t << "\t" << fwhm << "\t" << Threshold << endl;
	// cout << amp << "\t" << charge << "\t" << rise_t << "\t" << fall_t << "\t" << fwhm << endl;
	for (int j = 1; j <= 500; j ++)
		out << Charge->GetBinContent(j) << endl;
	/*
	for (auto hist : {Amplitude, Charge, Rise, Fall, FWHM})
		format_h(hist);
	format_h(risevamp);
	TCanvas* fig1 = new TCanvas("Rise Time vs. Amplitude", "Rise Time vs. Amplitude", 1920, 1080);// "Rise vs. Amplitude", "Rise vs. Amplitude", 800, 600
	risevamp->Draw("colz");
	fig1->Print("Rise_vs_amp.pdf");
	// delete(fig1);
	TCanvas* fig2 = new TCanvas("Amplitude", "Amplitude", 1920, 1080);// "Amplitude", "Amplitude", 800, 600
	fig2->SetLogy();
	Amplitude->Draw();
	fig2->Print("Amplitude.pdf");
	// delete(fig2);
	TCanvas* fig3 = new TCanvas();
	fig3->Divide(1, 2, 0, 0);
	fig3->cd(1);
	TArrow arrow1(1, 0, 1, 10000, 0, "|>");
	TLatex text(1.2, 1000, "Amplitude Threshold");
	arrow1.SetLineWidth(2);
	arrow1.SetLineColor(kRed);
	fig3->GetPad(1)->SetLogy();
	Amplitude->SetStats(kFALSE);
	Amplitude->Draw();
	arrow1.DrawClone();
	text.DrawClone();

	fig3->cd(2);
	TArrow arrow2(1, 0, 1, 30, 0, "|>");
	arrow2.SetLineWidth(2);
	arrow2.SetLineColor(kRed);
	risevamp->SetStats(kFALSE);
	risevamp->Draw("colz");
	arrow2.DrawClone();
	fig3->Print("How_I_choose_amp_threshold.pdf");
	// delete(fig3);
	TCanvas* fig4 = new TCanvas("Fall Time vs. Amplitude", "Fall Time vs. Amplitude", 1920, 1080);// "Rise vs. Amplitude", "Rise vs. Amplitude", 800, 600
	fallvamp->Draw("colz");
	fig4->Print("Fall_vs_amp.pdf");
	// delete(fig4)
	TCanvas* fig5 = new TCanvas("Charge vs. Amplitude", "Charge vs. Amplitude", 1920, 1080);// "Rise vs. Amplitude", "Rise vs. Amplitude", 800, 600
	ampvcharge->Draw("colz");
	fig4->Print("Q_vs_amp.pdf");
	// delete(fig5)
	*/
	name += nn;
	TFile *file_out = new TFile(name + "_out.root", "recreate");
	file_out->cd();
	Amplitude->Write();
	Charge->Write();
	ampvcharge->Write();
	Charge_100ns->Write();
	Charge_120ns->Write();
	Charge_140ns->Write();
	Rise->Write();
	risevamp->Write();
	risevcharge->Write();
	Fall->Write();
	fallvamp->Write();
	fallvcharge->Write();
	FWHM->Write();
	InFFT_Coll->Write();
	Raw_samp->Write();
	InFFT_samp->Write();
	AveWav->Write();
	file_out->Close();
	in.close();
	ampcharge.close();
	file->Close();
	return 1;
}
