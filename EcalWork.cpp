#include "MpdDataConverter.h"
#include <cstddef>
#include <ostream>
#include <iostream>
#include <ctime>
#include <string>
#include <cstdlib>
#include <vector>

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMathBase.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"

#include <filesystem>
#include "TSystem.h"
#include "TF1.h"
#include "TGraphErrors.h"

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 0, 0)
R__LOAD_LIBRARY(./build/libMPDDataConverter.dylib)
#endif

template <typename T>
void SetWaveformHistogram(
    const std::vector<T> &data,
    TH1F *hist,
    const std::string &name = "waveform",
    const std::string &title = "",
    const std::string &xTitle = "Sample Number",
    const std::string &yTitle = "ADC")
{
    if (data.empty())
        return;

    int nSamples = data.size();

    // TH1F* hist = new TH1F(name.c_str(), title.c_str(), nSamples, 0, nSamples);

    hist->SetName(name.c_str());
    hist->SetTitle(title.c_str());

    for (int i = 0; i < nSamples; i++)
    {
        hist->SetBinContent(i + 1, static_cast<Float_t>(data[i]));
    }

    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());

    hist->SetLineColor(kBlue);
    hist->SetLineWidth(2);
    hist->SetDrawOption("L");
    hist->SetBit(TH1::kNoStats);
}

void DrawAllWaveformsOnOneCanvas(std::vector<ChannelData> _ChData,
                                 std::vector<TH1F *> &histograms,
                                 TH2D *h2_histo,
                                 std::vector<TLine *> &_userGrid,
                                 TDirectory *eventDir)
{

    if (_ChData.empty())
        return;

    int eventNumber = _ChData[0].eventNum;

    std::vector<Float_t> axis_y_min;
    std::vector<Float_t> axis_y_max;

    for (size_t i = 0; i < _ChData.size(); i++)
    {
        int chNum = _ChData[i].channelNums;
        int chPhi = _ChData[i].channel_Phi;
        int chZ = _ChData[i].channel_Z;
        int chInt = _ChData[i].integral;
        std::string h_name = Form("h_ev_%i_ch_%i_phi_%i_z_%i", eventNumber, chNum, chPhi, chZ);
        std::string h_title = Form("Basket 38, Event %i, channel %i, Phi %i, Z %i", eventNumber, chNum, chPhi, chZ);

        SetWaveformHistogram(_ChData[i].adcValues, histograms[chNum - 1], h_name, h_title);

        h2_histo->SetBinContent(h2_histo->FindBin(chZ, chPhi), chInt);

        axis_y_max.push_back(histograms[chNum - 1]->GetMaximum());
        axis_y_min.push_back(histograms[chNum - 1]->GetMinimum());
    }

    Float_t AxisYmax = *std::max_element(axis_y_max.begin(), axis_y_max.end());
    Float_t AxisYmin = *std::min_element(axis_y_min.begin(), axis_y_min.end());

    TCanvas *canvas = new TCanvas(Form("wf_event_%d", eventNumber), Form("Event %d - All Waveforms Overlay", eventNumber), 3 * 720, 2 * 720);

    TPad *topPad = new TPad("topPad", "Top Pad", 0, 0.35, 1, 1);
    TPad *bottomPad = new TPad("bottomPad", "Bottom Pad", 0, 0, 1, 0.35);

    topPad->SetFillColor(0);
    topPad->SetFrameFillColor(0);

    bottomPad->SetFillColor(0);
    bottomPad->SetFrameFillColor(0);
    bottomPad->SetTopMargin(0);
    bottomPad->SetBottomMargin(0.2);

    canvas->cd();
    topPad->Draw();
    bottomPad->Draw();

    std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan, kYellow,
                               kTeal, kPink, kViolet, kSpring, kAzure, kOrange};
    std::vector<int> colorTone = {-2, -1, 0, 1, 2, 3};
    std::vector<int> LineColor;

    LineColor.reserve(colors.size() * colorTone.size());
    for (int i = 0; i < colorTone.size(); i++)
    {
        for (int j = 0; j < colors.size(); j++)
        {
            LineColor.push_back(colors[j] + colorTone[i]);
        }
    }

    TH2F *AxisHisto = new TH2F(Form("h_axis_%i", eventNumber), "",
                               2, 0, 60, 2, AxisYmin - 0.1 * TMath::Abs(AxisYmin - AxisYmax), AxisYmax + 0.1 * TMath::Abs(AxisYmin - AxisYmax));
    AxisHisto->SetLineColor(colors[0 % colors.size()]);
    AxisHisto->SetBit(TH1::kNoStats);
    AxisHisto->SetLineWidth(2);
    topPad->cd();
    AxisHisto->Draw();

    AxisHisto->GetXaxis()->SetTitle("Sample Number");
    AxisHisto->GetYaxis()->SetTitle("ADC Value");
    AxisHisto->SetTitle(Form("Event %d - All %d Channels", eventNumber, (int)_ChData.size()));

    TLegend *legend = new TLegend(0.91, 0.1, 1.0, 0.9);
    legend->SetBorderSize(0);
    legend->SetNColumns(2);
    legend->SetFillStyle(1001);

    for (size_t i = 0; i < _ChData.size(); i++)
    {
        int iCh = _ChData[i].channelNums - 1;
        if (histograms[iCh]->GetEntries() < 1)
            continue;
        // histograms[iCh]->SetLineColor(colors[i % colors.size()]);
        histograms[iCh]->SetLineColor(LineColor[(iCh + 1) % 64]);
        histograms[iCh]->SetLineWidth(2);
        histograms[iCh]->Draw("L same");
        std::string legendEntry = Form("%i", _ChData[i].channelNums);
        legend->AddEntry(histograms[iCh], legendEntry.c_str(), "l");
    }

    legend->Draw("same");

    bottomPad->cd();
    // h2_histo->GetXaxis()->SetTitleOffset(AxisHisto->GetXaxis()->GetTitleOffset());
    h2_histo->GetXaxis()->SetTitleSize(2 * AxisHisto->GetXaxis()->GetTitleSize());
    h2_histo->GetXaxis()->SetLabelSize(2 * AxisHisto->GetXaxis()->GetLabelSize());
    h2_histo->GetYaxis()->SetTitleSize(2 * AxisHisto->GetYaxis()->GetTitleSize());
    h2_histo->GetYaxis()->SetLabelSize(2 * AxisHisto->GetYaxis()->GetLabelSize());
    h2_histo->GetYaxis()->SetTitleOffset(0.4);
    h2_histo->GetYaxis()->SetTickSize(0.005);
    // h2_histo->GetXaxis()->SetTitleSize(0.5);
    h2_histo->Draw("COLZ");

    for (int i = 0; i < _userGrid.size(); i++)
    {
        _userGrid[i]->SetLineColor(kBlack);
        _userGrid[i]->SetLineStyle(2);
        _userGrid[i]->SetLineWidth(1);
        _userGrid[i]->Draw("same");
    }

    bottomPad->Modified();
    bottomPad->Update();

    eventDir->cd();

    // canvas->Write();
    //  canvas->SaveAs(Form("/home/aleksandr/ecal_work/pict/%s.png",canvas->GetName()));
    canvas->SaveAs(Form("pict_%s.png", canvas->GetName()));

    for (size_t i = 0; i < _ChData.size(); i++)
    {
        int chNum = _ChData[i].channelNums;
        histograms[chNum - 1]->Reset();
    }

    h2_histo->Reset();

    delete canvas;
    delete legend;
    delete AxisHisto;

    canvas = nullptr;
    legend = nullptr;
    AxisHisto = nullptr;
}

Int_t SetPedestal(std::vector<Int_t> _adc_value, Int_t _Skip, Int_t _count)
{
    Int_t SumAdc = 0;
    Int_t nCount = 0;
    Int_t result = 0;
    if (_adc_value.size() == 0)
        return 0;

    for (int i = _Skip; i < (_Skip + _count); i++)
    {
        SumAdc += _adc_value[i];
        nCount++;
    }
    result = (Int_t)(SumAdc / nCount);
    return result;
}

void PedestalSubtraction(std::vector<Int_t> &_ch_Data, Int_t _pedestal)
{
    for (int i = 0; i < _ch_Data.size(); i++)
    {
        _ch_Data[i] = _ch_Data[i] - _pedestal;
    }
}

Int_t GetPedestalAmpl(std::vector<Int_t> _adc_value, Int_t _Skip, Int_t _count)
{
    Int_t SumAdc = 0;
    Int_t nCount = 0;
    Int_t result = 0;
    if (_adc_value.size() == 0)
        return 0;

    auto [minIt, maxIt] = std::minmax_element(_adc_value.begin(), _adc_value.begin() + _Skip + _count);

    return TMath::Abs(*maxIt - *minIt);
}

void SetChannelIntegralAmp(ChannelData &_chData, TH1F *hist, int binStart, int binEnd, int nBinLeft, int nBinRight)
{
    if (!hist)
        return;
    if (hist->GetEntries() == 0)
        return;
    // Пик
    Float_t maxValue = 0;
    Int_t peakBin = 0;

    for (int i = binStart; i <= binEnd; i++)
    {
        Float_t currentValue = TMath::Abs(hist->GetBinContent(i));
        if (currentValue > maxValue)
        {
            maxValue = currentValue;
            peakBin = i;
        }
    }

    Int_t int_left = (peakBin - nBinLeft) <= 0 ? 20 : peakBin - nBinLeft;
    Int_t int_right = (peakBin + nBinRight) >= 60 ? 46 : peakBin + nBinRight;

    double peakX = hist->GetBinCenter(peakBin);
    double peakY = hist->GetBinContent(peakBin);

    // Интеграл
    _chData.integral = hist->Integral(int_left, int_right);
    _chData.amplitude = peakY;
    //_chData.Print();
}

void ResetTH1Fvector(std::vector<TH1F *> &hists, std::vector<int> iNumClear)
{
    for (int i = 0; i < iNumClear.size(); i++)
    {
        hists[iNumClear[i]]->Reset();
    }
}

void FitTargetChannelHistograms(TFile *outFile,
                                const std::string &inputDirName = "target_channels",
                                const std::string &outputDirName = "fitted_target_channels",
                                int minEntriesToFit = 50,
                                int rebinFactor = 2)
{
    if (!outFile || outFile->IsZombie())
    {
        std::cout << "FitTargetChannelHistograms: output file is null or zombie" << std::endl;
        return;
    }

    TDirectory *dirIn = (TDirectory *)outFile->Get(inputDirName.c_str());
    if (!dirIn)
    {
        std::cout << "FitTargetChannelHistograms: directory " << inputDirName << " not found" << std::endl;
        return;
    }

    outFile->cd();
    TDirectory *dirOut = (TDirectory *)outFile->Get(outputDirName.c_str());
    if (!dirOut)
    {
        dirOut = outFile->mkdir(outputDirName.c_str());
    }
    dirOut->cd();

    const int NCHANNELS = 768;

    std::vector<double> v_ch;
    std::vector<double> v_mpv;
    std::vector<double> v_mpv_err;
    std::vector<double> v_width;
    std::vector<double> v_width_err;
    std::vector<double> v_entries;

    int nFitted = 0;

    for (int ch = 1; ch <= NCHANNELS; ++ch)
    {
        TH1F *hIn = (TH1F *)dirIn->Get(Form("h_int_ch_%d", ch));
        if (!hIn)
            continue;

        if (hIn->GetEntries() < minEntriesToFit)
            continue;

        TH1F *hFit = (TH1F *)hIn->Clone(Form("h_fit_ch_%d", ch));
        hFit->SetDirectory(nullptr);

        if (rebinFactor > 1)
            hFit->Rebin(rebinFactor);

        int nBins = hFit->GetNbinsX();
        if (nBins < 10)
        {
            delete hFit;
            continue;
        }

        // Find the global maximum bin after rebinning
        int maxBin = hFit->GetMaximumBin();
        double peakX = hFit->GetBinCenter(maxBin);
        double peakY = hFit->GetBinContent(maxBin);

        if (peakY <= 0)
        {
            delete hFit;
            continue;
        }

        // Find left/right bins above 35% of peak to estimate the signal hump width
        double frac = 0.35 * peakY;
        int leftBin = maxBin;
        int rightBin = maxBin;

        while (leftBin > 1 && hFit->GetBinContent(leftBin) > frac)
            --leftBin;
        while (rightBin < nBins && hFit->GetBinContent(rightBin) > frac)
            ++rightBin;

        double fitMin = hFit->GetBinCenter(leftBin);
        double fitMax = hFit->GetBinCenter(rightBin);

        // Expand the window a bit to accommodate the Landau tail
        double widthGuess = peakX - fitMin;
        if (widthGuess <= 0)
            widthGuess = hFit->GetRMS() / 4.0;
        if (widthGuess <= 0)
            widthGuess = 500.0;

        fitMin = std::max(0.0, peakX - 1.5 * widthGuess);
        fitMax = peakX + 4.0 * widthGuess;

        // Avoid fitting the low-integral pedestal if the peak is well away from zero
        if (peakX > 1500.0 && fitMin < 0.5 * peakX)
            fitMin = 0.5 * peakX;

        if (fitMax <= fitMin)
        {
            delete hFit;
            continue;
        }

        TF1 *fLan = new TF1(Form("f_landau_ch_%d", ch), "landau", fitMin, fitMax);
        fLan->SetParameters(peakY, peakX, widthGuess);
        fLan->SetParLimits(1, 4000.0, 9000.0); // MPV between 3k and 9k
        fLan->SetLineColor(kRed);
        fLan->SetLineWidth(2);

        int fitStatus = hFit->Fit(fLan, "QRS");
        if (fitStatus != 0)
        {
            delete hFit;
            delete fLan;
            continue;
        }

        double mpv = fLan->GetParameter(1);
        double mpvErr = fLan->GetParError(1);
        double width = fLan->GetParameter(2);
        double widthErr = fLan->GetParError(2);
        double entries = hFit->GetEntries();

        // Quality cuts on fit
        bool badFit = false;
        if (mpv < 4000.0 || mpv > 9000.0)
            badFit = true;
        if (width <= 0.0 || width > 6000.0)
            badFit = true;
        if (mpv <= 0.0 || mpvErr / mpv > 0.5)
            badFit = true;

        if (badFit)
        {
            delete hFit;
            delete fLan;
            continue;
        }

        dirOut->cd();
        hFit->Write();
        fLan->Write();

        v_ch.push_back(ch);
        v_mpv.push_back(mpv);
        v_mpv_err.push_back(mpvErr);
        v_width.push_back(width);
        v_width_err.push_back(widthErr);
        v_entries.push_back(entries);

        nFitted++;
    }

    if (nFitted > 0)
    {
        TGraphErrors *grMPV = new TGraphErrors(
            nFitted,
            v_ch.data(), v_mpv.data(),
            nullptr, v_mpv_err.data());
        grMPV->SetName("gr_mpv_vs_channel");
        grMPV->SetTitle("Landau MPV vs channel;Channel;Landau MPV");
        grMPV->Write();

        TGraphErrors *grWidth = new TGraphErrors(
            nFitted,
            v_ch.data(), v_width.data(),
            nullptr, v_width_err.data());
        grWidth->SetName("gr_width_vs_channel");
        grWidth->SetTitle("Landau width vs channel;Channel;Landau width");
        grWidth->Write();

        TGraph *grEntries = new TGraph(
            nFitted,
            v_ch.data(), v_entries.data());
        grEntries->SetName("gr_entries_vs_channel");
        grEntries->SetTitle("Entries vs channel;Channel;Entries");
        grEntries->Write();
    }

    outFile->cd();

    std::cout << "FitTargetChannelHistograms: fitted " << nFitted
              << " channel histograms with Landau into directory "
              << outputDirName << std::endl;
}

// Define what "neighbor" means in (Phi, Z)
inline bool AreNeighbors(const ChannelData &a, const ChannelData &b,
                         int maxDeltaPhi = 1, int maxDeltaZ = 1)
{
    int dPhi = std::abs(a.channel_Phi - b.channel_Phi);
    int dZ = std::abs(a.channel_Z - b.channel_Z);
    return (dPhi <= maxDeltaPhi && dZ <= maxDeltaZ && !(dPhi == 0 && dZ == 0));
}

// Compute contamination around a given "track cluster"
// Here all channels are taken in ChannelDataInEvent as "selected",
// and everything in neighboring cells that is not selected is treated as contamination.
double ComputeNeighborContamination(const std::vector<ChannelData> &eventHits,
                                    double &sumSelected,
                                    double &sumNeighbors,
                                    int maxDeltaPhi = 1,
                                    int maxDeltaZ = 1)
{
    sumSelected = 0.0;
    sumNeighbors = 0.0;

    if (eventHits.empty())
        return 0.0;

    // First, sum integral in all selected hits ("target + selecting channels").
    // For now, treat the entire ChannelDataInEvent as the selected cluster.
    for (const auto &ch : eventHits)
        sumSelected += static_cast<double>(ch.integral);

    // Now, estimate contamination from neighbors that *should not* belong to the same track.
    // If later identified a core subset as the track, it can be treated the others differently.
    const int nHits = static_cast<int>(eventHits.size());

    for (int i = 0; i < nHits; ++i)
    {
        const ChannelData &ci = eventHits[i];

        // For each hit, look at all other hits that are "neighbor" in Phi/Z
        for (int j = 0; j < nHits; ++j)
        {
            if (i == j)
                continue;

            const ChannelData &cj = eventHits[j];

            if (AreNeighbors(ci, cj, maxDeltaPhi, maxDeltaZ))
            {
                // For a simple implementation, treat neighbor integrals as contamination
                // if they are not part of a minimal 1D track (e.g. same Phi).
                // Here: any neighbor with different Phi is "suspicious".
                if (cj.channel_Phi != ci.channel_Phi)
                {
                    sumNeighbors += static_cast<double>(cj.integral);
                }
            }
        }
    }

    // Avoid double counting contamination: it just summed many neighbors multiple times.
    // A simple fix is to normalize by number of contributing pairs, or better:
    // recompute sumNeighbors via a unique set. Here is a safer version:

    std::set<int> neighborIndices;
    for (int i = 0; i < nHits; ++i)
    {
        const auto &ci = eventHits[i];
        for (int j = 0; j < nHits; ++j)
        {
            if (i == j)
                continue;
            const auto &cj = eventHits[j];
            if (AreNeighbors(ci, cj, maxDeltaPhi, maxDeltaZ) &&
                cj.channel_Phi != ci.channel_Phi)
            {
                neighborIndices.insert(j);
            }
        }
    }

    sumNeighbors = 0.0;
    for (int idx : neighborIndices)
        sumNeighbors += static_cast<double>(eventHits[idx].integral);

    const double denom = sumSelected + sumNeighbors;
    if (denom <= 0.0)
        return 0.0;

    return sumNeighbors / denom;
}

void EcalWork(std::string inputData = "../run_rc-hs1_088.data",
              std::string outputData = "out2.root",
              int targetEvent = -1)
{

    TStopwatch timer1;
    timer1.Start();

    // Число первых бинов для пропуска и число последующих бинов для подсчета подложки
    Int_t iBinPedestalSkip = 2;
    Int_t iBinPedestalCount = 10;
    // Для поиска пика
    Int_t iBinStart = 20;
    Int_t iBinStop = 50;
    // Интеграл вокруг пика
    Int_t iBinLeft = 4;
    Int_t iBinRight = 16;

    const int NCHANNELS = 768;
    std::vector<TH1F *> h_int_per_channel;
    h_int_per_channel.reserve(NCHANNELS);
    for (int i = 0; i < NCHANNELS; ++i)
    {
        h_int_per_channel.push_back(
            new TH1F(Form("h_int_ch_%d", i + 1),
                     Form("Integral channel %d;Integral;Entries", i + 1),
                     100, 0, 15000));
    }

    InitChannelMap("basket_channel_map_phiZ.csv");

    MpdDataConverter converter;

    converter.OpenInputFile(inputData);
    converter.OpenOutRootFile(outputData);

    // Base directory = directory of the output ROOT file
    std::filesystem::path outRootPath(outputData);
    std::filesystem::path outDir = outRootPath.parent_path();
    if (outDir.empty())
    {
        outDir = ".";
    }

    // Create pict subfolder next to the ROOT file
    std::filesystem::path pictDir = outDir / "pict";
    std::filesystem::create_directories(pictDir);

    // Change process working directory so relative paths point there
    gSystem->ChangeDirectory(pictDir.string().c_str());

    std::cout << "Pictures will be saved under: " << pictDir << std::endl;

    // converter.WriteChannelSamplesVector();
    converter.SetTypeWaveForms(TypeWaveForms::invert);
    converter.InitTree();

    Int_t nBinsSamples = 60;
    Float_t SamplesMin = 0 - 0.5;
    Float_t SamplesMax = 60 - 0.5;

    TH1F *h1_wf = new TH1F("hist", "", nBinsSamples, SamplesMin, SamplesMax);

    Int_t MaxNChannels = 64 * 12;
    std::vector<TH1F *> eventHistograms;
    eventHistograms.reserve(MaxNChannels);
    for (int i = 0; i < MaxNChannels; i++)
    {
        eventHistograms.push_back(new TH1F(Form("histo_%i", i), "", nBinsSamples, SamplesMin, SamplesMax));
    }
    TH2D *h2_int_ch = new TH2D(Form("ch_integral"), ";channel Num;Integral", MaxNChannels, 0.5, MaxNChannels + 0.5, 700, 0, 700000);

    std::vector<ChannelData> ChannelDataInEvent;
    ChannelDataInEvent.reserve(MaxNChannels);

    TH2D *h2_integral_z_phi = new TH2D("h2", ";Z;Phi", 64, 0.5, 64.5, 12, 0.5, 12.5);
    h2_integral_z_phi->GetXaxis()->SetNdivisions(13, 5, kTRUE);
    h2_integral_z_phi->GetYaxis()->SetNdivisions(12, 0, kTRUE);
    h2_integral_z_phi->SetBit(TH1::kNoStats);

    std::vector<TLine *> UserGridXY;
    for (int j = 1; j <= 12; j++)
    {
        UserGridXY.push_back(new TLine(0.5, j + 0.5, 64.5, j + 0.5));
    }
    for (int i = 1; i <= 64; i++)
    {
        UserGridXY.push_back(new TLine(i + 0.5, 0.5, i + 0.5, 12.5));
    }

    TDirectory *WfDir = converter.outFile->mkdir("channel_wf");
    WfDir->cd();

    // while (converter.ReadEvent() && converter.EventNumber<1000){
    while (converter.ReadEvent())
    {

        if (targetEvent > 0 && converter.EventNumber != (uint32_t)targetEvent)
        {
            continue;
        }

        ChannelDataInEvent.clear();
        std::vector<int> _PhiCount(12);

        while (converter.ReadADC())
        {
            while (converter.ReadChannel())
            {

                // Считаем подложку и вычитем ее
                converter.channelEvent.pedestal = SetPedestal(converter.channelEvent.adcValues, iBinPedestalSkip, iBinPedestalCount);
                PedestalSubtraction(converter.channelEvent.adcValues, converter.channelEvent.pedestal);

                // Далее гистограмм для WF и подсчета интеграла
                int evNum = converter.channelEvent.eventNum;
                int chNum = converter.channelEvent.channelNums;
                int chPhi = converter.channelEvent.channel_Phi;
                int chZ = converter.channelEvent.channel_Z;
                std::string h_name = Form("h_ev_%i_ch_%i_phi_%i_z_%i", evNum, chNum, chPhi, chZ);
                std::string h_title = Form("Basket 38, Event %i, channel %i, Phi %i, Z %i", evNum, chNum, chPhi, chZ);

                SetWaveformHistogram(converter.channelEvent.adcValues, h1_wf, h_name, h_title);
                SetChannelIntegralAmp(converter.channelEvent, h1_wf, iBinStart, iBinStop, iBinLeft, iBinRight);

                if (converter.channelEvent.amplitude < 100)
                    continue;
                if (converter.channelEvent.integral < 500)
                    continue;
                _PhiCount.at(chPhi - 1)++;

                ChannelDataInEvent.push_back(converter.channelEvent);

                // int ch = converter.channelEvent.channelNums;
                // h_int_per_channel[ch - 1]->Fill(converter.channelEvent.integral);

                // converter.outFile->cd();
                // converter.FillTreeData();
            } // end ReadChannel
        } // end ReadADC

        WfDir->cd();

        if (ChannelDataInEvent.size() < 5)
            continue;
        if (ChannelDataInEvent.size() > 128)
            continue;

        int max_PhiCount = *std::max_element(_PhiCount.begin(), _PhiCount.end());
        if (max_PhiCount < 5)
            continue;

        // --- New 4.3: contamination check from neighboring channels ---
        double sumSelected = 0.0;
        double sumNeighbors = 0.0;

        // For now, treat all hits in ChannelDataInEvent as "target + selecting channels"
        double contamination = ComputeNeighborContamination(ChannelDataInEvent,
                                                            sumSelected,
                                                            sumNeighbors,
                                                            /*maxDeltaPhi=*/1,
                                                            /*maxDeltaZ=*/1);

        // User-configurable threshold_3 (e.g. 0.20 = 20%)
        const double threshold_3 = 0.20;

        if (contamination > threshold_3)
        {
            // Reject this event as multi-muon (or heavily contaminated) in neighboring cells
            continue;
        }

        for (int i = 0; i < ChannelDataInEvent.size(); i++)
        {
            converter.channelEvent = ChannelDataInEvent[i];
            converter.FillTreeData();
            h2_int_ch->Fill(converter.channelEvent.channelNums, converter.channelEvent.integral);

            int ch = converter.channelEvent.channelNums;
            if (ch >= 1 && ch <= NCHANNELS)
            {
                h_int_per_channel[ch - 1]->Fill(converter.channelEvent.integral);
            }
        }
        // DrawAllWaveformsOnOneCanvas(ChannelDataInEvent, eventHistograms,h2_integral_z_phi, UserGridXY, WfDir);

        // int ch = converter.channelEvent.channelNums;
        // h_int_per_channel[ch - 1]->Fill(converter.channelEvent.integral);
        if (targetEvent > 0)
        {
            break;
        }
    }

    converter.outFile->cd();

    // Create (or get) subdirectory for target-channel histograms
    TDirectory *dirTargets = nullptr;
    if (!(dirTargets = (TDirectory *)converter.outFile->Get("target_channels")))
    {
        dirTargets = converter.outFile->mkdir("target_channels");
    }

    // Move into that subdirectory
    dirTargets->cd();

    // Write all per-channel histograms there
    for (auto *h : h_int_per_channel)
    {
        if (!h)
            continue;
        h->Write();
    }

    h2_int_ch->Write();

    converter.outFile->cd();
    FitTargetChannelHistograms(converter.outFile);

    // ProgressBar(converter.inFile.tellg(), converter.inFile.tellg());
    std::cout << "\n\nOutput saved to: " << outputData << "\n\n"
              << std::endl;

    converter.WriteTreeAndClose();

    timer1.Stop();
    timer1.Print();
}