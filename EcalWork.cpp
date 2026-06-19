#include "MpdDataConverter.h"
#include <algorithm>
#include <cstddef>
#include <ostream>
#include <iostream>
#include <ctime>
#include <string>
#include <cstdlib>
#include <vector>
#include <set>
#include <filesystem>

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
#include "TTreeReader.h"
#include "TSystem.h"
#include "TF1.h"
#include "TGraphErrors.h"

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 0, 0)
    #ifdef __APPLE__
        R__LOAD_LIBRARY(./build/libMPDDataConverter.dylib)
    #else
        R__LOAD_LIBRARY(./build/libMPDDataConverter.so)
    #endif
#endif

const int NCHANNELS = 768;
constexpr int MAX_CH = 65;

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

double ComputeStripContamination(const std::vector<ChannelData> &data,
                                 const std::vector<int> &stripIndices,
                                 double &sumSelected,
                                 double &sumNeighbors,
                                 int maxDeltaPhi = 1,
                                 int maxDeltaZ = 1)
{
    sumSelected = 0.0;
    sumNeighbors = 0.0;

    if (data.empty() || stripIndices.empty())
        return 0.0;

    std::vector<bool> isStrip(data.size(), false);
    for (int idx : stripIndices)
    {
        if (idx >= 0 && idx < (int)data.size())
        {
            isStrip[idx] = true;
            sumSelected += data[idx].integral;
        }
    }

    std::set<int> neighborIndices;

    for (int iStrip : stripIndices)
    {
        const auto &ci = data[iStrip];

        for (int j = 0; j < (int)data.size(); ++j)
        {
            if (isStrip[j])
                continue;

            const auto &cj = data[j];

            int dPhi = std::abs(ci.channel_Phi - cj.channel_Phi);
            int dZ = std::abs(ci.channel_Z - cj.channel_Z);

            if (dPhi <= maxDeltaPhi && dZ <= maxDeltaZ)
                neighborIndices.insert(j);
        }
    }

    for (int idx : neighborIndices)
        sumNeighbors += data[idx].integral;

    const double denom = sumSelected + sumNeighbors;
    if (denom <= 0.0)
        return 0.0;

    return sumNeighbors / denom;
}
// Given all hits in the event (data), check if there exists one phi/Z row that contains
// at least minLen neighboring Z cells (e.g. Z = 10,11,12,13,14).
// If yes, return the indices of those hits
bool FindOneStrip(const std::vector<ChannelData> &data,
                  std::vector<int> &stripIndices,
                  bool alongPhi,
                  int minLen = 5)
{
    stripIndices.clear();

    // int outerMax = alongPhi ? 12 : 64;
    int outerMax = alongPhi ? 12 : 65;

    for (int outer = 1; outer <= outerMax; ++outer) // loop over each phi/Z row
    {
        std::vector<std::pair<int, int>> lineHits; // {innerCoord, index in data}

        for (int i = 0; i < (int)data.size(); ++i) // iterate through all hits in the event
        {
            int fixedCoord = alongPhi ? data[i].channel_Phi : data[i].channel_Z;
            int scanCoord = alongPhi ? data[i].channel_Z : data[i].channel_Phi;

            if (fixedCoord == outer) //  If a hit belongs to the current phi/Z row, add {phi/Z, i} to lineHits
                lineHits.push_back({scanCoord, i});
        }

        if ((int)lineHits.size() < minLen)
            continue;

        std::sort(lineHits.begin(), lineHits.end()); // sort the data to get ordered list like { (6, idxA), (8, idxB), (12, idxC), ... }

        std::vector<int> current; // current will store the indices in data of the current run of neighbouring phi/Z cells
        current.push_back(lineHits[0].second);

        for (int i = 1; i < (int)lineHits.size(); ++i) // looping through the rest of the hits in that row
        {
            if (lineHits[i].first == lineHits[i - 1].first + 1) //  ch_current == ch_previous + 1
            {
                current.push_back(lineHits[i].second);
            }
            else
            {
                current.clear();
                current.push_back(lineHits[i].second);
            }

            if ((int)current.size() >= minLen)
            {
                stripIndices = current;
                return true;
            }
        }
    }

    return false;
}

bool TransverseAnalysis(std::vector<ChannelData> &eventCh, std::vector<TH1F *> h_int_per_channel, TH1D *hCut)
{
    std::vector<int> _PhiCount(12);
    std::vector<int> _ZCount(64);
    std::vector<ChannelData> data;
    data.reserve(eventCh.size());

    for (int i = 0; i < eventCh.size(); i++)
    {

        if (eventCh[i].amplitude < 100)
            continue; // apm < 100 is a basic selection fro all channels in event
        hCut->Fill(2);

        if (eventCh[i].integral < 500)
            continue; // int < 500 is a basic selection  fro all channels in event
        hCut->Fill(3);

        data.push_back(eventCh[i]);              // we have got the 'target' channels in event
        _PhiCount[eventCh[i].channel_Phi - 1]++; // count the number of phi triggered cells in event
        _ZCount[eventCh[i].channel_Z - 1]++;     // count the number of Z triggered cells in event
    }

    if (data.size() < 5)
    { // if there are less tna five channels (doesn't matter which ones), then reject the event
        eventCh.clear();
        return false;
    }
    hCut->Fill(4, data.size());

    // This ones only check the number of triggered channels in row
    //--->
    /*
    int max_PhiCount = *std::max_element(_PhiCount.begin(), _PhiCount.end());
    int max_ZCount = *std::max_element(_ZCount.begin(), _ZCount.end());
    if (max_PhiCount < 5 && max_ZCount < 5) return false; //check if there a less than five channels in phi/Z row
    */
    //<---
    // Add the neighbors check here ---> poop from below
    // /*
    std::vector<int> phiStrip, zStrip, stripIndices;

    bool hasPhiStrip = FindOneStrip(data, phiStrip, true, 5);
    bool hasZStrip = FindOneStrip(data, zStrip, false, 5);

    if (!hasPhiStrip && !hasZStrip)
        return false;

    if (hasPhiStrip && hasZStrip)
        stripIndices = (phiStrip.size() >= zStrip.size()) ? phiStrip : zStrip;
    else if (hasPhiStrip)
        stripIndices = phiStrip;
    else
        stripIndices = zStrip;
    // */
    // if(data[0].eventNum<120)eventCh.back().Print();
    hCut->Fill(5, data.size());

    //  User-configurable threshold_3 (e.g. 0.20 = 20%)
    const double threshold_3 = 0.20;

    // /*
    // --- New 4.3: contamination check from neighboring channels ---
    double sumSelected = 0.0;
    double sumNeighbors = 0.0;
    double contamination = ComputeStripContamination(data,
                                                     stripIndices,
                                                     sumSelected,
                                                     sumNeighbors,
                                                     1, 1);
    if (contamination > threshold_3)
    {
        // Reject this event as multi-muon (or heavily contaminated) in neighboring cells
        eventCh.clear();
        return false;
    }
    // */
    hCut->Fill(6, data.size());

    for (int i = 0; i < data.size(); i++)
    {
        int ch = data[i].channelNums;
        if (ch >= 1 && ch <= NCHANNELS)
        {
            h_int_per_channel[ch - 1]->Fill(data[i].integral);
        }
    }


    // for (int i = 0; i < stripIndices.size(); i++)
    // {
    //     int ch = data[stripIndices[i]].channelNums;
    //     if (ch >= 1 && ch <= NCHANNELS)
    //     {
    //         h_int_per_channel[ch - 1]->Fill(data[stripIndices[i]].integral);
    //     }
    // }

    eventCh.clear();
    eventCh = std::move(data);

    return true;
}

// Selection 3x3 Neighbourhood window
// Finding the second strongest cell
bool CheckThreeOnThreeWindow(std::vector<ChannelData> &data, ChannelData hottestCell, float &adc_max_2)
{
    adc_max_2 = -1;
    float adc_max = hottestCell.integral;
    int phi_0 = hottestCell.channel_Phi;
    int z_0 = hottestCell.channel_Z;

    float min_ratio = 0.70; // dummy value
    float max_ratio = 0.95; // dummy value

    for (int phi = phi_0 - 1; phi <= phi_0 + 1; ++phi)
    {
        for (int z = z_0 - 1; z <= z_0 + 1; ++z)
        {

            // Skip the center cell itself
            if (phi == phi_0 && z == z_0)
                continue;

            // Search for any channel with these indices in this event
            for (const auto &ch : data)
            {
                if (ch.channel_Phi == phi && ch.channel_Z == z)
                {
                    if (ch.integral > adc_max_2)
                        adc_max_2 = ch.integral;
                }
            }
        }
    }

    float signifOfMax = adc_max / (adc_max + adc_max_2);
    if (signifOfMax > max_ratio || signifOfMax < min_ratio)
        return false;

    return true;
}

bool CalcAreaFiveOnFiveWindow(std::vector<ChannelData> &data,
                              const ChannelData &hottestCell,
                              float &sum5x5,
                              float &sum3x3,
                              float &ratio_cut5x5)
{
    sum5x5 = 0.0f;
    sum3x3 = 0.0f;
    ratio_cut5x5 = -1.0f;

    int phi_0 = hottestCell.channel_Phi;
    int z_0 = hottestCell.channel_Z;

    for (const auto &ch : data)
    {
        int dPhi = std::abs(ch.channel_Phi - phi_0);
        int dZ = std::abs(ch.channel_Z - z_0);

        // inside 5x5 window
        if (dPhi == 2 && dZ == 2){
            sum5x5 += ch.integral;
        }
        // inside 5x5 window
        if (dPhi == 1 && dZ == 1){
            sum3x3 += ch.integral;
        }
    }



    if (sum5x5 <= 0.0f)
        return false;

    // energy fraction outside the compact 3x3 core
    ratio_cut5x5 = (sum5x5 - sum3x3) / sum5x5;

    float max_diffusivity = 0.2f; // dummy value, tune from data
    if (ratio_cut5x5 > max_diffusivity)
        return false;

    return true;
}
bool LongAnalysis(std::vector<ChannelData> &eventCh, std::vector<TH1F *> h_int_per_channel, TH1D *hCut)
{
    std::vector<ChannelData> data;
    data.reserve(eventCh.size());

    for (int i = 0; i < eventCh.size(); i++)
    {
        if (eventCh[i].amplitude < 100)
            continue;
        hCut->Fill(2);
        if (eventCh[i].integral < 500)
            continue; // int < 500 is a basic selection  fro all channels in event
        hCut->Fill(3);

        data.push_back(eventCh[i]);
    }

    if (data.empty())
        return false; // nothing passed cuts, no maximum

    // find element with maximum integral
    // maximum by amplitude
    auto itMax = std::max_element(
        data.begin(), data.end(),
        [](const ChannelData &a, const ChannelData &b)
        {
            return a.integral < b.integral;
        });

    if (itMax == eventCh.end())
        return false;

    const ChannelData &maxCell = *itMax;
    float adc_max = maxCell.integral;
    // float adc_max_2 = -1;
    float adc_max_2;

    //int maxThreshold = 1000; // dummy value
    //int secondNoise = 50;    // dummy value

    float ration3on3 = 0;

    // 2. Reject if hottest cell is too small
    // if (adc_max < maxThreshold)
    //     return false;
    hCut->Fill(4,data.size());

    if (!CheckThreeOnThreeWindow(data, maxCell, adc_max_2))
        return false;
    hCut->Fill(5,data.size());

    if (adc_max_2 < 500)
        return false;
    hCut->Fill(6,data.size());

    float sum5x5, sum3x3, ratio_cut5x5;

    if(!CalcAreaFiveOnFiveWindow(data, maxCell, sum5x5, sum3x3, ratio_cut5x5)) 
        return false;

    hCut->Fill(7,data.size());

    for (int i = 0; i < data.size(); i++)
    {
        int ch = data[i].channelNums;
        if (ch >= 1 && ch <= NCHANNELS)
        {
            int phi_0 = maxCell.channel_Phi;
            int z_0 = maxCell.channel_Z;
            if(abs(phi_0-data[i].channel_Phi)<=1 && abs(z_0-data[i].channel_Z)<=1 ){
                h_int_per_channel[ch - 1]->Fill(data[i].integral);
                }    
            }
    }

    // std::cout << "Approx deposed energy = " << adc_max_2+adc_max  << "[ADC?]" <<std::endl;
    eventCh.clear();
    eventCh = std::move(data);

    return true;
}

void EcalWork(std::string inputDataTree = "out2.root", std::string outputData = "long_x_38.root")
{

    TStopwatch timer1;
    timer1.Start();

    bool TransverAnalysis = false;

    std::vector<TH1F *> h_int_per_channel;
    h_int_per_channel.reserve(NCHANNELS);
    for (int i = 0; i < NCHANNELS; ++i)
    {
        h_int_per_channel.push_back(
            new TH1F(Form("h_int_ch_%d", i + 1),
                     Form("Integral channel %d;Integral;Entries", i + 1),
                     100, 0, 15000));
    }

    TFile *outFile = new TFile(Form("%s", outputData.c_str()), "RECREATE");

    TFile *iFile = TFile::Open(inputDataTree.c_str());
    TTreeReader reader("events", iFile);

    TTreeReaderValue<UInt_t> eventNumber(reader, "eventNum");
    TTreeReaderValue<Int_t> chInt(reader, "integral");
    TTreeReaderValue<Int_t> chAmp(reader, "amplitude");
    TTreeReaderValue<unsigned short> chNum(reader, "channelNums");
    TTreeReaderValue<unsigned short> chZ(reader, "channel_Z");
    TTreeReaderValue<unsigned short> chPhi(reader, "channel_Phi");

    TH1D *h1_evNum = new TH1D("event_number", ";event;count", 2.0e6, 0, 2.0e6);
    TH1D *h1_chNum = new TH1D("channel_number", ";Num;count", NCHANNELS, 0.5, (Float_t)NCHANNELS + 0.5);
    TH1D *h1_chZ = new TH1D("channel_Z", ";Z;count", 64, 0.5, 64.5);
    TH1D *h1_chPhi = new TH1D("channel_Phi", ";Phi;count", 12, 0.5, 12.5);
    TH1D *h1_chAmp = new TH1D("channel_amplitude", ";Amp;count", 100, 0, 10000);
    TH1D *h1_chInt = new TH1D("channel_integral", ";Amp;count", 300, 0, 15000);

    std::filesystem::path outRootPath(outputData);
    std::filesystem::path outDir = outRootPath.parent_path();
    if (outDir.empty())
    {
        outDir = ".";
    }
    std::filesystem::path pictDir = outDir / "pict";
    std::filesystem::create_directories(pictDir);
    gSystem->ChangeDirectory(pictDir.string().c_str());
    std::cout << "Pictures will be saved under: " << pictDir << std::endl;

    TDirectory *WfDir = outFile->mkdir("channel_wf");
    WfDir->cd();

    UInt_t currentEventNum = 0;
    bool hasData = false;
    std::vector<ChannelData> ChannelDataInEvent;
    ChannelDataInEvent.reserve(NCHANNELS);

    TH1D *h_CountCut = new TH1D("hCountCut", "Number of events after applying cut;;count", 20, -0.5, 19.5);

    while (reader.Next())
    {
        ChannelData data;
        data.eventNum = *eventNumber;
        data.integral = *chInt;
        data.amplitude = *chAmp;
        data.channelNums = *chNum;
        data.channel_Z = *chZ;
        data.channel_Phi = *chPhi;

        h_CountCut->Fill(0.);

        // Проверка на конец дерева
        if (data.eventNum == 0)
        {
            data.eventNum = currentEventNum + 1; // Новый номер, чтобы вызвать обработку
        }

        if (!hasData)
        {
            currentEventNum = data.eventNum;
            hasData = true;
        }
        else if (data.eventNum != currentEventNum)
        { // Тут все действия с одним событием и каналами сработавшими в событии
            if (ChannelDataInEvent.size() <= MAX_CH)
            {
                h_CountCut->Fill(1, ChannelDataInEvent.size());

                if (TransverAnalysis)
                {
                    // std::cout << "Roflan pominki ... " << std::endl;

                    if (TransverseAnalysis(ChannelDataInEvent, h_int_per_channel, h_CountCut))
                    {
                        for (int i = 0; i < ChannelDataInEvent.size(); i++)
                        {
                            h1_evNum->Fill(ChannelDataInEvent[i].eventNum);
                            h1_chAmp->Fill(ChannelDataInEvent[i].amplitude);
                            h1_chInt->Fill(ChannelDataInEvent[i].integral);
                            h1_chNum->Fill(ChannelDataInEvent[i].channelNums);
                            h1_chPhi->Fill(ChannelDataInEvent[i].channel_Phi);
                            h1_chZ->Fill(ChannelDataInEvent[i].channel_Z);
                        }
                    }
                }
                else
                {
                    if (LongAnalysis(ChannelDataInEvent, h_int_per_channel, h_CountCut))
                    {
                        for (int i = 0; i < ChannelDataInEvent.size(); i++)
                        {
                            h1_evNum->Fill(ChannelDataInEvent[i].eventNum);
                            h1_chAmp->Fill(ChannelDataInEvent[i].amplitude);
                            h1_chInt->Fill(ChannelDataInEvent[i].integral);
                            h1_chNum->Fill(ChannelDataInEvent[i].channelNums);
                            h1_chPhi->Fill(ChannelDataInEvent[i].channel_Phi);
                            h1_chZ->Fill(ChannelDataInEvent[i].channel_Z);
                        }
                    }
                }
            }
            ChannelDataInEvent.clear();
            currentEventNum = data.eventNum;
        } // конец магии над событием

        // Если это была сторожевая запись - выходим (не добавляем её)
        if (*eventNumber == 0)
        {
            break;
        }

        ChannelDataInEvent.push_back(data);
    }

    outFile->cd();
    h1_evNum->Write();
    h_CountCut->Write();
    h1_chAmp->Write();
    h1_chInt->Write();
    h1_chNum->Write();
    h1_chPhi->Write();
    h1_chZ->Write();

    // Create (or get) subdirectory for target-channel histograms
    TDirectory *dirTargets = nullptr;
    if (!(dirTargets = (TDirectory *)outFile->Get("target_channels")))
    {
        dirTargets = outFile->mkdir("target_channels");
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

    outFile->cd();
    FitTargetChannelHistograms(outFile);

    // ProgressBar(converter.inFile.tellg(), converter.inFile.tellg());
    std::cout << "\n\nOutput saved to: " << outputData << "\n\n"
              << std::endl;

    timer1.Stop();
    timer1.Print();
    outFile->Close();
}

//Convert iFile.data to oFile.root without cuts
void ConvertToRoot( int targetEvent = -1,
                    std::vector<int> targetEvents={},
                    std::string inputData = "../run_rc-hs1_088.data",
                    std::string outputData = "out_all.root")
{

    TStopwatch timer1;
    timer1.Start();

    // Для поиска пика
    Int_t iBinStart = 20;
    Int_t iBinStop = 50;
    // Интеграл вокруг пика
    Int_t iBinLeft = 4;
    Int_t iBinRight = 16;

    InitChannelMap("basket_channel_map_phiZ.csv");

    MpdDataConverter converter;

    converter.OpenInputFile(inputData);
    converter.OpenOutRootFile(outputData);
    // converter.WriteChannelSamplesVector();
    converter.SetTypeWaveForms(TypeWaveForms::invert);
    converter.SetPedestalPar(/*pedestal_Skip*/ 2, /*pedestal_count*/ 10); // Число первых бинов для пропуска и число последующих бинов для подсчета подложки
    converter.InitTree();

    std::filesystem::path outRootPath(outputData);
    std::filesystem::path outDir = outRootPath.parent_path();
    if (outDir.empty())
    {
        outDir = ".";
    }
    std::filesystem::path pictDir = outDir / "pict";
    std::filesystem::create_directories(pictDir);
    gSystem->ChangeDirectory(pictDir.string().c_str());
    std::cout << "Pictures will be saved under: " << pictDir << std::endl;

    Int_t nBinsSamples = 60;
    Float_t SamplesMin = 0 - 0.5;
    Float_t SamplesMax = 60 - 0.5;

    TH1F *h1_wf = new TH1F("hist", "", nBinsSamples, SamplesMin, SamplesMax);

    std::vector<TH1F *> eventHistograms;
    eventHistograms.reserve(NCHANNELS);
    for (int i = 0; i < NCHANNELS; i++)
    {
        eventHistograms.push_back(new TH1F(Form("histo_%i", i), "", nBinsSamples, SamplesMin, SamplesMax));
    }

    std::vector<ChannelData> ChannelDataInEvent;
    ChannelDataInEvent.reserve(NCHANNELS);

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

        if (!targetEvents.empty())
        {
            if (std::find(targetEvents.begin(), targetEvents.end(), (uint32_t)targetEvent) != targetEvents.end())
            {
                continue;
            }
        }

        ChannelDataInEvent.clear();

        while (converter.ReadADC())
        {
            while (converter.ReadChannel())
            {
                // Далее гистограмм для WF и подсчета интеграла
                int evNum = converter.channelEvent.eventNum;
                int chNum = converter.channelEvent.channelNums;
                int chPhi = converter.channelEvent.channel_Phi;
                int chZ = converter.channelEvent.channel_Z;
                std::string h_name = Form("h_ev_%i_ch_%i_phi_%i_z_%i", evNum, chNum, chPhi, chZ);
                std::string h_title = Form("Basket 38, Event %i, channel %i, Phi %i, Z %i", evNum, chNum, chPhi, chZ);

                SetWaveformHistogram(converter.channelEvent.adcValues, h1_wf, h_name, h_title);
                SetChannelIntegralAmp(converter.channelEvent, h1_wf, iBinStart, iBinStop, iBinLeft, iBinRight);
                ChannelDataInEvent.push_back(converter.channelEvent);
                converter.FillTreeData();
            } // end ReadChannel
        } // end ReadADC

        WfDir->cd();

        //if(ChannelDataInEvent.size()>64)

        if(converter.EventNumber>1500 && converter.EventNumber<2000)
            DrawAllWaveformsOnOneCanvas(ChannelDataInEvent, eventHistograms,h2_integral_z_phi, UserGridXY, WfDir);

        if(targetEvent > 0 || targetEvents.empty()==false){
            DrawAllWaveformsOnOneCanvas(ChannelDataInEvent, eventHistograms,h2_integral_z_phi, UserGridXY, WfDir);
            if (converter.EventNumber >= targetEvents.back()){
                break;
            }
            if (targetEvent > 0)
            {
                break;
            }
        }
    }

    // Запись последнего ивента со всеми значенями равными 0.
    converter.channelEvent.Clear();
    converter.FillTreeData();

    converter.outFile->cd();

    std::cout << "\nOutput Tree saved to: " << outputData << "\n"
              << std::endl;

    converter.WriteTreeAndClose();

    timer1.Stop();
    timer1.Print();
}