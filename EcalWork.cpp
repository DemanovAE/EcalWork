#include "EcalDrawClass.h"
#include "MpdDataConverter.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMathBase.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTreeReader.h"

// clang-format off
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 0, 0)
#ifdef __APPLE__
R__LOAD_LIBRARY(./build/libMPDDataConverter.dylib)
#else
R__LOAD_LIBRARY(./build/libMPDDataConverter.so)
#endif
#endif
// clang-format on

const int NCHANNELS = 768;
constexpr int MAX_CH = 65;

const int MAX_Z = 65;
const int MAX_PHI = 12;

bool SetBinNameCutHisto = true;

// For Long Analysis
bool TransverAnalysis = false;
const int gl_long_max_1st_Integral = 1200;
const int gl_long_max_2nd_Integral = 500;
const float gl_long_min_3x3_ratio = 0.70;       // dummy value
const float gl_long_max_3x3_ratio = 0.99;       // dummy value
const float gl_long_max_diffusivity_5x5 = 0.45; // dummy value, tune from data

// For Trans Analysis

template <typename T>
void SetWaveformHistogram(const std::vector<T> &data, TH1F *hist,
                          const std::string &name = "waveform",
                          const std::string &title = "",
                          const std::string &xTitle = "Sample Number",
                          const std::string &yTitle = "ADC") {
  if (data.empty())
    return;

  int nSamples = data.size();

  // TH1F* hist = new TH1F(name.c_str(), title.c_str(), nSamples, 0, nSamples);

  hist->SetName(name.c_str());
  hist->SetTitle(title.c_str());

  for (int i = 0; i < nSamples; i++) {
    hist->SetBinContent(i + 1, static_cast<Float_t>(data[i]));
  }

  hist->GetXaxis()->SetTitle(xTitle.c_str());
  hist->GetYaxis()->SetTitle(yTitle.c_str());

  hist->SetLineColor(kBlue);
  hist->SetLineWidth(2);
  hist->SetDrawOption("L");
  hist->SetBit(TH1::kNoStats);
}

void SetChannelIntegralAmp(ChannelData &_chData, TH1F *hist, int binStart,
                           int binEnd, int nBinLeft, int nBinRight) {
  if (!hist)
    return;
  if (hist->GetEntries() == 0)
    return;
  // Пик
  Float_t maxValue = 0;
  Int_t peakBin = 0;

  for (int i = binStart; i <= binEnd; i++) {
    Float_t currentValue = TMath::Abs(hist->GetBinContent(i));
    if (currentValue > maxValue) {
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

void FitTargetChannelHistograms(
    TFile *outFile, const std::string &inputDirName = "target_channels",
    const std::string &outputDirName = "fitted_target_channels",
    int minEntriesToFit = 50, int rebinFactor = 2) {
  if (!outFile || outFile->IsZombie()) {
    std::cout << "FitTargetChannelHistograms: output file is null or zombie"
              << std::endl;
    return;
  }

  TDirectory *dirIn = (TDirectory *)outFile->Get(inputDirName.c_str());
  if (!dirIn) {
    std::cout << "FitTargetChannelHistograms: directory " << inputDirName
              << " not found" << std::endl;
    return;
  }

  outFile->cd();
  TDirectory *dirOut = (TDirectory *)outFile->Get(outputDirName.c_str());
  if (!dirOut) {
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

  for (int ch = 1; ch <= NCHANNELS; ++ch) {
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
    if (nBins < 10) {
      delete hFit;
      continue;
    }

    // Find the global maximum bin after rebinning
    int maxBin = hFit->GetMaximumBin();
    double peakX = hFit->GetBinCenter(maxBin);
    double peakY = hFit->GetBinContent(maxBin);

    if (peakY <= 0) {
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

    // Avoid fitting the low-integral pedestal if the peak is well away from
    // zero
    if (peakX > 1500.0 && fitMin < 0.5 * peakX)
      fitMin = 0.5 * peakX;

    if (fitMax <= fitMin) {
      delete hFit;
      continue;
    }

    TF1 *fLan = new TF1(Form("f_landau_ch_%d", ch), "landau", fitMin, fitMax);
    fLan->SetParameters(peakY, peakX, widthGuess);
    fLan->SetParLimits(1, 4000.0, 9000.0); // MPV between 3k and 9k
    fLan->SetLineColor(kRed);
    fLan->SetLineWidth(2);

    int fitStatus = hFit->Fit(fLan, "QRS");
    if (fitStatus != 0) {
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

    if (badFit) {
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

  if (nFitted > 0) {
    TGraphErrors *grMPV = new TGraphErrors(nFitted, v_ch.data(), v_mpv.data(),
                                           nullptr, v_mpv_err.data());
    grMPV->SetName("gr_mpv_vs_channel");
    grMPV->SetTitle("Landau MPV vs channel;Channel;Landau MPV");
    grMPV->Write();

    TGraphErrors *grWidth = new TGraphErrors(
        nFitted, v_ch.data(), v_width.data(), nullptr, v_width_err.data());
    grWidth->SetName("gr_width_vs_channel");
    grWidth->SetTitle("Landau width vs channel;Channel;Landau width");
    grWidth->Write();

    TGraph *grEntries = new TGraph(nFitted, v_ch.data(), v_entries.data());
    grEntries->SetName("gr_entries_vs_channel");
    grEntries->SetTitle("Entries vs channel;Channel;Entries");
    grEntries->Write();
  }

  outFile->cd();

  std::cout << "FitTargetChannelHistograms: fitted " << nFitted
            << " channel histograms with Landau into directory "
            << outputDirName << std::endl;
}

bool CheckContamination(const std::vector<ChannelData> &data,
                        const std::vector<ChannelData> &strip, bool alongPhi,
                        double minIntegralForNoise, int maxContaminatingCells,
                        int neighborhoodRadiusZ = 1,
                        int neighborhoodRadiusPhi = 1) {
  // 1. Build a set of strip coordinates
  std::set<std::pair<int, int>> stripCoords;
  for (const auto &ch : strip) {
    stripCoords.insert({(int)ch.channel_Z, (int)ch.channel_Phi});
  }

  int contaminationCount = 0;

  // 2. For each channel in the event
  for (const auto &ch : data) {
    int z = (int)ch.channel_Z;
    int phi = (int)ch.channel_Phi;

    // Ignore tiny hits
    if (ch.integral < minIntegralForNoise)
      continue;

    // Skip if it's part of the strip
    if (stripCoords.count({z, phi}))
      continue;

    // 3. Check if it is close to ANY strip cell
    bool nearStrip = false;
    for (const auto &s : stripCoords) {
      int sz = s.first;
      int sphi = s.second;
      if (std::abs(z - sz) <= neighborhoodRadiusZ &&
          std::abs(phi - sphi) <= neighborhoodRadiusPhi) {
        nearStrip = true;
        break;
      }
    }

    if (nearStrip) {
      ++contaminationCount;
      if (contaminationCount > maxContaminatingCells)
        return false; // too much contamination
    }
  }

  // 4. Accept if contamination below threshold
  return true;
}

// Given all hits in the event (data), check if there exists one phi/Z row that
// contains at least minLen neighboring Z cells (e.g. Z = 10,11,12,13,14). If
// yes, return the indices of those hits
bool FindOneStrip(const std::vector<ChannelData> &data,
                  std::vector<int> &stripIndices, int &foundOuter,
                  bool alongPhi, double threshold_2, int minLen = 5) {

  stripIndices.clear();
  foundOuter = -1;
  int outerMax = alongPhi ? MAX_PHI : MAX_Z;
  int innerMax = alongPhi ? MAX_Z : MAX_PHI;

  for (int outer = 1; outer <= outerMax; ++outer) {
    std::vector<int> acceptedCells;
    acceptedCells.reserve(innerMax);

    for (int inner = 1; inner <= innerMax; ++inner) {
      for (int i = 0; i < (int)data.size(); ++i) {
        if (alongPhi ? (inner == (int)data[i].channel_Z &&
                        outer == (int)data[i].channel_Phi)
                     : (inner == (int)data[i].channel_Phi &&
                        outer == (int)data[i].channel_Z)) {
          acceptedCells.push_back(alongPhi ? (int)data[i].channel_Z
                                           : (int)data[i].channel_Phi);
          break;
        }
      }
    }

    if ((int)acceptedCells.size() < minLen)
      continue;

    int runLen = 1;
    int bestRunLen = 0;
    int bestRunStart = 0;

    for (int i = 1; i < (int)acceptedCells.size(); ++i) {
      if (acceptedCells[i] == acceptedCells[i - 1] + 1) {
        ++runLen;
      } else {
        runLen = 1;
      }
      if (runLen > bestRunLen) {
        bestRunLen = runLen;
        bestRunStart = i - runLen + 1;
      }
    }
    if (bestRunLen >= minLen) {
      stripIndices.clear();
      for (int k = 0; k < bestRunLen; ++k)
        stripIndices.push_back(acceptedCells[bestRunStart + k]);
      foundOuter = outer;
      return true;
      // }
    }
  }
  return false;
}

// Find a target channel in a strip: center of 5 consecutive channels
bool FindTargetChannelInStrip(const std::vector<ChannelData> &strip,
                              bool alongPhi, double threshold_2,
                              ChannelData &target) {
  // Need at least 5 channels to form a consecutive group
  if (strip.size() < 5)
    return false;

  // Work on a copy of the input strip
  std::vector<ChannelData> ordered = strip;

  // Sort channels either by Z or by Phi coordinate
  std::sort(ordered.begin(), ordered.end(),
            [alongPhi](const ChannelData &a, const ChannelData &b) {
              return alongPhi ? (a.channel_Z < b.channel_Z)
                              : (a.channel_Phi < b.channel_Phi);
            });

  // Try each channel as the center of a group of 5
  for (size_t i = 2; i + 2 < ordered.size(); ++i) {
    // Check that 5 channels around i are consecutive in Z or Phi
    bool consecutive = (alongPhi ? ((int)ordered[i - 2].channel_Z + 1 ==
                                        (int)ordered[i - 1].channel_Z &&
                                    (int)ordered[i - 1].channel_Z + 1 ==
                                        (int)ordered[i].channel_Z &&
                                    (int)ordered[i].channel_Z + 1 ==
                                        (int)ordered[i + 1].channel_Z &&
                                    (int)ordered[i + 1].channel_Z + 1 ==
                                        (int)ordered[i + 2].channel_Z)
                                 : ((int)ordered[i - 2].channel_Phi + 1 ==
                                        (int)ordered[i - 1].channel_Phi &&
                                    (int)ordered[i - 1].channel_Phi + 1 ==
                                        (int)ordered[i].channel_Phi &&
                                    (int)ordered[i].channel_Phi + 1 ==
                                        (int)ordered[i + 1].channel_Phi &&
                                    (int)ordered[i + 1].channel_Phi + 1 ==
                                        (int)ordered[i + 2].channel_Phi));

    // If not consecutive, skip this center
    if (!consecutive)
      continue;

    // Check that the four neighbors around the center are above the threshold
    bool neighborsOK = ordered[i - 2].integral > threshold_2 &&
                       ordered[i - 1].integral > threshold_2 &&
                       ordered[i + 1].integral > threshold_2 &&
                       ordered[i + 2].integral > threshold_2;

    // If neighbors are strong enough, select the center as target
    if (neighborsOK) {
      target = ordered[i];
      return true;
    }
  }

  // No suitable target found
  return false;
}

bool TransverseAnalysis(std::vector<ChannelData> &eventCh,
                        std::vector<TH1F *> h_int_per_channel, TH1D *hCut) {
  std::vector<int> _PhiCount(12);
  std::vector<int> _ZCount(64);
  std::vector<ChannelData> data;
  data.reserve(eventCh.size());

  const double threshold_1 = 100.0;
  const double threshold_2 = 500.0;

  for (int i = 0; i < eventCh.size(); i++) {

    if (eventCh[i].amplitude < threshold_1)
      continue; // apm < 100 is a basic selection fro all channels in event
    hCut->Fill(2);

    data.push_back(eventCh[i]); // we have got the 'target' channels in event
    _PhiCount[eventCh[i].channel_Phi -
              1]++; // count the number of phi triggered cells in event
    _ZCount[eventCh[i].channel_Z -
            1]++; // count the number of Z triggered cells in event
  }

  if (data.size() < 5) { // if there are less tna five channels (doesn't
                         // matter which ones), then reject the event
    eventCh.clear();
    return false;
  }
  hCut->Fill(3, data.size());

  std::vector<int> phiStrip, zStrip;

  int phiStripOuter = -1, zStripOuter = -1;
  bool hasPhiStrip =
      FindOneStrip(data, phiStrip, phiStripOuter, true, threshold_2, 5);
  bool hasZStrip =
      FindOneStrip(data, zStrip, zStripOuter, false, threshold_2, 5);

  bool accepted = false;
  std::vector<ChannelData> filteredData;
  if (hasPhiStrip) {
    accepted = true;
    // оставляем только ячейки с нужным phi=phiStripOuter И Z из phiStrip
    for (auto &ch : data) {
      bool outerOk = (int)ch.channel_Phi == phiStripOuter;
      bool innerOk = std::find(phiStrip.begin(), phiStrip.end(),
                               (int)ch.channel_Z) != phiStrip.end();
      if (outerOk && innerOk)
        filteredData.push_back(ch);
    }
  }

  if (hasZStrip && !accepted) {
    accepted = true;
    // оставляем только ячейки с нужным Z=zStripOuter И Phi из zStrip
    for (auto &ch : data) {
      bool outerOk = (int)ch.channel_Z == zStripOuter;
      bool innerOk = std::find(zStrip.begin(), zStrip.end(),
                               (int)ch.channel_Phi) != zStrip.end();
      if (outerOk && innerOk)
        filteredData.push_back(ch);
    }
  }
  if (!accepted) {
    eventCh.clear();
    return false;
  }

  // if(data[0].eventNum<120)eventCh.back().Print();
  hCut->Fill(4, data.size());

  //  User-configurable threshold_3 (e.g. 0.20 = 20%)
  const double threshold_3 = 0.20;

  // 4. Contamination check
  double minIntegralForNoise =
      threshold_2 * 0.2; // example: weaker threshold for counting contamination
  int maxContaminatingCells = 2; // allow up to 2 near off-strip hits
  int neighborhoodRadiusZ = 1;
  int neighborhoodRadiusPhi = 1;
  if (hasZStrip && !hasPhiStrip) {
    if (!CheckContamination(data, filteredData, false, minIntegralForNoise,
                            maxContaminatingCells, neighborhoodRadiusZ,
                            neighborhoodRadiusPhi)) {
      eventCh.clear();
      return false;
    }
  }
  hCut->Fill(5, data.size());

  for (int i = 0; i < data.size(); i++) {
    int ch = data[i].channelNums;
    if (ch >= 1 && ch <= NCHANNELS) {
      h_int_per_channel[ch - 1]->Fill(data[i].integral);
    }
  }

  if (data[0].eventNum > 1700 && data[0].eventNum < 2000) {
    // if (data[0].eventNum > 1700 && data[0].eventNum < 1800) {

    EcalDrawClass drawObj;
    drawObj.InitCanvas1Pad();
    drawObj.UpdateAndSaveCanvas("", data);
    // drawObj.UpdateAndSaveCanvas("_cut", data, stripIndices);
    drawObj.UpdateAndSaveCanvas(
        "_cut", filteredData); // ← передаём filteredData, не stripIndices
  }

  ChannelData targetChannel;
  bool hasTarget = false;

  if (hasPhiStrip) {
    hasTarget = FindTargetChannelInStrip(filteredData, true, threshold_2,
                                         targetChannel);
  } else if (hasZStrip) {
    hasTarget = FindTargetChannelInStrip(filteredData, false, threshold_2,
                                         targetChannel);
  }

  if (!hasTarget) {
    eventCh.clear();
    return false;
  }

  eventCh.clear();
  // eventCh = std::move(data);
  // eventCh = std::move(filteredData);
  eventCh.push_back(targetChannel);

  return true;
}

// Selection 3x3 Neighbourhood window
// Finding the second strongest cell
bool CheckThreeOnThreeWindow(std::vector<ChannelData> &data,
                             ChannelData hottestCell,
                             ChannelData &secondHotCell, int &num_i,
                             float &signifOfMax) {
  // adc_max_2 = -1;
  num_i = -1;
  // to prevent unpredictable behavior
  secondHotCell = hottestCell;
  secondHotCell.integral = -999999;

  int phi_0 = hottestCell.channel_Phi;
  int z_0 = hottestCell.channel_Z;

  for (int ch = 0; ch < static_cast<int>(data.size()); ++ch) {
    int dPhi = std::abs(data[ch].channel_Phi - phi_0);
    int dZ = std::abs(data[ch].channel_Z - z_0);
    if (data[ch].channel_Phi == phi_0 && data[ch].channel_Z == z_0)
      continue;
    if (dPhi <= 1 && dZ <= 1) {
      if (data[ch].integral > secondHotCell.integral) {
        secondHotCell = data[ch];
        num_i = ch;
      }
    }
  }

  float adc_max = hottestCell.integral;
  float adc_max_2 = secondHotCell.integral;
  // std::cout<<"adc_max_2 "<<adc_max_2<<std::endl;
  if (num_i < 0 || adc_max_2 <= 0.0f)
    return false;

  // float signifOfMax = adc_max / (adc_max + adc_max_2);
  signifOfMax = adc_max / (adc_max + adc_max_2);

  return true;
}

bool CalcAreaFiveOnFiveWindow(std::vector<ChannelData> &data,
                              const ChannelData &hottestCell, float &sum5x5,
                              float &sum3x3, float &ratio_cut5x5) {
  sum5x5 = 0.0f;
  sum3x3 = 0.0f;
  ratio_cut5x5 = -1.0f;

  int phi_0 = hottestCell.channel_Phi;
  int z_0 = hottestCell.channel_Z;

  for (const auto &ch : data) {
    int dPhi = std::abs(ch.channel_Phi - phi_0);
    int dZ = std::abs(ch.channel_Z - z_0);

    if (dPhi <= 2 && dZ <= 2)
      sum5x5 += ch.integral;

    if (dPhi <= 1 && dZ <= 1)
      sum3x3 += ch.integral;
  }

  if (sum5x5 <= 0.0f)
    return false;

  float e_outside_core = sum5x5 - sum3x3;
  if (e_outside_core < 0.0f)
    e_outside_core = 0.0f;

  ratio_cut5x5 = e_outside_core / sum5x5;
  return true;
}

bool LongAnalysis(std::vector<ChannelData> &eventCh,
                  std::vector<TH1F *> h_int_per_channel, TH1D *hCut) {

  std::vector<ChannelData> data;
  data.reserve(eventCh.size());
  int itMax2 = -1;

  for (int i = 0; i < eventCh.size(); i++) {
    // if (eventCh[i].amplitude < 200)
    //   continue;
    hCut->Fill(2);
    if (SetBinNameCutHisto)
      hCut->GetXaxis()->SetBinLabel(3, "Nope");
    // if (eventCh[i].integral < 500)
    //   continue; // int < 500 is a basic selection  fro all channels in event
    // hCut->Fill(3);

    data.push_back(eventCh[i]);
  }

  if (data.empty())
    return false; // nothing passed cuts, no maximum

  hCut->Fill(3, data.size());
  if (SetBinNameCutHisto)
    hCut->GetXaxis()->SetBinLabel(4, "empty event");

  int MaxIntegral_1st = -9999999;
  int itMax = -1;
  for (int iCh = 0; iCh < data.size(); iCh++) {
    if (data[iCh].integral > MaxIntegral_1st) {
      MaxIntegral_1st = data[iCh].integral;
      itMax = iCh;
    }
  }

  const ChannelData &maxCell = data[itMax];
  float adc_max = maxCell.integral;
  // float adc_max_2;

  if (adc_max < gl_long_max_1st_Integral)
    return false;

  hCut->Fill(4, data.size());
  if (SetBinNameCutHisto)
    hCut->GetXaxis()->SetBinLabel(
        5, Form("1stMaxInt>%i", gl_long_max_1st_Integral));

  if (maxCell.channel_Phi == 1 || maxCell.channel_Phi == 12)
    return false;

  // Need full 5x5 around hottest cell
  // if (maxCell.channel_Phi <= 2 || maxCell.channel_Phi >= 11)
  //   return false;

  hCut->Fill(5, data.size());
  if (SetBinNameCutHisto)
    hCut->GetXaxis()->SetBinLabel(6, "Phi==1;12");

  if (maxCell.channel_Z == 1 || maxCell.channel_Z == 64)
    return false;

  // if (maxCell.channel_Z <= 2 || maxCell.channel_Z >= 63)
  //   return false;

  hCut->Fill(6, data.size());
  if (SetBinNameCutHisto)
    hCut->GetXaxis()->SetBinLabel(7, "Z==1;64");

  float ratio_3x3 = 0;
  ChannelData secMaxCell;
  if (!CheckThreeOnThreeWindow(data, maxCell, secMaxCell, itMax2, ratio_3x3))
    return false;

  if (ratio_3x3 >= gl_long_max_3x3_ratio || ratio_3x3 < gl_long_min_3x3_ratio)
    return false;

  hCut->Fill(7, data.size());
  if (SetBinNameCutHisto)
    hCut->GetXaxis()->SetBinLabel(
        8, Form("1st/(1st+2nd)>%.2f", gl_long_min_3x3_ratio));

  float MaxIntegral_2st = secMaxCell.integral;
  if (MaxIntegral_2st < gl_long_max_2nd_Integral)
    return false;

  hCut->Fill(8, data.size());
  if (SetBinNameCutHisto)
    hCut->GetXaxis()->SetBinLabel(
        9, Form("2ndMaxInt>%i", gl_long_max_2nd_Integral));

  float sum5x5, sum3x3, ratio_cut5x5;

  if (!CalcAreaFiveOnFiveWindow(data, maxCell, sum5x5, sum3x3, ratio_cut5x5))
    return false;

  if (ratio_cut5x5 > gl_long_max_diffusivity_5x5)
    return false;
  hCut->Fill(9, data.size());
  if (SetBinNameCutHisto) {
    hCut->GetXaxis()->SetBinLabel(
        10, Form("(5x5)/(3x3+5x5)<%.2f", gl_long_max_diffusivity_5x5));
    SetBinNameCutHisto = false;
  }

  // Fill target channels directly
  // int ch_max = maxCell.channelNums;
  // if (ch_max >= 1 && ch_max <= NCHANNELS && h_int_per_channel[ch_max - 1]) {
  //   h_int_per_channel[ch_max - 1]->Fill(maxCell.integral);
  // }

  // int ch_second = secMaxCell.channelNums;
  // if (ch_second >= 1 && ch_second <= NCHANNELS &&
  //     h_int_per_channel[ch_second - 1]) {
  //   h_int_per_channel[ch_second - 1]->Fill(secMaxCell.integral);
  // }

  // Longitudinal observable: sum of hottest and strongest neighbour
  const float coreEnergy = maxCell.integral + secMaxCell.integral;

  // Store in histogram of the hottest channel
  const int ch_core = maxCell.channelNums;
  if (ch_core >= 1 && ch_core <= NCHANNELS && h_int_per_channel[ch_core - 1]) {
    h_int_per_channel[ch_core - 1]->Fill(coreEnergy);
  }
  std::vector<ChannelData> FinalCh = {maxCell, secMaxCell};

  if (data[0].eventNum > 10000 && data[0].eventNum < 10100) {
    EcalDrawClass drawObj;
    drawObj.InitCanvas1Pad();
    drawObj.UpdateAndSaveCanvas("", data);
    drawObj.UpdateAndSaveCanvas("_cut", FinalCh);
  }
  // std::cout << "Approx deposed energy = " << adc_max_2+adc_max  << "[ADC?]"
  // <<std::endl;
  // eventCh.clear();
  // FinalCh.pop_back();
  // FinalCh.back().integral = coreEnergy;
  // eventCh = std::move(FinalCh);

  // eventCh.clear();
  // eventCh = std::move(FinalCh);

    eventCh.clear();
  ChannelData outCell = maxCell;
  outCell.integral = coreEnergy;
  // put that single cell into eventCh
  eventCh.push_back(outCell);

  return true;
}

void EcalWork(std::string inputDataTree = "out_all.root",
              std::string outputData = "basket_38.root") {

  TStopwatch timer1;
  timer1.Start();

  std::vector<TH1F *> h_int_per_channel;
  h_int_per_channel.reserve(NCHANNELS);
  for (int i = 0; i < NCHANNELS; ++i) {
    h_int_per_channel.push_back(
        new TH1F(Form("h_int_ch_%d", i + 1),
                 Form("Integral channel %d;Integral;Entries", i + 1), 100, 0,
                 TransverAnalysis == true ? 15000 : 50000));
  }

  TFile *outFile =
      new TFile(Form("%s_%s", TransverAnalysis == true ? "trans" : "long",
                     outputData.c_str()),
                "RECREATE");

  TFile *iFile = TFile::Open(inputDataTree.c_str());
  TTreeReader reader("events", iFile);

  std::cout << "The number of entries in decoded file '"
            << inputDataTree.c_str() << "' " << reader.GetEntries()
            << std::endl;
  TTreeReaderValue<UInt_t> eventNumber(reader, "eventNum");
  TTreeReaderValue<Int_t> chInt(reader, "integral");
  TTreeReaderValue<Int_t> chAmp(reader, "amplitude");
  TTreeReaderValue<unsigned short> chNum(reader, "channelNums");
  TTreeReaderValue<unsigned short> chZ(reader, "channel_Z");
  TTreeReaderValue<unsigned short> chPhi(reader, "channel_Phi");

  TH1D *h1_evNum = new TH1D("event_number", ";event;count", 2.0e6, 0, 2.0e6);
  TH1D *h1_chNum = new TH1D("channel_number", ";Num;count", NCHANNELS, 0.5,
                            (Float_t)NCHANNELS + 0.5);
  TH1D *h1_chZ = new TH1D("channel_Z", ";Z;count", 64, 0.5, 64.5);
  TH1D *h1_chPhi = new TH1D("channel_Phi", ";Phi;count", 12, 0.5, 12.5);
  TH1D *h1_chAmp = new TH1D("channel_amplitude", ";Amp;count", 100, 0,
                            TransverAnalysis == true ? 10000 : 30000);
  TH1D *h1_chInt = new TH1D("channel_integral", ";Amp;count", 500, 0,
                            // TransverAnalysis == true ? 15000 : 50000);
                            TransverAnalysis == true ? 15000 : 1.e5);

  std::filesystem::path outRootPath(outputData);
  std::filesystem::path outDir = outRootPath.parent_path();
  if (outDir.empty()) {
    outDir = ".";
  }

  std::filesystem::path pictDir =
      outDir / "pict" / (TransverAnalysis ? "transverse" : "longitudinal");
  std::filesystem::create_directories(pictDir);
  gSystem->ChangeDirectory(pictDir.string().c_str());
  std::cout << "Pictures will be saved under: " << pictDir << std::endl;

  TDirectory *WfDir = outFile->mkdir("channel_wf");
  WfDir->cd();

  UInt_t currentEventNum = 0;
  bool hasData = false;
  std::vector<ChannelData> ChannelDataInEvent;
  ChannelDataInEvent.reserve(NCHANNELS);

  TH1D *h_CountCut =
      new TH1D("hCountCut", "Number of events after applying cut;;count", 20,
               -0.5, 19.5);

  while (reader.Next()) {
    ChannelData data;
    data.eventNum = *eventNumber;
    data.integral = *chInt;
    data.amplitude = *chAmp;
    data.channelNums = *chNum;
    data.channel_Z = *chZ;
    data.channel_Phi = *chPhi;

    h_CountCut->Fill(0.);
    h_CountCut->GetXaxis()->SetBinLabel(1, "def");

    // Проверка на конец дерева
    if (data.eventNum == 0) {
      data.eventNum =
          currentEventNum + 1; // Новый номер, чтобы вызвать обработку
    }

    if (!hasData) {
      currentEventNum = data.eventNum;
      hasData = true;
    } else if (data.eventNum !=
               currentEventNum) { // Тут все действия с одним событием и
                                  // каналами сработавшими в событии
      if (ChannelDataInEvent.size() <= MAX_CH) {
        h_CountCut->Fill(1, ChannelDataInEvent.size());
        h_CountCut->GetXaxis()->SetBinLabel(2, "Ch.Size()<65");

        bool KeyChData = false;

        if (TransverAnalysis) {
          KeyChData = TransverseAnalysis(ChannelDataInEvent, h_int_per_channel,
                                         h_CountCut);
        } else {
          KeyChData =
              LongAnalysis(ChannelDataInEvent, h_int_per_channel, h_CountCut);
        }

        if (KeyChData == true) {
          for (int i = 0; i < ChannelDataInEvent.size(); i++) {
            h1_evNum->Fill(ChannelDataInEvent[i].eventNum);
            h1_chAmp->Fill(ChannelDataInEvent[i].amplitude);
            h1_chInt->Fill(ChannelDataInEvent[i].integral);
            h1_chNum->Fill(ChannelDataInEvent[i].channelNums);
            h1_chPhi->Fill(ChannelDataInEvent[i].channel_Phi);
            h1_chZ->Fill(ChannelDataInEvent[i].channel_Z);
          }
        }
      }
      ChannelDataInEvent.clear();
      currentEventNum = data.eventNum;
    } // конец магии над событием

    // Если это была сторожевая запись - выходим (не добавляем её)
    if (*eventNumber == 0) {
      break;
    }

    ChannelDataInEvent.push_back(data);
  }

  // h_CountCut->Scale(100./h_CountCut->GetBinContent(1));

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
  if (!(dirTargets = (TDirectory *)outFile->Get("target_channels"))) {
    dirTargets = outFile->mkdir("target_channels");
  }

  // Move into that subdirectory
  dirTargets->cd();

  // Write all per-channel histograms there
  for (auto *h : h_int_per_channel) {
    if (!h)
      continue;
    h->Write();
  }

  outFile->cd();
  FitTargetChannelHistograms(outFile);

  // ProgressBar(converter.inFile.tellg(), converter.inFile.tellg());
  std::cout << "\n\nOutput saved to: " << outFile->GetName() << "\n\n"
            << std::endl;

  timer1.Stop();
  timer1.Print();
  outFile->Close();
}

// Convert iFile.data to oFile.root without cuts
void ConvertToRoot(int targetEvent = -1, std::vector<int> targetEvents = {},
                   std::string inputData = "../run_rc-hs1_088.data",
                   std::string outputData = "out_all.root") {

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
  converter.SetPedestalPar(
      /*pedestal_Skip*/ 2,
      /*pedestal_count*/ 10); // Число первых бинов для пропуска и число
                              // последующих бинов для подсчета подложки
  converter.InitTree();

  std::filesystem::path outRootPath(outputData);
  std::filesystem::path outDir = outRootPath.parent_path();
  if (outDir.empty()) {
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

  std::vector<ChannelData> ChannelDataInEvent;
  ChannelDataInEvent.reserve(NCHANNELS);

  TDirectory *WfDir = converter.outFile->mkdir("channel_wf");
  WfDir->cd();

  EcalDrawClass drawObj;
  drawObj.InitCanvas2Pad(100);

  // while (converter.ReadEvent() && converter.EventNumber<1000){
  while (converter.ReadEvent()) {

    if (targetEvent > 0 && converter.EventNumber != (uint32_t)targetEvent) {
      continue;
    }

    if (!targetEvents.empty()) {
      if (std::find(targetEvents.begin(), targetEvents.end(),
                    (uint32_t)targetEvent) != targetEvents.end()) {
        continue;
      }
    }

    ChannelDataInEvent.clear();

    while (converter.ReadADC()) {
      while (converter.ReadChannel()) {
        // Далее гистограмм для WF и подсчета интеграла
        int evNum = converter.channelEvent.eventNum;
        int chNum = converter.channelEvent.channelNums;
        int chPhi = converter.channelEvent.channel_Phi;
        int chZ = converter.channelEvent.channel_Z;
        std::string h_name =
            Form("h_ev_%i_ch_%i_phi_%i_z_%i", evNum, chNum, chPhi, chZ);
        std::string h_title =
            Form("Basket 38, Event %i, channel %i, Phi %i, Z %i", evNum, chNum,
                 chPhi, chZ);

        SetWaveformHistogram(converter.channelEvent.adcValues, h1_wf, h_name,
                             h_title);
        SetChannelIntegralAmp(converter.channelEvent, h1_wf, iBinStart,
                              iBinStop, iBinLeft, iBinRight);
        ChannelDataInEvent.push_back(converter.channelEvent);
        converter.FillTreeData();
      } // end ReadChannel
    } // end ReadADC

    WfDir->cd();

    if (converter.EventNumber > 1700 && converter.EventNumber < 2000)
      drawObj.UpdateAndSaveCanvas("", ChannelDataInEvent);

    if (targetEvent > 0 || targetEvents.empty() == false) {
      if (!targetEvents.empty()) {
        if (converter.EventNumber >= targetEvents.back()) {
          break;
        }
      }
      if (targetEvent > 0) {
        break;
      }
    }
  }

  // Запись последнего ивента со всеми значенями равными 0.
  converter.channelEvent.Clear();
  converter.FillTreeData();

  converter.outFile->cd();

  std::cout << "\nOutput Tree saved to: " << outputData << "\n" << std::endl;

  converter.WriteTreeAndClose();

  timer1.Stop();
  timer1.Print();
}