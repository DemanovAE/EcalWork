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
#include "TDirectory.h"
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

// bool UseLongQA = true;
// bool UseTransQA = false;

bool SetBinNameCutHisto = true;

struct EcalConfig {
  bool transAnalysis;
  bool useLongQA;
  bool useTransQA;

  // Longitudinal cuts
  int long_max_1st_Integral;
  int long_max_2nd_Integral;
  float long_min_3x3_ratio;
  float long_max_3x3_ratio;
  float long_max_diffusivity_5x5;

  // Transverse cuts (steering)
  double trans_amp_thr1;    // threshold_1: basic amp cut (ADC)
  double trans_amp_thr2;    // threshold_2: strong amp cut (ADC)
  int trans_min_strip_len;  // minLen: minimal strip length (cells)
  double trans_contam_frac; // threshold_3: fraction (0–1), reserved
};

// For Trans Analysis
struct LongSelHists {
  TH1D *hCutFlow = nullptr;

  TH1F *hMaxInt_all = nullptr;
  TH1F *hMaxInt_pass = nullptr;

  TH1F *hSecondInt_all = nullptr;
  TH1F *hSecondInt_pass = nullptr;

  TH1F *hRatio3x3_all = nullptr;
  TH1F *hRatio3x3_pass = nullptr;

  TH1F *hDiff5x5_all = nullptr;
  TH1F *hDiff5x5_pass = nullptr;

  TH2F *hRatio3x3_vs_max_all = nullptr;
  TH2F *hDiff5x5_vs_core_all = nullptr;
};

struct TransSelHists {
  TH1D *hCutFlow = nullptr;

  TH1F *hNamp100_all = nullptr;
  TH1F *hNamp100_pass = nullptr;

  TH1F *hAmp_all = nullptr;
  TH1F *hAmp_pass = nullptr;

  TH1F *hStripLen_all = nullptr;
  TH1F *hStripLen_pass = nullptr;

  TH1F *hTargetIntegral_all = nullptr;
  TH1F *hTargetIntegral_pass = nullptr;

  TH1F *hContamCount_all = nullptr;
  TH1F *hContamCount_pass = nullptr;
};

inline void FillIf(TH1 *h, double x) {
  if (h)
    h->Fill(x);
}

inline void FillIf(TH2 *h, double x, double y) {
  if (h)
    h->Fill(x, y);
}

LongSelHists BookLongSelHists(TDirectory *dir) {
  LongSelHists h;
  if (!dir)
    return h;
  dir->cd();

  h.hCutFlow = new TH1D(
      "hLongCutFlow", "Longitudinal selection cut flow;;Events", 10, -0.5, 9.5);
  h.hCutFlow->GetXaxis()->SetBinLabel(1, "entered");
  h.hCutFlow->GetXaxis()->SetBinLabel(2, "non-empty");
  h.hCutFlow->GetXaxis()->SetBinLabel(3, "1st max");
  h.hCutFlow->GetXaxis()->SetBinLabel(4, "phi edge");
  h.hCutFlow->GetXaxis()->SetBinLabel(5, "z edge");
  h.hCutFlow->GetXaxis()->SetBinLabel(6, "3x3 ratio");
  h.hCutFlow->GetXaxis()->SetBinLabel(7, "2nd max");
  h.hCutFlow->GetXaxis()->SetBinLabel(8, "5x5 diff");

  h.hMaxInt_all =
      new TH1F("hMaxInt_all", "1st max integral before cut;Integral;Events",
               200, 0, 20000);
  h.hMaxInt_pass =
      new TH1F("hMaxInt_pass", "1st max integral after cut;Integral;Events",
               200, 0, 20000);

  h.hSecondInt_all =
      new TH1F("hSecondInt_all", "2nd max in 3x3 before cut;Integral;Events",
               200, 0, 5000);
  h.hSecondInt_pass =
      new TH1F("hSecondInt_pass", "2nd max in 3x3 after cut;Integral;Events",
               200, 0, 5000);

  h.hRatio3x3_all = new TH1F(
      "hRatio3x3_all", "max/(max+2nd) before cut;Ratio;Events", 100, 0.0, 1.0);
  h.hRatio3x3_pass = new TH1F(
      "hRatio3x3_pass", "max/(max+2nd) after cut;Ratio;Events", 100, 0.0, 1.0);

  h.hDiff5x5_all = new TH1F(
      "hDiff5x5_all", "(5x5-3x3)/(5x5) before cut;Ratio;Events", 100, 0.0, 1.0);
  h.hDiff5x5_pass = new TH1F(
      "hDiff5x5_pass", "(5x5-3x3)/(5x5) after cut;Ratio;Events", 100, 0.0, 1.0);

  h.hRatio3x3_vs_max_all = new TH2F(
      "hRatio3x3_vs_max_all", "ratio3x3 vs 1st max;1st max integral;ratio3x3",
      120, 0, 20000, 100, 0.0, 1.0);

  h.hDiff5x5_vs_core_all =
      new TH2F("hDiff5x5_vs_core_all",
               "diffusivity vs core energy;core energy;diffusivity", 120, 0,
               30000, 100, 0.0, 1.0);

  return h;
}

TransSelHists BookTransSelHists(TDirectory *dir) {
  TransSelHists h;
  if (!dir)
    return h;
  dir->cd();

  h.hCutFlow = new TH1D("hTransCutFlow",
                        "Transverse selection cut flow;;Events", 10, -0.5, 9.5);
  h.hCutFlow->GetXaxis()->SetBinLabel(1, "entered");
  h.hCutFlow->GetXaxis()->SetBinLabel(2, "amp>thr");
  h.hCutFlow->GetXaxis()->SetBinLabel(3, "N>=5");
  h.hCutFlow->GetXaxis()->SetBinLabel(4, "strip found");
  h.hCutFlow->GetXaxis()->SetBinLabel(5, "contam");
  h.hCutFlow->GetXaxis()->SetBinLabel(6, "target found");

  h.hNamp100_all = new TH1F(
      "hNamp100_all", "N channels with amp>threshold before cuts;N;Events", 70,
      -0.5, 69.5);
  h.hNamp100_pass = new TH1F(
      "hNamp100_pass", "N channels with amp>threshold after N>=5;N;Events", 70,
      -0.5, 69.5);

  h.hAmp_all = new TH1F(
      "hAmp_all", "Channel amplitude entering transverse;Amplitude;Counts", 200,
      0, 5000);
  h.hAmp_pass =
      new TH1F("hAmp_pass", "Channel amplitude after amp cut;Amplitude;Counts",
               200, 0, 5000);

  h.hStripLen_all = new TH1F(
      "hStripLen_all", "Found strip length before strip cut;Length;Events", 20,
      -0.5, 19.5);
  h.hStripLen_pass =
      new TH1F("hStripLen_pass", "Found strip length accepted;Length;Events",
               20, -0.5, 19.5);

  h.hTargetIntegral_all = new TH1F(
      "hTargetIntegral_all",
      "Target integral before final accept;Integral;Events", 200, 0, 15000);
  h.hTargetIntegral_pass =
      new TH1F("hTargetIntegral_pass",
               "Target integral accepted;Integral;Events", 200, 0, 15000);

  h.hContamCount_all = new TH1F(
      "hContamCount_all", "Nearby contamination count before cut;N;Events", 20,
      -0.5, 19.5);
  h.hContamCount_pass =
      new TH1F("hContamCount_pass",
               "Nearby contamination count accepted;N;Events", 20, -0.5, 19.5);

  return h;
}
void WriteLongSelHists(const LongSelHists *h) {
  if (!h)
    return;
  if (h->hCutFlow)
    h->hCutFlow->Write();
  if (h->hMaxInt_all)
    h->hMaxInt_all->Write();
  if (h->hMaxInt_pass)
    h->hMaxInt_pass->Write();
  if (h->hSecondInt_all)
    h->hSecondInt_all->Write();
  if (h->hSecondInt_pass)
    h->hSecondInt_pass->Write();
  if (h->hRatio3x3_all)
    h->hRatio3x3_all->Write();
  if (h->hRatio3x3_pass)
    h->hRatio3x3_pass->Write();
  if (h->hDiff5x5_all)
    h->hDiff5x5_all->Write();
  if (h->hDiff5x5_pass)
    h->hDiff5x5_pass->Write();
  if (h->hRatio3x3_vs_max_all)
    h->hRatio3x3_vs_max_all->Write();
  if (h->hDiff5x5_vs_core_all)
    h->hDiff5x5_vs_core_all->Write();
}

void WriteTransSelHists(const TransSelHists *h) {
  if (!h)
    return;
  if (h->hCutFlow)
    h->hCutFlow->Write();
  if (h->hNamp100_all)
    h->hNamp100_all->Write();
  if (h->hNamp100_pass)
    h->hNamp100_pass->Write();
  if (h->hAmp_all)
    h->hAmp_all->Write();
  if (h->hAmp_pass)
    h->hAmp_pass->Write();
  if (h->hStripLen_all)
    h->hStripLen_all->Write();
  if (h->hStripLen_pass)
    h->hStripLen_pass->Write();
  if (h->hTargetIntegral_all)
    h->hTargetIntegral_all->Write();
  if (h->hTargetIntegral_pass)
    h->hTargetIntegral_pass->Write();
  if (h->hContamCount_all)
    h->hContamCount_all->Write();
  if (h->hContamCount_pass)
    h->hContamCount_pass->Write();
}

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
                  bool alongPhi, int minLen = 5) {

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
// Example: stricter contamination check around strip edges.
// This is NOT used in current production; call is commented out in
// TransverseAnalysis. The idea: count non-strip hits in a band of width 1
// around the strip, especially near the ends, and reject if too many.

bool CheckEdgeContamination(const std::vector<ChannelData> &eventHits,
                            const std::vector<ChannelData> &stripHits,
                            bool alongPhi, double noiseIntegralThreshold,
                            int maxEdgeHits,
                            int &edgeCount, // ← added
                            int bandRadiusZ = 1, int bandRadiusPhi = 1) {
  edgeCount = 0;

  std::set<std::pair<int, int>> stripCoords;
  for (const auto &ch : stripHits)
    stripCoords.insert({(int)ch.channel_Z, (int)ch.channel_Phi});

  int minInner = +9999, maxInner = -9999;
  for (const auto &ch : stripHits) {
    int inner = alongPhi ? (int)ch.channel_Z : (int)ch.channel_Phi;
    if (inner < minInner)
      minInner = inner;
    if (inner > maxInner)
      maxInner = inner;
  }

  for (const auto &hit : eventHits) {
    if (hit.integral < noiseIntegralThreshold)
      continue;
    int z = (int)hit.channel_Z;
    int phi = (int)hit.channel_Phi;
    if (stripCoords.count({z, phi}))
      continue;

    int dz = 0, dphi = 0;
    if (alongPhi) {
      // Phi-strip: fixed Phi row, extended in Z
      // dz   = distance beyond strip Z endpoints
      // dphi = distance from strip's Phi row (sideways)
      int stripPhi = (int)stripHits.front().channel_Phi;
      dz = z - std::max(std::min(z, maxInner), minInner);
      dphi = phi - stripPhi;
    } else {
      // Z-strip: fixed Z column, extended in Phi
      // dphi = distance beyond strip Phi endpoints
      // dz   = distance from strip's Z column (sideways)
      int stripZ = (int)stripHits.front().channel_Z;
      dz = z - stripZ;
      dphi = phi - std::max(std::min(phi, maxInner), minInner);
    }

    if (std::abs(dz) <= bandRadiusZ && std::abs(dphi) <= bandRadiusPhi) {
      ++edgeCount;
      if (edgeCount > maxEdgeHits)
        return false;
    }
  }
  return true;
}

// Find a target channel in a strip: center of 5 consecutive channels
bool FindTargetChannelInStrip(const std::vector<ChannelData> &strip,
                              bool alongPhi, double threshold_2,
                              ChannelData &target, const EcalConfig &cfg) {

  // std::cout << "[FindTarget] strip.size=" << strip.size()
  //           << " alongPhi=" << alongPhi << " thr2=" << threshold_2
  //           << "   cfg.thr2 " << cfg.trans_amp_thr2 << "\n";

  // for (const auto &ch : strip) {
  //   std::cout << "    [strip cell] Phi=" << ch.channel_Phi
  //             << " Z=" << ch.channel_Z << " int=" << ch.integral << "\n";
  // }

  // Need at least 5 channels to form a consecutive group
  if (strip.size() < cfg.trans_min_strip_len)
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

    // std::cout << "  center i=" << i << " coord="
    //           << (alongPhi ? (int)ordered[i].channel_Z
    //                        : (int)ordered[i].channel_Phi)
    //           << " consecutive=" << consecutive << "\n";

    // If not consecutive, skip this center
    if (!consecutive)
      continue;

    // std::cout << "    neighbors integrals=" << ordered[i - 2].integral << ","
    //           << ordered[i - 1].integral << "," << ordered[i + 1].integral
    //           << "," << ordered[i + 2].integral << "\n";

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
  std::cout << "[FindTarget] no suitable center found\n";

  // No suitable target found
  return false;
}

bool TransverseAnalysis(std::vector<ChannelData> &eventCh,
                        std::vector<TH1F *> h_int_per_channel, TH1D *hCut,
                        const EcalConfig &cfg, TransSelHists *hs = nullptr) {
  std::vector<int> _PhiCount(12);
  std::vector<int> _ZCount(64);
  std::vector<ChannelData> data;
  data.reserve(eventCh.size());
  // FillIf(hs ? hs->hCutFlow : nullptr, 0);

  const double threshold_1 = cfg.trans_amp_thr1;
  const double threshold_2 = cfg.trans_amp_thr2;

  for (int i = 0; i < eventCh.size(); i++) {
    FillIf(hs ? hs->hAmp_all : nullptr, eventCh[i].amplitude);
    if (eventCh[i].amplitude < threshold_1)
      continue; // apm < 100 is a basic selection fro all channels in event
    hCut->Fill(2);
    FillIf(hs ? hs->hAmp_pass : nullptr, eventCh[i].amplitude);

    data.push_back(eventCh[i]); // we have got the 'target' channels in event
    _PhiCount[eventCh[i].channel_Phi -
              1]++; // count the number of phi triggered cells in event
    _ZCount[eventCh[i].channel_Z -
            1]++; // count the number of Z triggered cells in event
  }
  FillIf(hs ? hs->hNamp100_all : nullptr, data.size());
  if (data.size() < 5) { // if there are less tna five channels (doesn't
                         // matter which ones), then reject the event
    eventCh.clear();
    return false;
  }
  FillIf(hs ? hs->hNamp100_pass : nullptr, data.size());
  FillIf(hs ? hs->hCutFlow : nullptr, 1);
  hCut->Fill(3, data.size());

  std::vector<int> phiStrip, zStrip;

  int phiStripOuter = -1, zStripOuter = -1;
  bool hasPhiStrip = FindOneStrip(data, phiStrip, phiStripOuter, true,
                                  (int)cfg.trans_min_strip_len);
  bool hasZStrip = FindOneStrip(data, zStrip, zStripOuter, false,
                                (int)cfg.trans_min_strip_len);

  FillIf(hs ? hs->hStripLen_all : nullptr,
         hasPhiStrip ? (int)phiStrip.size()
                     : (hasZStrip ? (int)zStrip.size() : 0));

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
  FillIf(hs ? hs->hStripLen_pass : nullptr, filteredData.size());
  FillIf(hs ? hs->hCutFlow : nullptr, 2);

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
  const double threshold_3 = cfg.trans_contam_frac; // reserved, not used yet

  // FillIf(hs ? hs->hContamCount_pass : nullptr, contaminationCount);
  FillIf(hs ? hs->hCutFlow : nullptr, 3);

  hCut->Fill(5, data.size());

  // Edge contamination: reject too-diagonal tracks
  // noiseThr: only count hits above this integral as contamination
  // maxDiagHits=99 for QA-only pass, lower later based on hContamCount_all
  double noiseThr = threshold_2 * 0.3;
  int maxDiagHits = 99; // set to 99 first to just fill QA, then tune
  int edgeCount = 0;

  if (hasPhiStrip) {
    bool contamOk = CheckEdgeContamination(data, filteredData,
                                           /*alongPhi=*/true, noiseThr,
                                           maxDiagHits, edgeCount,
                                           /*bandRadiusZ=*/1,
                                           /*bandRadiusPhi=*/1);
    FillIf(hs ? hs->hContamCount_all : nullptr, edgeCount);
    if (!contamOk) {
      eventCh.clear();
      return false;
    }
    FillIf(hs ? hs->hContamCount_pass : nullptr, edgeCount);

  } else if (hasZStrip) {
    bool contamOk = CheckEdgeContamination(data, filteredData,
                                           /*alongPhi=*/false, noiseThr,
                                           maxDiagHits, edgeCount,
                                           /*bandRadiusZ=*/1,
                                           /*bandRadiusPhi=*/1);
    FillIf(hs ? hs->hContamCount_all : nullptr, edgeCount);
    if (!contamOk) {
      eventCh.clear();
      return false;
    }
    FillIf(hs ? hs->hContamCount_pass : nullptr, edgeCount);
  }

  FillIf(hs ? hs->hCutFlow : nullptr, 3);
  hCut->Fill(5, data.size());

  // for (int i = 0; i < data.size(); i++) {
  //   int ch = data[i].channelNums;
  //   if (ch >= 1 && ch <= NCHANNELS) {
  //     h_int_per_channel[ch - 1]->Fill(data[i].integral);
  //   }
  // }

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
                                         targetChannel, cfg);
  } else if (hasZStrip) {
    hasTarget = FindTargetChannelInStrip(filteredData, false, threshold_2,
                                         targetChannel, cfg);
  }
  if (hasTarget)
    FillIf(hs ? hs->hTargetIntegral_all : nullptr, targetChannel.integral);

  if (!hasTarget) {
    eventCh.clear();
    return false;
  }

  if (hasTarget) {
    int ch = targetChannel.channelNums;
    if (ch >= 1 && ch <= NCHANNELS)
      h_int_per_channel[ch - 1]->Fill(targetChannel.integral);
  }
  
  FillIf(hs ? hs->hTargetIntegral_pass : nullptr, targetChannel.integral);
  FillIf(hs ? hs->hCutFlow : nullptr, 4);

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
                  std::vector<TH1F *> h_int_per_channel, const EcalConfig &cfg,
                  TH1D *hCut, LongSelHists *hs = nullptr) {

  std::vector<ChannelData> data;
  data.reserve(eventCh.size());
  int itMax2 = -1;
  FillIf(hs ? hs->hCutFlow : nullptr, 0);

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
  FillIf(hs ? hs->hCutFlow : nullptr, 1);

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

  FillIf(hs ? hs->hMaxInt_all : nullptr, adc_max);

  if (adc_max < cfg.long_max_1st_Integral)
    return false;

  FillIf(hs ? hs->hMaxInt_pass : nullptr, adc_max);
  FillIf(hs ? hs->hCutFlow : nullptr, 2);

  hCut->Fill(4, data.size());
  if (SetBinNameCutHisto)
    hCut->GetXaxis()->SetBinLabel(
        5, Form("1stMaxInt>%i", cfg.long_max_1st_Integral));

  if (maxCell.channel_Phi == 1 || maxCell.channel_Phi == 12)
    return false;

  FillIf(hs ? hs->hCutFlow : nullptr, 3);

  hCut->Fill(5, data.size());
  if (SetBinNameCutHisto)
    hCut->GetXaxis()->SetBinLabel(6, "Phi==1;12");

  if (maxCell.channel_Z == 1 || maxCell.channel_Z == 64)
    return false;

  FillIf(hs ? hs->hCutFlow : nullptr, 4);

  hCut->Fill(6, data.size());
  if (SetBinNameCutHisto)
    hCut->GetXaxis()->SetBinLabel(7, "Z==1;64");

  float ratio_3x3 = 0;
  ChannelData secMaxCell;
  if (!CheckThreeOnThreeWindow(data, maxCell, secMaxCell, itMax2, ratio_3x3))
    return false;

  FillIf(hs ? hs->hSecondInt_all : nullptr, secMaxCell.integral);
  FillIf(hs ? hs->hRatio3x3_all : nullptr, ratio_3x3);
  FillIf(hs ? hs->hRatio3x3_vs_max_all : nullptr, adc_max, ratio_3x3);

  if (ratio_3x3 >= cfg.long_max_3x3_ratio || ratio_3x3 < cfg.long_min_3x3_ratio)
    return false;

  FillIf(hs ? hs->hRatio3x3_pass : nullptr, ratio_3x3);
  FillIf(hs ? hs->hCutFlow : nullptr, 5);

  hCut->Fill(7, data.size());
  if (SetBinNameCutHisto)
    hCut->GetXaxis()->SetBinLabel(
        8, Form("1st/(1st+2nd)>%.2f", cfg.long_min_3x3_ratio));

  float MaxIntegral_2st = secMaxCell.integral;
  if (MaxIntegral_2st < cfg.long_max_2nd_Integral)
    return false;

  FillIf(hs ? hs->hSecondInt_pass : nullptr, MaxIntegral_2st);
  FillIf(hs ? hs->hCutFlow : nullptr, 6);

  hCut->Fill(8, data.size());
  if (SetBinNameCutHisto)
    hCut->GetXaxis()->SetBinLabel(
        9, Form("2ndMaxInt>%i", cfg.long_max_2nd_Integral));

  float sum5x5, sum3x3, ratio_cut5x5;

  if (!CalcAreaFiveOnFiveWindow(data, maxCell, sum5x5, sum3x3, ratio_cut5x5))
    return false;

  const float coreEnergy = maxCell.integral + secMaxCell.integral;

  FillIf(hs ? hs->hDiff5x5_all : nullptr, ratio_cut5x5);
  FillIf(hs ? hs->hDiff5x5_vs_core_all : nullptr, coreEnergy, ratio_cut5x5);

  if (ratio_cut5x5 > cfg.long_max_diffusivity_5x5)
    return false;
  hCut->Fill(9, data.size());
  if (SetBinNameCutHisto) {
    hCut->GetXaxis()->SetBinLabel(
        10, Form("(5x5)/(3x3+5x5)<%.2f", cfg.long_max_diffusivity_5x5));
    SetBinNameCutHisto = false;
  }

  FillIf(hs ? hs->hDiff5x5_pass : nullptr, ratio_cut5x5);
  FillIf(hs ? hs->hCutFlow : nullptr, 7);

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

  // Store in histogram of the hottest channel
  const int ch_core = maxCell.channelNums;
  if (ch_core >= 1 && ch_core <= NCHANNELS && h_int_per_channel[ch_core - 1]) {
    h_int_per_channel[ch_core - 1]->Fill(coreEnergy);
  }
  std::vector<ChannelData> FinalCh = {maxCell, secMaxCell};

  // if (data[0].eventNum > 10000 && data[0].eventNum < 10100) {
  //   EcalDrawClass drawObj;
  //   drawObj.InitCanvas1Pad();
  //   drawObj.UpdateAndSaveCanvas("", data);
  //   drawObj.UpdateAndSaveCanvas("_cut", FinalCh);
  // }
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

// void EcalWork(
//     std::string inputDataTree =
//         "/nica/mpd1/demanov/ecal_mpd/run_rc1-hs4_133_basket5_1616162.root",
//     std::string outputData = "hs4_133_basket5_1616162.root",
//     Long64_t firstEntry = 0, Long64_t lastEntry = 20.e6, const EcalConfig
//     &cfg) {

void EcalWork(std::string inputDataTree, std::string outputData,
              Long64_t firstEntry, Long64_t lastEntry, const EcalConfig &cfg) {

  TStopwatch timer1;
  timer1.Start();

  // Per-channel histograms for target-cell integral / core energy
  std::vector<TH1F *> h_int_per_channel;
  h_int_per_channel.reserve(NCHANNELS);
  // for (int i = 0; i < NCHANNELS; ++i) {
  //   h_int_per_channel.push_back(
  //       new TH1F(Form("h_int_ch_%d", i + 1),
  //                Form("Integral channel %d;Integral;Entries", i + 1), 100, 0,
  //                cfg.transAnalysis == true ? 15000 : 1.e5));
  // }
  for (int i = 0; i < NCHANNELS; ++i) {
    const char *hname = Form("h_int_ch_%d", i + 1);

    const char *htitle =
        cfg.transAnalysis
            ? Form("Target integral channel %d;Integral;Entries", i + 1)
            : Form("Core energy (1st+2nd) channel %d;Core energy;Entries",
                   i + 1);

    h_int_per_channel.push_back(
        new TH1F(hname, htitle, 100, 0, cfg.transAnalysis ? 15000 : 1.e5));
  }

  // Use outputData exactly as given (full path from script)
  TFile *outFile = new TFile(outputData.c_str(), "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    std::cerr << "Cannot create output file " << outputData << std::endl;
    return;
  }

  // Input file and tree
  TFile *iFile = TFile::Open(inputDataTree.c_str());
  if (!iFile || iFile->IsZombie()) {
    std::cerr << "Cannot open input file " << inputDataTree << std::endl;
    return;
  }

  TTree *tree = dynamic_cast<TTree *>(iFile->Get("events"));
  if (!tree) {
    std::cerr << "Cannot find TTree 'events' in file " << inputDataTree
              << std::endl;
    return;
  }

  Long64_t nEntries = tree->GetEntries();
  std::cout << "The number of entries in decoded file '"
            << inputDataTree.c_str() << "' " << nEntries << std::endl;

  // Clamp entry range
  if (lastEntry < 0 || lastEntry > nEntries) {
    lastEntry = nEntries;
  }
  if (firstEntry < 0) {
    firstEntry = 0;
  }
  if (firstEntry >= lastEntry) {
    std::cerr << "Empty entry range: firstEntry=" << firstEntry
              << " lastEntry=" << lastEntry << std::endl;
    return;
  }

  std::cout << "Processing entries in range [" << firstEntry << ", "
            << lastEntry << ")" << std::endl;

  // Branch variables
  UInt_t b_eventNum = 0;
  Int_t b_integral = 0;
  Int_t b_amplitude = 0;
  unsigned short b_chNum = 0;
  unsigned short b_chZ = 0;
  unsigned short b_chPhi = 0;

  tree->SetBranchAddress("eventNum", &b_eventNum);
  tree->SetBranchAddress("integral", &b_integral);
  tree->SetBranchAddress("amplitude", &b_amplitude);
  tree->SetBranchAddress("channelNums", &b_chNum);
  tree->SetBranchAddress("channel_Z", &b_chZ);
  tree->SetBranchAddress("channel_Phi", &b_chPhi);

  // Summary histograms
  TH1D *h1_evNum = new TH1D("event_number", ";event;count", 2.0e6, 0, 2.0e6);
  TH1D *h1_chNum = new TH1D("channel_number", ";Num;count", NCHANNELS, 0.5,
                            (Float_t)NCHANNELS + 0.5);
  TH1D *h1_chZ = new TH1D("channel_Z", ";Z;count", 64, 0.5, 64.5);
  TH1D *h1_chPhi = new TH1D("channel_Phi", ";Phi;count", 12, 0.5, 12.5);
  TH1D *h1_chAmp = new TH1D("channel_amplitude", ";Amp;count", 100, 0,
                            cfg.transAnalysis == true ? 10000 : 30000);
  TH1D *h1_chInt = new TH1D("channel_integral", ";Amp;count", 500, 0,
                            cfg.transAnalysis == true ? 15000 : 1.e5);

  // Pictures directory
  std::filesystem::path outRootPath(outputData);
  std::filesystem::path outDir = outRootPath.parent_path();
  if (outDir.empty()) {
    outDir = ".";
  }

  std::filesystem::path pictDir =
      outDir / "pict" / (cfg.transAnalysis ? "transverse" : "longitudinal");
  std::filesystem::create_directories(pictDir);
  gSystem->ChangeDirectory(pictDir.string().c_str());
  std::cout << "Pictures will be saved under: " << pictDir << std::endl;

  // WF directory (kept for consistency, though not used now)
  TDirectory *WfDir = outFile->mkdir("channel_wf");
  WfDir->cd();

  LongSelHists longHsLocal;
  TransSelHists transHsLocal;
  LongSelHists *pLongHs = nullptr;
  TransSelHists *pTransHs = nullptr;

  outFile->cd();
  if (!cfg.transAnalysis && cfg.useLongQA) {
    TDirectory *dirLongQA = outFile->mkdir("longitudinal_QA");
    longHsLocal = BookLongSelHists(dirLongQA);
    pLongHs = &longHsLocal;
  }

  if (cfg.transAnalysis && cfg.useTransQA) {
    TDirectory *dirTransQA = outFile->mkdir("transverse_QA");
    transHsLocal = BookTransSelHists(dirTransQA);
    pTransHs = &transHsLocal;
  }

  // Per-event handling
  UInt_t currentEventNum = 0;
  bool hasData = false;
  std::vector<ChannelData> ChannelDataInEvent;
  ChannelDataInEvent.reserve(NCHANNELS);

  TH1D *h_CountCut =
      new TH1D("hCountCut", "Number of events after applying cut;;count", 20,
               -0.5, 19.5);

  // Main loop over the selected entry range
  for (Long64_t iEntry = firstEntry; iEntry < lastEntry; ++iEntry) {
    tree->GetEntry(iEntry);

    ChannelData data;
    data.eventNum = b_eventNum;
    data.integral = b_integral;
    data.amplitude = b_amplitude;
    data.channelNums = b_chNum;
    data.channel_Z = b_chZ;
    data.channel_Phi = b_chPhi;

    h_CountCut->Fill(0.);
    h_CountCut->GetXaxis()->SetBinLabel(1, "def");

    // Event building based on eventNum, same logic as original code
    if (!hasData) {
      currentEventNum = data.eventNum;
      hasData = true;
    } else if (data.eventNum != currentEventNum) {
      // Process finished event
      if (ChannelDataInEvent.size() <= MAX_CH) {
        h_CountCut->Fill(1, ChannelDataInEvent.size());
        h_CountCut->GetXaxis()->SetBinLabel(2, "Ch.Size()<65");

        bool KeyChData = false;

        if (cfg.transAnalysis) {
          KeyChData = TransverseAnalysis(ChannelDataInEvent, h_int_per_channel,
                                         h_CountCut, cfg, pTransHs);
        } else {
          KeyChData = LongAnalysis(ChannelDataInEvent, h_int_per_channel, cfg,
                                   h_CountCut, pLongHs);
        }

        if (KeyChData) {
          for (int i = 0; i < (int)ChannelDataInEvent.size(); i++) {
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
    }

    ChannelDataInEvent.push_back(data);
  }

  // Process the last event in the range (if any)
  if (!ChannelDataInEvent.empty() && ChannelDataInEvent.size() <= MAX_CH) {
    bool KeyChData = false;

    if (cfg.transAnalysis) {
      KeyChData = TransverseAnalysis(ChannelDataInEvent, h_int_per_channel,
                                     h_CountCut, cfg, pTransHs);
    } else {
      KeyChData = KeyChData = LongAnalysis(
          ChannelDataInEvent, h_int_per_channel, cfg, h_CountCut, pLongHs);
    }

    if (KeyChData) {
      for (int i = 0; i < (int)ChannelDataInEvent.size(); i++) {
        h1_evNum->Fill(ChannelDataInEvent[i].eventNum);
        h1_chAmp->Fill(ChannelDataInEvent[i].amplitude);
        h1_chInt->Fill(ChannelDataInEvent[i].integral);
        h1_chNum->Fill(ChannelDataInEvent[i].channelNums);
        h1_chPhi->Fill(ChannelDataInEvent[i].channel_Phi);
        h1_chZ->Fill(ChannelDataInEvent[i].channel_Z);
      }
    }
  }

  // Write basic histograms
  outFile->cd();

  if (pLongHs) {
    TDirectory *dirLongQA = (TDirectory *)outFile->Get("longitudinal_QA");
    if (dirLongQA) {
      dirLongQA->cd();
      WriteLongSelHists(pLongHs);
    }
    outFile->cd();
  }

  if (pTransHs) {
    TDirectory *dirTransQA = (TDirectory *)outFile->Get("transverse_QA");
    if (dirTransQA) {
      dirTransQA->cd();
      WriteTransSelHists(pTransHs);
    }
    outFile->cd();
  }

  h1_evNum->Write();
  h_CountCut->Write();
  h1_chAmp->Write();
  h1_chInt->Write();
  h1_chNum->Write();
  h1_chPhi->Write();
  h1_chZ->Write();

  // Subdirectory for target-channel histograms
  TDirectory *dirTargets = (TDirectory *)outFile->Get("target_channels");
  if (!dirTargets) {
    dirTargets = outFile->mkdir("target_channels");
  }
  dirTargets->cd();

  // Write per-channel target histograms
  for (auto *h : h_int_per_channel) {
    if (!h)
      continue;
    h->Write();
  }

  outFile->cd();

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