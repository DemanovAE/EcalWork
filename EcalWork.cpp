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

double ComputeStripContamination(const std::vector<ChannelData> &data,
                                 const std::vector<int> &stripIndices,
                                 double &sumSelected, double &sumNeighbors,
                                 int maxDeltaPhi = 1, int maxDeltaZ = 1) {
  sumSelected = 0.0;
  sumNeighbors = 0.0;

  if (data.empty() || stripIndices.empty())
    return 0.0;

  std::vector<bool> isStrip(data.size(), false);
  for (int idx : stripIndices) {
    if (idx >= 0 && idx < (int)data.size()) {
      isStrip[idx] = true;
      sumSelected += data[idx].integral;
    }
  }

  std::set<int> neighborIndices;

  for (int iStrip : stripIndices) {
    const auto &ci = data[iStrip];

    for (int j = 0; j < (int)data.size(); ++j) {
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

    // if (eventCh[i].integral < threshold_2)
    //     continue; // int < 500 is a basic selection  fro all channels in
    //     event
    hCut->Fill(3);

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
  hCut->Fill(4, data.size());

  std::vector<int> phiStrip, zStrip, stripIndices;

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
  hCut->Fill(5, data.size());

  //  User-configurable threshold_3 (e.g. 0.20 = 20%)
  const double threshold_3 = 0.20;

  /*
  // --- New 4.3: contamination check from neighboring channels ---
  double sumSelected = 0.0;
  double sumNeighbors = 0.0;
  double contamination = ComputeStripContamination(
      data, stripIndices, sumSelected, sumNeighbors, 1, 1);
  if (contamination > threshold_3) {
    // Reject this event as multi-muon (or heavily contaminated) in
  neighboring
    // cells
    eventCh.clear();
    return false;
  }
  */
  hCut->Fill(6, data.size());

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
  eventCh = std::move(filteredData);

  return true;
}

// Selection 3x3 Neighbourhood window
// Finding the second strongest cell
bool CheckThreeOnThreeWindow(std::vector<ChannelData> &data,
                             ChannelData hottestCell, float &adc_max_2) {
  adc_max_2 = -1;
  float adc_max = hottestCell.integral;
  int phi_0 = hottestCell.channel_Phi;
  int z_0 = hottestCell.channel_Z;

  float min_ratio = 0.70; // dummy value
  float max_ratio = 0.95; // dummy value

  for (int phi = phi_0 - 1; phi <= phi_0 + 1; ++phi) {
    for (int z = z_0 - 1; z <= z_0 + 1; ++z) {

      // Skip the center cell itself
      if (phi == phi_0 && z == z_0)
        continue;

      // Search for any channel with these indices in this event
      for (const auto &ch : data) {
        if (ch.channel_Phi == phi && ch.channel_Z == z) {
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

    // inside 5x5 window
    if (dPhi == 2 && dZ == 2) {
      sum5x5 += ch.integral;
    }
    // inside 3x3 window
    if (dPhi == 1 && dZ == 1) {
      sum3x3 += ch.integral;
    }
  }
  if (sum5x5 <= 0.0f)
    return false;

  sum3x3 = sum3x3 + hottestCell.integral;
  // energy fraction outside the compact 3x3 core
  // ratio_cut5x5 = (sum5x5 - sum3x3) / sum5x5;
  ratio_cut5x5 = sum5x5 / (sum3x3+sum5x5);
  //std::cout << "5x5: " << sum5x5 << " 3x3: " << sum3x3 << " hot: " << hottestCell.integral << " r: " << ratio_cut5x5 << std::endl;


  float max_diffusivity = 0.5f; // dummy value, tune from data
  if (ratio_cut5x5 > max_diffusivity)
    return false;

  return true;
}
bool LongAnalysis(std::vector<ChannelData> &eventCh,
                  std::vector<TH1F *> h_int_per_channel, TH1D *hCut) {
  std::vector<ChannelData> data;
  data.reserve(eventCh.size());

  for (int i = 0; i < eventCh.size(); i++) {
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
  auto itMax = std::max_element(data.begin(), data.end(),
                                [](const ChannelData &a, const ChannelData &b) {
                                  return a.integral < b.integral;
                                });

  if (itMax == eventCh.end())
    return false;

  const ChannelData &maxCell = *itMax;
  float adc_max = maxCell.integral;
  // float adc_max_2 = -1;
  float adc_max_2;

  // int maxThreshold = 1000; // dummy value
  // int secondNoise = 50;    // dummy value

  float ration3on3 = 0;

  // 2. Reject if hottest cell is too small
  // if (adc_max < maxThreshold)
  //     return false;
  hCut->Fill(4, data.size());

  if (!CheckThreeOnThreeWindow(data, maxCell, adc_max_2))
    return false;
  hCut->Fill(5, data.size());

  if (adc_max_2 < 500)
    return false;
  hCut->Fill(6, data.size());

  float sum5x5, sum3x3, ratio_cut5x5;

  if (!CalcAreaFiveOnFiveWindow(data, maxCell, sum5x5, sum3x3, ratio_cut5x5))
    return false;

  hCut->Fill(7, data.size());

  for (int i = 0; i < data.size(); i++) {
    int ch = data[i].channelNums;
    if (ch >= 1 && ch <= NCHANNELS) {
      int phi_0 = maxCell.channel_Phi;
      int z_0 = maxCell.channel_Z;
      if (abs(phi_0 - data[i].channel_Phi) <= 1 &&
          abs(z_0 - data[i].channel_Z) <= 1) {
        h_int_per_channel[ch - 1]->Fill(data[i].integral);
      }
    }
  }
  if (data[0].eventNum > 1 && data[0].eventNum < 100) {

    EcalDrawClass drawObj;
    drawObj.InitCanvas1Pad();
    drawObj.UpdateAndSaveCanvas("", data);
  }
  // std::cout << "Approx deposed energy = " << adc_max_2+adc_max  << "[ADC?]"
  // <<std::endl;
  eventCh.clear();
  eventCh = std::move(data);

  return true;
}

// void EcalWork(std::string inputDataTree = "out_all.root",
//               std::string outputData = "tran_x_38.root") {
void EcalWork(std::string inputDataTree = "out_all.root",
              std::string outputData = "trans_38.root",
              bool TransverAnalysis = true) {

  TStopwatch timer1;
  timer1.Start();

  // bool TransverAnalysis = true;
  // bool TransverAnalysis = false;

  std::vector<TH1F *> h_int_per_channel;
  h_int_per_channel.reserve(NCHANNELS);
  for (int i = 0; i < NCHANNELS; ++i) {
    h_int_per_channel.push_back(new TH1F(
        Form("h_int_ch_%d", i + 1),
        Form("Integral channel %d;Integral;Entries", i + 1), 100, 0, 15000));
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
  TH1D *h1_chNum = new TH1D("channel_number", ";Num;count", NCHANNELS, 0.5,
                            (Float_t)NCHANNELS + 0.5);
  TH1D *h1_chZ = new TH1D("channel_Z", ";Z;count", 64, 0.5, 64.5);
  TH1D *h1_chPhi = new TH1D("channel_Phi", ";Phi;count", 12, 0.5, 12.5);
  TH1D *h1_chAmp = new TH1D("channel_amplitude", ";Amp;count", 100, 0, 10000);
  TH1D *h1_chInt = new TH1D("channel_integral", ";Amp;count", 300, 0, 15000);

  std::filesystem::path outRootPath(outputData);
  std::filesystem::path outDir = outRootPath.parent_path();
  if (outDir.empty()) {
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
  std::cout << "\n\nOutput saved to: " << outputData << "\n\n" << std::endl;

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