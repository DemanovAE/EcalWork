#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TObject.h"
#include "TString.h"
#include "TTree.h"

const int NCHANNELS = 768;
bool TransverAnalysis = true;

// Find the bin with maximum content in [xmin, xmax]
int FindMaxBinInRange(TH1 *h, double xmin, double xmax) {
  if (!h)
    return -1;

  int bmin = h->GetXaxis()->FindBin(xmin);
  int bmax = h->GetXaxis()->FindBin(xmax);

  bmin = std::max(1, bmin);
  bmax = std::min(h->GetNbinsX(), bmax);
  if (bmin > bmax)
    return -1;

  int    bestBin = -1;
  double bestY   = -1.0;

  for (int b = bmin; b <= bmax; ++b) {
    double y = h->GetBinContent(b);
    if (y > bestY) {
      bestY  = y;
      bestBin = b;
    }
  }

  return bestBin;
}

// ----------------------------------------------------------------------
// Transverse: fit Landau to per-channel target histos
// ----------------------------------------------------------------------
void FitTransverseTargetChannelHistograms(
    TFile *outFile,
    const std::string &inputDirName  = "target_channels",
    const std::string &outputDirName = "fitted_target_channels_transverse",
    int  minEntriesToFit             = 50,
    int  rebinFactor                 = 2)
{
  if (!outFile || outFile->IsZombie()) {
    std::cout
        << "FitTransverseTargetChannelHistograms: output file is null or zombie"
        << std::endl;
    return;
  }

  TDirectory *dirIn = (TDirectory *)outFile->Get(inputDirName.c_str());
  if (!dirIn) {
    std::cout << "FitTransverseTargetChannelHistograms: directory "
              << inputDirName << " not found" << std::endl;
    return;
  }

  outFile->cd();
  TDirectory *dirOut = (TDirectory *)outFile->Get(outputDirName.c_str());
  if (!dirOut)
    dirOut = outFile->mkdir(outputDirName.c_str());

  if (!dirOut) {
    std::cout << "FitTransverseTargetChannelHistograms: cannot create output "
                 "directory "
              << outputDirName << std::endl;
    return;
  }
  dirOut->cd();
  dirOut->Delete("*;*"); // remove all keys in this directory

  std::vector<double> v_ch, v_mpv, v_mpv_err, v_width, v_width_err, v_entries;
  int nFitted = 0;

  // Optional: per-channel tree for calibration
  TTree *tCal = new TTree("trans_calib", "Transverse per-channel calibration");
  double t_ch = 0, t_mpv = 0, t_mpvErr = 0, t_width = 0, t_widthErr = 0, t_entries = 0;
  tCal->Branch("channel", &t_ch, "channel/D");
  tCal->Branch("mpv",     &t_mpv,     "mpv/D");
  tCal->Branch("mpvErr",  &t_mpvErr,  "mpvErr/D");
  tCal->Branch("width",   &t_width,   "width/D");
  tCal->Branch("widthErr",&t_widthErr,"widthErr/D");
  tCal->Branch("entries", &t_entries, "entries/D");

  for (int ch = 1; ch <= NCHANNELS; ++ch) {
    TH1F *hIn = (TH1F *)dirIn->Get(Form("h_int_ch_%d", ch));
    if (!hIn)
      continue;
    if (hIn->GetEntries() < minEntriesToFit)
      continue;

    TH1F *hFit = (TH1F *)hIn->Clone(Form("h_fit_trans_ch_%d", ch));
    hFit->SetDirectory(nullptr);

    if (rebinFactor > 1)
      hFit->Rebin(rebinFactor);

    int nBins = hFit->GetNbinsX();
    if (nBins < 10) {
      delete hFit;
      continue;
    }

    int maxBin   = hFit->GetMaximumBin();
    double peakX = hFit->GetBinCenter(maxBin);
    double peakY = hFit->GetBinContent(maxBin);

    if (peakY <= 0) {
      delete hFit;
      continue;
    }

    double frac = 0.35 * peakY;
    int leftBin  = maxBin;
    int rightBin = maxBin;

    while (leftBin > 1 && hFit->GetBinContent(leftBin) > frac)
      --leftBin;
    while (rightBin < nBins && hFit->GetBinContent(rightBin) > frac)
      ++rightBin;

    double fitMin = hFit->GetBinCenter(leftBin);
    double fitMax = hFit->GetBinCenter(rightBin);

    double widthGuess = peakX - fitMin;
    if (widthGuess <= 0)
      widthGuess = hFit->GetRMS() / 4.0;
    if (widthGuess <= 0)
      widthGuess = 500.0;

    fitMin = std::max(0.0, peakX - 1.5 * widthGuess);
    fitMax = peakX + 4.0 * widthGuess;

    if (peakX > 1500.0 && fitMin < 0.5 * peakX)
      fitMin = 0.5 * peakX;

    if (fitMax <= fitMin) {
      delete hFit;
      continue;
    }

    TF1 *fLan =
        new TF1(Form("f_landau_trans_ch_%d", ch), "landau", fitMin, fitMax);
    fLan->SetParameters(peakY, peakX, widthGuess);
    fLan->SetParLimits(1, 4000.0, 9000.0);
    fLan->SetLineColor(kRed);
    fLan->SetLineWidth(2);

    int fitStatus = hFit->Fit(fLan, "QRS");
    if (fitStatus != 0) {
      delete hFit;
      delete fLan;
      continue;
    }

    double mpv    = fLan->GetParameter(1);
    double mpvErr = fLan->GetParError(1);
    double width  = fLan->GetParameter(2);
    double widthErr = fLan->GetParError(2);
    double entries  = hFit->GetEntries();

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
    hFit->Write("", TObject::kOverwrite);
    fLan->Write("", TObject::kOverwrite);

    v_ch.push_back(ch);
    v_mpv.push_back(mpv);
    v_mpv_err.push_back(mpvErr);
    v_width.push_back(width);
    v_width_err.push_back(widthErr);
    v_entries.push_back(entries);
    ++nFitted;

    // Fill tree
    t_ch      = ch;
    t_mpv     = mpv;
    t_mpvErr  = mpvErr;
    t_width   = width;
    t_widthErr= widthErr;
    t_entries = entries;
    tCal->Fill();
  }

  if (nFitted > 0) {
    TGraphErrors *grMPV =
        new TGraphErrors(nFitted, v_ch.data(), v_mpv.data(), nullptr,
                         v_mpv_err.data());
    grMPV->SetName("gr_trans_mpv_vs_channel");
    grMPV->SetTitle("Transverse Landau MPV vs channel;Channel;Landau MPV");
    grMPV->Write("", TObject::kOverwrite);

    TGraphErrors *grWidth =
        new TGraphErrors(nFitted, v_ch.data(), v_width.data(), nullptr,
                         v_width_err.data());
    grWidth->SetName("gr_trans_width_vs_channel");
    grWidth->SetTitle(
        "Transverse Landau width vs channel;Channel;Landau width");
    grWidth->Write("", TObject::kOverwrite);

    TGraph *grEntries =
        new TGraph(nFitted, v_ch.data(), v_entries.data());
    grEntries->SetName("gr_trans_entries_vs_channel");
    grEntries->SetTitle("Transverse entries vs channel;Channel;Entries");
    grEntries->Write("", TObject::kOverwrite);

    dirOut->cd();
    tCal->Write("", TObject::kOverwrite);
  } else {
    delete tCal;
  }

  outFile->cd();
  std::cout << "FitTransverseTargetChannelHistograms: fitted " << nFitted
            << " channel histograms into directory " << outputDirName
            << std::endl;
}

// ----------------------------------------------------------------------
// Longitudinal: fit gaus+gaus to per-channel target histos
// ----------------------------------------------------------------------
void FitLongitudinalTargetChannelHistograms(
    TFile *outFile,
    const std::string &inputDirName  = "target_channels",
    const std::string &outputDirName = "fitted_target_channels_longitudinal",
    int  minEntriesToFit             = 50,
    int  rebinFactor                 = 2)
{
  if (!outFile || outFile->IsZombie()) {
    std::cout << "FitLongitudinalTargetChannelHistograms: output file is null "
                 "or zombie"
              << std::endl;
    return;
  }

  TDirectory *dirIn = (TDirectory *)outFile->Get(inputDirName.c_str());
  if (!dirIn) {
    std::cout << "FitLongitudinalTargetChannelHistograms: directory "
              << inputDirName << " not found" << std::endl;
    return;
  }

  outFile->cd();
  TDirectory *dirOut = (TDirectory *)outFile->Get(outputDirName.c_str());
  if (!dirOut)
    dirOut = outFile->mkdir(outputDirName.c_str());

  if (!dirOut) {
    std::cout << "FitLongitudinalTargetChannelHistograms: cannot create output "
                 "directory "
              << outputDirName << std::endl;
    return;
  }
  dirOut->cd();

  std::vector<double> v_ch, v_mean, v_mean_err, v_sigma, v_sigma_err, v_entries;
  int nFitted = 0;

  // Optional: per-channel tree
  TTree *tCal = new TTree("long_calib", "Longitudinal per-channel calibration");
  double t_ch = 0, t_mean = 0, t_meanErr = 0, t_sigma = 0, t_sigmaErr = 0, t_entries = 0;
  tCal->Branch("channel", &t_ch, "channel/D");
  tCal->Branch("mean",    &t_mean,    "mean/D");
  tCal->Branch("meanErr", &t_meanErr, "meanErr/D");
  tCal->Branch("sigma",   &t_sigma,   "sigma/D");
  tCal->Branch("sigmaErr",&t_sigmaErr,"sigmaErr/D");
  tCal->Branch("entries", &t_entries, "entries/D");

  for (int ch = 1; ch <= NCHANNELS; ++ch) {
    TH1F *hIn = (TH1F *)dirIn->Get(Form("h_int_ch_%d", ch));
    if (!hIn)
      continue;
    if (hIn->GetEntries() < minEntriesToFit)
      continue;

    TH1F *hFit = (TH1F *)hIn->Clone(Form("h_fit_long_ch_%d", ch));
    hFit->SetDirectory(nullptr);

    if (rebinFactor > 1)
      hFit->Rebin(rebinFactor);

    int nBins = hFit->GetNbinsX();
    if (nBins < 10) {
      delete hFit;
      continue;
    }

    int bkgBin = FindMaxBinInRange(hFit, 0.0, 15000.0);
    int sigBin = FindMaxBinInRange(hFit, 22000.0, 45000.0);

    if (bkgBin < 1 || sigBin < 1) {
      delete hFit;
      continue;
    }

    double bkgAmp0  = hFit->GetBinContent(bkgBin);
    double bkgMean0 = hFit->GetBinCenter(bkgBin);
    double sigAmp0  = hFit->GetBinContent(sigBin);
    double sigMean0 = hFit->GetBinCenter(sigBin);

    if (bkgAmp0 <= 0 || sigAmp0 <= 0) {
      delete hFit;
      continue;
    }

    TF1 *fTot =
        new TF1(Form("f_long_2g_ch_%d", ch), "gaus(0)+gaus(3)", 0.0, 50000.0);
    fTot->SetParameters(bkgAmp0, bkgMean0, 3000.0, sigAmp0, sigMean0, 6000.0);

    fTot->SetParLimits(0, 0.0, 1e9);
    fTot->SetParLimits(1, 0.0, 18000.0);
    fTot->SetParLimits(2, 300.0, 12000.0);

    fTot->SetParLimits(3, 0.0, 1e9);
    fTot->SetParLimits(4, 22000.0, 45000.0);
    fTot->SetParLimits(5, 500.0, 15000.0);

    fTot->SetLineColor(kRed);
    fTot->SetLineWidth(2);

    int fitStatus = hFit->Fit(fTot, "QRS");
    if (fitStatus != 0) {
      delete hFit;
      delete fTot;
      continue;
    }

    double mean    = fTot->GetParameter(4);
    double meanErr = fTot->GetParError(4);
    double sigma   = fTot->GetParameter(5);
    double sigmaErr= fTot->GetParError(5);
    double entries = hFit->GetEntries();

    bool badFit = false;
    if (mean < 22000.0 || mean > 45000.0)
      badFit = true;
    if (sigma <= 0.0 || sigma > 15000.0)
      badFit = true;
    if (mean <= 0.0 || meanErr / mean > 0.5)
      badFit = true;

    if (badFit) {
      delete hFit;
      delete fTot;
      continue;
    }

    dirOut->cd();
    hFit->Write("", TObject::kOverwrite);
    fTot->Write("", TObject::kOverwrite);

    v_ch.push_back(ch);
    v_mean.push_back(mean);
    v_mean_err.push_back(meanErr);
    v_sigma.push_back(sigma);
    v_sigma_err.push_back(sigmaErr);
    v_entries.push_back(entries);
    ++nFitted;

    // Fill tree
    t_ch       = ch;
    t_mean     = mean;
    t_meanErr  = meanErr;
    t_sigma    = sigma;
    t_sigmaErr = sigmaErr;
    t_entries  = entries;
    tCal->Fill();
  }

  if (nFitted > 0) {
    TGraphErrors *grMean =
        new TGraphErrors(nFitted, v_ch.data(), v_mean.data(), nullptr,
                         v_mean_err.data());
    grMean->SetName("gr_long_signal_mean_vs_channel");
    grMean->SetTitle("Longitudinal signal mean vs channel;Channel;Signal mean");
    grMean->Write("", TObject::kOverwrite);

    TGraphErrors *grSigma =
        new TGraphErrors(nFitted, v_ch.data(), v_sigma.data(), nullptr,
                         v_sigma_err.data());
    grSigma->SetName("gr_long_signal_sigma_vs_channel");
    grSigma->SetTitle(
        "Longitudinal signal sigma vs channel;Channel;Signal sigma");
    grSigma->Write("", TObject::kOverwrite);

    TGraph *grEntries =
        new TGraph(nFitted, v_ch.data(), v_entries.data());
    grEntries->SetName("gr_long_entries_vs_channel");
    grEntries->SetTitle("Longitudinal entries vs channel;Channel;Entries");
    grEntries->Write("", TObject::kOverwrite);

    dirOut->cd();
    tCal->Write("", TObject::kOverwrite);
  } else {
    delete tCal;
  }

  outFile->cd();
  std::cout << "FitLongitudinalTargetChannelHistograms: fitted " << nFitted
            << " channel histograms with gaus+gaus into directory "
            << outputDirName << std::endl;
}

// ----------------------------------------------------------------------
// Wrapper functions for merged files
// ----------------------------------------------------------------------
void FitMergedTargetHistograms(
    const char *mergedFileName,
    bool        doTransverse = true,
    const char *inputDirName  = "target_channels",
    const char *outputDirName = "fitted_target_channels",
    int         minEntriesToFit = 50,
    int         rebinFactor     = 2)
{
  TransverAnalysis = doTransverse;

  TFile *f = TFile::Open(mergedFileName, "UPDATE");
  if (!f || f->IsZombie()) {
    std::cout << "FitMergedTargetHistograms: cannot open file "
              << mergedFileName << std::endl;
    return;
  }

  if (TransverAnalysis) {
    FitTransverseTargetChannelHistograms(f,
                                         inputDirName,
                                         outputDirName,
                                         minEntriesToFit,
                                         rebinFactor);
  } else {
    FitLongitudinalTargetChannelHistograms(f,
                                           inputDirName,
                                           outputDirName,
                                           minEntriesToFit,
                                           rebinFactor);
  }

  f->Write("", TObject::kOverwrite);
  f->Close();
}

void FitMergedTargetHistogramsTransverse(const char *mergedFileName,
                                         int  minEntriesToFit = 50,
                                         int  rebinFactor     = 2)
{
  FitMergedTargetHistograms(mergedFileName,
                            true,
                            "target_channels",
                            "fitted_target_channels_trans",
                            minEntriesToFit,
                            rebinFactor);
}

void FitMergedTargetHistogramsLongitudinal(const char *mergedFileName,
                                           int  minEntriesToFit = 50,
                                           int  rebinFactor     = 2)
{
  FitMergedTargetHistograms(mergedFileName,
                            false,
                            "target_channels",
                            "fitted_target_channels_long",
                            minEntriesToFit,
                            rebinFactor);
}