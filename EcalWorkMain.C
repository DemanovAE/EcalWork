#include "EcalDrawClass.h"
#include "MpdDataConverter.h"
#include "Rtypes.h"
#include <cstdlib>
#include <iostream>

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

void EcalWork(std::string inputDataTree, std::string outputData,
              Long64_t firstEntry, Long64_t lastEntry, const EcalConfig &cfg);

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
              << " input.root output.root [firstEntry] [lastEntry]"
              << " [transFlag] [useLongQA] [useTransQA]"
              << " [longMax1] [longMax2] [min3x3] [max3x3] [maxDiff5x5]"
              << " [transAmp1] [transAmp2] [transMinLen] [transContam]\n";
    return 1;
  }

  std::string input = argv[1];
  std::string output = argv[2];

  Long64_t firstEntry = 0;
  Long64_t lastEntry = -1;

  if (argc >= 4)
    firstEntry = std::stoll(argv[3]);
  if (argc >= 5)
    lastEntry = std::stoll(argv[4]);

  // Defaults (same as now)
  EcalConfig cfg;
  cfg.transAnalysis = false;
  cfg.useLongQA = true;
  cfg.useTransQA = false;
  cfg.long_max_1st_Integral = 3500;
  cfg.long_max_2nd_Integral = 1500;
  cfg.long_min_3x3_ratio = 0.45f;
  cfg.long_max_3x3_ratio = 0.99f;
  cfg.long_max_diffusivity_5x5 = 0.35f;

  // Transverse defaults
  cfg.trans_amp_thr1 = 100.0;
  cfg.trans_amp_thr2 = 500.0;
  cfg.trans_min_strip_len = 5;
  cfg.trans_contam_frac = 0.20; // reserved, not yet used

  int arg = 5;
  if (argc > arg)
    cfg.transAnalysis = (std::atoi(argv[arg]) != 0);
  ++arg;
  if (argc > arg)
    cfg.useLongQA = (std::atoi(argv[arg]) != 0);
  ++arg;
  if (argc > arg)
    cfg.useTransQA = (std::atoi(argv[arg]) != 0);
  ++arg;
  if (argc > arg)
    cfg.long_max_1st_Integral = std::atoi(argv[arg]);
  ++arg;
  if (argc > arg)
    cfg.long_max_2nd_Integral = std::atoi(argv[arg]);
  ++arg;
  if (argc > arg)
    cfg.long_min_3x3_ratio = std::atof(argv[arg]);
  ++arg;
  if (argc > arg)
    cfg.long_max_3x3_ratio = std::atof(argv[arg]);
  ++arg;
  if (argc > arg)
    cfg.long_max_diffusivity_5x5 = std::atof(argv[arg]);
  ++arg;
  if (argc > arg)
    cfg.trans_amp_thr1 = std::atof(argv[arg]);
  ++arg;
  if (argc > arg)
    cfg.trans_amp_thr2 = std::atof(argv[arg]);
  ++arg;
  if (argc > arg)
    cfg.trans_min_strip_len = std::atoi(argv[arg]);
  ++arg;
  if (argc > arg)
    cfg.trans_contam_frac = std::atof(argv[arg]);
  ++arg;

  EcalWork(input, output, firstEntry, lastEntry, cfg);
  return 0;
}