#ifndef ECAL_DRAW_CLASS_H
#define ECAL_DRAW_CLASS_H

#include "MpdDataConverter.h"

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMathBase.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TTreeReader.h"
#include "TSystem.h"
#include "TF1.h"
#include "TGraphErrors.h"

class TFile;
class TTree;

// Класс для чтения бинарного файла и записи в root tree
class EcalDrawClass {

    public:

    bool gInitCanPad = false;
    Int_t gNumberDrawPict = 0;
    Int_t gTotalDrawPict = -1;
    Int_t gMaxPhi = 12;
    Int_t gMaxZ = 64;
    Int_t NChannels = 768;
    
    Int_t xCanvasSize = 3 * 720;
    Int_t yCanvasSize = 2 * 720;

    std::vector<Float_t> axis_y_min;
    std::vector<Float_t> axis_y_max;

    TCanvas *gCanvas = nullptr;
    TPad *gPadUp = nullptr;
    TPad *gPadBottom = nullptr;
    TH2D *gAxisHistoUp = nullptr;
    TH2D *gAxisHistoBottom = nullptr;
    TH2D *h2_ChIntPhiZ = nullptr;
    TLegend *gLegend = nullptr;

    std::vector<TLine *> UserGridXY;

    std::vector<int> gColors = {kRed, kBlue, kGreen, kMagenta, kCyan, kYellow, kTeal, kPink, kViolet, kSpring, kAzure, kOrange};
    std::vector<int> gColorTone = {-2, -1, 0, 1, 2, 3};
    std::vector<int> gLineColor;

    std::vector<TH1F *> eventHistograms;

    void InitCanvas2Pad();
    void InitCanvas2Pad(int TotalPict);
    void InitCanvas1Pad();
    void InitCanvas1Pad(int TotalPict);

    void UpdateAndSaveCanvas(std::string suf, std::vector<ChannelData> &data);
    void UpdateAndSaveCanvas(std::string suf, std::vector<ChannelData> &data, std::vector<int> numDraw);

    void SetLineColorsWF();
    void InitWfHisto(std::string xTitle, std::string yTitle);
    void FillAdcHisto(ChannelData &data);
    void DrawUserGrid();

    ~EcalDrawClass();

};

#endif