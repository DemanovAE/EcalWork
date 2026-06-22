#include "EcalDrawClass.h"
#include "MpdDataConverter.h"
#include "Rtypes.h"
#include "RtypesCore.h"
#include <string>
#include <vector>

void EcalDrawClass::SetLineColorsWF(){
    gLineColor.reserve(gColors.size() * gColorTone.size());
    for (int i = 0; i < gColorTone.size(); i++)
    {
        for (int j = 0; j < gColors.size(); j++)
        {
            gLineColor.push_back(gColors[j] + gColorTone[i]);
        }
    }
}

void EcalDrawClass::DrawUserGrid(){
    if(UserGridXY.empty()){
        for (int j = 1; j <= 12; j++){
            UserGridXY.push_back(new TLine(0.5, j + 0.5, 64.5, j + 0.5));
        }
        for (int i = 1; i <= 64; i++){
            UserGridXY.push_back(new TLine(i + 0.5, 0.5, i + 0.5, 12.5));
        }
        for (int i = 0; i < UserGridXY.size(); i++)
        {
            UserGridXY[i]->SetLineColor(kBlack);
            UserGridXY[i]->SetLineStyle(2);
            UserGridXY[i]->SetLineWidth(1);
            UserGridXY[i]->Draw("same");
        }
    }else{
        for (int i = 0; i < UserGridXY.size(); i++){
            UserGridXY[i]->Draw("same");
        }
    }
}

void EcalDrawClass::InitWfHisto(std::string xTitle, std::string yTitle){
    Int_t nBinsSamples = 60;
    Float_t SamplesMin = 0 - 0.5;
    Float_t SamplesMax = 60 - 0.5;

    eventHistograms.reserve(NChannels);
    
    for (int i = 0; i < NChannels; i++){
        eventHistograms.push_back(new TH1F(Form("h1_wf_ch%i", i+1), "", nBinsSamples, SamplesMin, SamplesMax));
        eventHistograms[i]->GetXaxis()->SetTitle(xTitle.c_str());
        eventHistograms[i]->GetYaxis()->SetTitle(yTitle.c_str());
        if(gColors.empty()){
            eventHistograms[i]->SetLineColor(kBlue);
        }else{
            eventHistograms[i]->SetLineColor(gLineColor[(i + 1) % 64]);
        }
        eventHistograms[i]->SetLineWidth(2);
        eventHistograms[i]->SetDrawOption("L");
        eventHistograms[i]->SetBit(TH1::kNoStats);
    }
}

void EcalDrawClass::FillAdcHisto(ChannelData &data){
    
    if(eventHistograms.empty())
        return;

    if(data.adcValues.empty())
        return;

    int ChNum = data.channelNums;
    int nSamples = data.adcValues.size();
    for (int i = 0; i < nSamples; i++)
    {
        eventHistograms[ChNum-1]->SetBinContent(i + 1, static_cast<Float_t>(data.adcValues[i]));
    }
}

void EcalDrawClass::InitCanvas1Pad(int TotalPict){
    gTotalDrawPict = TotalPict;
    InitCanvas1Pad();
}

void EcalDrawClass::InitCanvas2Pad(int TotalPict){
    gTotalDrawPict = TotalPict;
    InitCanvas2Pad();
}

void EcalDrawClass::InitCanvas1Pad()
{

    if(gInitCanPad==true) return;
    gInitCanPad=true;
    
    SetLineColorsWF();

    Int_t eventNumber = 0;
    
    gCanvas = new TCanvas(Form("c2"), Form("Event - All Waveforms Overlay"), 3 * 720, 1 * 720);
    gCanvas->SetBottomMargin(0.2);
    gCanvas->cd();

    gLegend = new TLegend(0.1, 0.02, 0.8, 0.12);
    gLegend->SetBorderSize(0);
    gLegend->SetNColumns(10);
    gLegend->SetFillStyle(1001);

    //legend->Draw("same");

    h2_ChIntPhiZ = new TH2D("", ";Z;Phi", 64, 0.5, 64.5, 12, 0.5, 12.5);
    h2_ChIntPhiZ->GetXaxis()->SetNdivisions(13, 5, kTRUE);
    h2_ChIntPhiZ->GetYaxis()->SetNdivisions(12, 0, kTRUE);
    h2_ChIntPhiZ->SetBit(TH1::kNoStats);
    h2_ChIntPhiZ->GetYaxis()->SetTickSize(0.005);
    h2_ChIntPhiZ->GetXaxis()->SetTitleSize(1.4 * h2_ChIntPhiZ->GetXaxis()->GetTitleSize());
    h2_ChIntPhiZ->GetXaxis()->SetLabelSize(1.4 * h2_ChIntPhiZ->GetXaxis()->GetLabelSize());
    h2_ChIntPhiZ->GetYaxis()->SetTitleSize(1.4 * h2_ChIntPhiZ->GetYaxis()->GetTitleSize());
    h2_ChIntPhiZ->GetYaxis()->SetLabelSize(1.4 * h2_ChIntPhiZ->GetYaxis()->GetLabelSize());
    h2_ChIntPhiZ->Draw("COLZ");
    // рисеут сетку
    DrawUserGrid();

    gCanvas->Modified();
    gCanvas->Update();

}

void EcalDrawClass::InitCanvas2Pad()
{

    if(gInitCanPad==true) return;
    gInitCanPad=true;

    SetLineColorsWF();
    InitWfHisto("","");

    Float_t AxisYmin = 0;
    Float_t AxisYmax = 64000;
    Int_t eventNumber = 0;
    
    gCanvas = new TCanvas(Form("c2"), Form("Event - All Waveforms Overlay"), 3 * 720, 2 * 720);

    gPadUp = new TPad("gPadUp", "Top Pad", 0, 0.35, 1, 1);
    gPadBottom = new TPad("gPadBottom", "Bottom Pad", 0, 0, 1, 0.35);

    gPadUp->SetFillColor(0);
    gPadUp->SetFrameFillColor(0);

    gPadBottom->SetFillColor(0);
    gPadBottom->SetFrameFillColor(0);
    gPadBottom->SetTopMargin(0);
    gPadBottom->SetBottomMargin(0.2);

    gCanvas->cd();
    gPadUp->Draw();
    gPadBottom->Draw();


    gAxisHistoUp = new TH2D(Form("h_axis_%i", eventNumber), "",
                               2, 0, 60, 64000, AxisYmin - 0.1 * TMath::Abs(AxisYmin - AxisYmax), AxisYmax + 0.1 * TMath::Abs(AxisYmin - AxisYmax));
    gAxisHistoUp->SetLineColor(kBlack);
    gAxisHistoUp->SetBit(TH1::kNoStats);
    gAxisHistoUp->SetLineWidth(2);
    gAxisHistoUp->GetXaxis()->SetTitle("Sample Number");
    gAxisHistoUp->GetYaxis()->SetTitle("ADC Value");
    gAxisHistoUp->SetTitle(Form("Event %d - All 0 Channels", eventNumber));

    gLegend = new TLegend(0.91, 0.1, 1.0, 0.9);
    gLegend->SetBorderSize(0);
    gLegend->SetNColumns(2);
    gLegend->SetFillStyle(1001);

    gPadBottom->cd();
    h2_ChIntPhiZ = new TH2D("h2", ";Z;Phi", 64, 0.5, 64.5, 12, 0.5, 12.5);
    h2_ChIntPhiZ->GetXaxis()->SetNdivisions(13, 5, kTRUE);
    h2_ChIntPhiZ->GetYaxis()->SetNdivisions(12, 0, kTRUE);
    h2_ChIntPhiZ->SetBit(TH1::kNoStats);
    // h2_ChIntPhiZ->GetXaxis()->SetTitleOffset(gAxisHistoUp->GetXaxis()->GetTitleOffset());
    h2_ChIntPhiZ->GetXaxis()->SetTitleSize(2 * gAxisHistoUp->GetXaxis()->GetTitleSize());
    h2_ChIntPhiZ->GetXaxis()->SetLabelSize(2 * gAxisHistoUp->GetXaxis()->GetLabelSize());
    h2_ChIntPhiZ->GetYaxis()->SetTitleSize(2 * gAxisHistoUp->GetYaxis()->GetTitleSize());
    h2_ChIntPhiZ->GetYaxis()->SetLabelSize(2 * gAxisHistoUp->GetYaxis()->GetLabelSize());
    h2_ChIntPhiZ->GetYaxis()->SetTitleOffset(0.4);
    h2_ChIntPhiZ->GetYaxis()->SetTickSize(0.005);
    // h2_ChIntPhiZ->GetXaxis()->SetTitleSize(0.5);
    h2_ChIntPhiZ->Draw("COLZ");
    // рисеут сетку
    DrawUserGrid();

    gPadBottom->Modified();
    gPadBottom->Update();

}

void EcalDrawClass::UpdateAndSaveCanvas(std::string suf, std::vector<ChannelData> &data, std::vector<int> numDraw){
    if(data.empty())return;
    std::vector<ChannelData> dataDraw;
    for (int i = 0; i < numDraw.size(); i++){
        dataDraw.push_back(data[i]);
    }
    UpdateAndSaveCanvas(suf,dataDraw);
}

void EcalDrawClass::UpdateAndSaveCanvas(std::string suf, std::vector<ChannelData> &data){
    
    if(data.empty())return;
    gLegend->Clear();

    if(gTotalDrawPict!=-1){
        if(gNumberDrawPict>=gTotalDrawPict){
            return;
        }
    }

    int eventNumber = data[0].eventNum; 

    gCanvas->SetTitle(Form("Event %i - All %i Channels",data[0].eventNum,(int)data.size()));

    if(gPadUp==nullptr){

        gCanvas->cd();
        gCanvas->SetName(Form("ChInt_PhiZ_%i",eventNumber));

        for (int i = 0; i < data.size(); i++)
        { 
            int chNum = data[i].channelNums;
            int chPhi = data[i].channel_Phi;
            int chZ = data[i].channel_Z;
            int chInt = data[i].integral;
            h2_ChIntPhiZ->SetName(Form("%s%s",h2_ChIntPhiZ->GetName(),suf.c_str()));
            h2_ChIntPhiZ->SetTitle(gCanvas->GetTitle());
    
            std::string legendEntry = Form("%i", chNum);
            gLegend->AddEntry(h2_ChIntPhiZ, legendEntry.c_str(), "l");
            h2_ChIntPhiZ->SetBinContent(h2_ChIntPhiZ->FindBin(chZ, chPhi), chInt);
        }
        gLegend->Draw("same");
        gCanvas->SaveAs(Form("pict_%s_%i%s.png", gCanvas->GetName(),eventNumber,suf.c_str()));
    }else{
        gCanvas->SetName(Form("wf_event_%i",eventNumber));
        gAxisHistoUp->SetTitle(gCanvas->GetTitle());
        gPadUp->cd();
        gPadUp->Clear();
        // Перерисовываем axis histogram
        gAxisHistoUp->Draw();

        axis_y_max.clear();
        axis_y_min.clear();
        
        h2_ChIntPhiZ->SetName(Form("%s%s",h2_ChIntPhiZ->GetName(),suf.c_str()));

        for (int i = 0; i < data.size(); i++)
        { 
            int chNum = data[i].channelNums;
            int chPhi = data[i].channel_Phi;
            int chZ = data[i].channel_Z;
            int chInt = data[i].integral;
            std::string h_name = Form("h_ev_%i_ch_%i_phi_%i_z_%i", eventNumber, chNum, chPhi, chZ);
            std::string h_title = Form("Basket 38, Event %i, channel %i, Phi %i, Z %i", eventNumber, chNum, chPhi, chZ);
            
            FillAdcHisto(data[i]);

            std::string legendEntry = Form("%i", chNum);
            gLegend->AddEntry(eventHistograms[chNum - 1], legendEntry.c_str(), "l");
            
            eventHistograms[chNum-1]->SetName(h_name.c_str());
            eventHistograms[chNum-1]->SetTitle(h_title.c_str());
            eventHistograms[chNum-1]->Draw("sameL");
    
            h2_ChIntPhiZ->SetBinContent(h2_ChIntPhiZ->FindBin(chZ, chPhi), chInt);

            axis_y_max.push_back(eventHistograms[chNum - 1]->GetMaximum());
            axis_y_min.push_back(eventHistograms[chNum - 1]->GetMinimum());
        }

        Float_t AxisYmax = *std::max_element(axis_y_max.begin(), axis_y_max.end());
        Float_t AxisYmin = *std::min_element(axis_y_min.begin(), axis_y_min.end());

        gAxisHistoUp->GetYaxis()->SetRangeUser(AxisYmin - 0.1 * TMath::Abs(AxisYmin - AxisYmax), AxisYmax + 0.1 * TMath::Abs(AxisYmin - AxisYmax));
        
        gLegend->Draw("same");
        
        gCanvas->SaveAs(Form("pict_%s_%s.png", gCanvas->GetName(),suf.c_str()));
    }

    if(gTotalDrawPict!=-1){
        gNumberDrawPict++;
    }

    if(!eventHistograms.empty()){
        for (int i = 0; i < eventHistograms.size(); i++){
            eventHistograms[i]->Reset();
        }
    }

    h2_ChIntPhiZ->Reset();
}

EcalDrawClass::~EcalDrawClass(){
    
    delete gCanvas;    
    delete gLegend;

    gCanvas = nullptr;
    gLegend = nullptr;
    gAxisHistoUp = nullptr;
    gAxisHistoBottom = nullptr;
    h2_ChIntPhiZ = nullptr;

    if(!eventHistograms.empty()){
        for (auto hist : eventHistograms) {
            if (hist) {
                delete hist;
                hist = nullptr;
            }
        }
        eventHistograms.clear();  // Очищаем вектор (опционально)
    }

    for (auto obj : UserGridXY) {
        if (obj) {
            delete obj;
            obj = nullptr;
        }
    }
    UserGridXY.clear();  // Очищаем вектор (опционально)

}