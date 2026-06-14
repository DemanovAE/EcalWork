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


#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(./build/libMPDDataConverter.so)
#endif

template<typename T>
void SetWaveformHistogram(
    const std::vector<T>& data, 
    TH1F *hist,
    const std::string& name = "waveform",
    const std::string& title = "",
    const std::string& xTitle = "Sample Number",
    const std::string& yTitle = "ADC")
{
    if (data.empty()) return;
    
    int nSamples = data.size();

    //TH1F* hist = new TH1F(name.c_str(), title.c_str(), nSamples, 0, nSamples);
    
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

void DrawAllWaveformsOnOneCanvas( std::vector<ChannelData> _ChData,
                                  std::vector<TH1F*>& histograms,
                                  TH2D * h2_histo,
                                  std::vector<TLine*>& _userGrid,
                                  TDirectory* eventDir) {
    
    if (_ChData.empty()) return;

    int eventNumber = _ChData[0].eventNum;

    std::vector<Float_t> axis_y_min;
    std::vector<Float_t> axis_y_max;

    for (size_t i = 0; i < _ChData.size(); i++) {
        int chNum=_ChData[i].channelNums;
        int chPhi=_ChData[i].channel_Phi;
        int chZ  =_ChData[i].channel_Z;
        int chInt =_ChData[i].integral;
        std::string h_name = Form("h_ev_%i_ch_%i_phi_%i_z_%i",eventNumber,chNum,chPhi,chZ);
        std::string h_title = Form("Basket 38, Event %i, channel %i, Phi %i, Z %i",eventNumber,chNum,chPhi,chZ);
        
        SetWaveformHistogram(_ChData[i].adcValues,histograms[chNum-1],h_name,h_title);
        
        h2_histo->SetBinContent(h2_histo->FindBin(chZ,chPhi),chInt);
    
        axis_y_max.push_back(histograms[chNum-1]->GetMaximum());
        axis_y_min.push_back(histograms[chNum-1]->GetMinimum());

    }

    Float_t AxisYmax= *std::max_element(axis_y_max.begin(), axis_y_max.end());
    Float_t AxisYmin= *std::min_element(axis_y_min.begin(), axis_y_min.end());

    TCanvas* canvas = new TCanvas(Form("wf_event_%d", eventNumber),Form("Event %d - All Waveforms Overlay", eventNumber), 3*720,2*720);

    TPad* topPad = new TPad("topPad", "Top Pad", 0, 0.35, 1, 1);
    TPad* bottomPad = new TPad("bottomPad", "Bottom Pad", 0, 0, 1, 0.35);

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
    std::vector<int> colorTone = {-2,-1,0,1,2,3};
    std::vector<int> LineColor;

    LineColor.reserve(colors.size()*colorTone.size());
    for(int i=0; i<colorTone.size();i++){
        for(int j=0; j<colors.size();j++){
            LineColor.push_back(colors[j]+colorTone[i]);
        }
    }

    TH2F* AxisHisto = new TH2F(Form("h_axis_%i",eventNumber),"",
        2,0,60,2,AxisYmin-0.1*TMath::Abs(AxisYmin-AxisYmax),AxisYmax+0.1*TMath::Abs(AxisYmin-AxisYmax));
    AxisHisto->SetLineColor(colors[0 % colors.size()]);
    AxisHisto->SetBit(TH1::kNoStats);
    AxisHisto->SetLineWidth(2);
    topPad->cd();
    AxisHisto->Draw();
    
    AxisHisto->GetXaxis()->SetTitle("Sample Number");
    AxisHisto->GetYaxis()->SetTitle("ADC Value");
    AxisHisto->SetTitle(Form("Event %d - All %d Channels", eventNumber, (int)_ChData.size()));
    
    TLegend* legend = new TLegend(0.91, 0.1, 1.0, 0.9);
    legend->SetBorderSize(0);
    legend->SetNColumns(2);
    legend->SetFillStyle(1001);

    for (size_t i = 0; i < _ChData.size(); i++) {
        int iCh = _ChData[i].channelNums-1;
        if(histograms[iCh]->GetEntries()<1)continue;
        //histograms[iCh]->SetLineColor(colors[i % colors.size()]);
        histograms[iCh]->SetLineColor(LineColor[(iCh+1)%64]);
        histograms[iCh]->SetLineWidth(2);
        histograms[iCh]->Draw("L same");
        std::string legendEntry = Form("%i", _ChData[i].channelNums);
        legend->AddEntry(histograms[iCh], legendEntry.c_str(), "l");
    }
        
    legend->Draw("same");
    
    bottomPad->cd();
    //h2_histo->GetXaxis()->SetTitleOffset(AxisHisto->GetXaxis()->GetTitleOffset());
    h2_histo->GetXaxis()->SetTitleSize(2*AxisHisto->GetXaxis()->GetTitleSize());
    h2_histo->GetXaxis()->SetLabelSize(2*AxisHisto->GetXaxis()->GetLabelSize());
    h2_histo->GetYaxis()->SetTitleSize(2*AxisHisto->GetYaxis()->GetTitleSize());
    h2_histo->GetYaxis()->SetLabelSize(2*AxisHisto->GetYaxis()->GetLabelSize());
    h2_histo->GetYaxis()->SetTitleOffset(0.4); 
    h2_histo->GetYaxis()->SetTickSize(0.005);
    //h2_histo->GetXaxis()->SetTitleSize(0.5);
    h2_histo->Draw("COLZ");
  
    for(int i=0; i<_userGrid.size();i++){
        _userGrid[i]->SetLineColor(kBlack);
        _userGrid[i]->SetLineStyle(2);
        _userGrid[i]->SetLineWidth(1);
        _userGrid[i]->Draw("same");
    }

    bottomPad->Modified();
    bottomPad->Update();

    eventDir->cd();
    //canvas->Write();
    canvas->SaveAs(Form("/home/aleksandr/ecal_work/pict/%s.png",canvas->GetName()));

    for (size_t i = 0; i < _ChData.size(); i++) {
        int chNum=_ChData[i].channelNums;
        histograms[chNum-1]->Reset();
    }

    h2_histo->Reset();

    delete canvas;
    delete legend;
    delete AxisHisto;

    canvas = nullptr;
    legend = nullptr;
    AxisHisto = nullptr;
}

Int_t SetPedestal(std::vector<Int_t> _adc_value, Int_t _Skip, Int_t _count){
    Int_t SumAdc = 0;
    Int_t nCount = 0;
    Int_t result = 0;
    if(_adc_value.size()==0)return 0;
    
    for(int i=_Skip; i<(_Skip+_count); i++){
        SumAdc+=_adc_value[i];
        nCount++;
    }
    result = (Int_t)(SumAdc/nCount);
    return result;
}

void PedestalSubtraction(std::vector<Int_t>& _ch_Data, Int_t _pedestal){
    for(int i=0; i<_ch_Data.size();i++){
        _ch_Data[i] = _ch_Data[i]-_pedestal;
    }
}

Int_t GetPedestalAmpl(std::vector<Int_t> _adc_value, Int_t _Skip, Int_t _count){
    Int_t SumAdc = 0;
    Int_t nCount = 0;
    Int_t result = 0;
    if(_adc_value.size()==0)return 0;
    
    auto [minIt, maxIt] = std::minmax_element(_adc_value.begin(), _adc_value.begin() + _Skip + _count);

    return TMath::Abs(*maxIt-*minIt);
}

void SetChannelIntegralAmp(ChannelData& _chData, TH1F* hist, int binStart, int binEnd, int nBinLeft, int nBinRight) {
    if (!hist) return;
    if(hist->GetEntries()==0) return;
    // Пик
    Float_t maxValue = 0;
    Int_t peakBin=0;

    for (int i = binStart; i <= binEnd; i++) {
        Float_t currentValue = TMath::Abs(hist->GetBinContent(i));
        if (currentValue > maxValue) {
            maxValue = currentValue;
            peakBin = i;
        }
    }
    
    Int_t int_left = (peakBin-nBinLeft)<=0 ? 20 : peakBin-nBinLeft;
    Int_t int_right = (peakBin+nBinRight)>=60 ? 46 : peakBin+nBinRight;

    double peakX = hist->GetBinCenter(peakBin);
    double peakY = hist->GetBinContent(peakBin);

    // Интеграл
    _chData.integral = hist->Integral(int_left, int_right);
    _chData.amplitude = peakY;
    //_chData.Print();
}

void ResetTH1Fvector(std::vector<TH1F*>& hists, std::vector<int> iNumClear) {
    for (int i=0; i<iNumClear.size();i++) {
        hists[iNumClear[i]]->Reset();
    }
}








void EcalWork(std::string inputData = "../run_rc-hs1_088.data",
              std::string outputData = "../test_all2.root", 
              int targetEvent = -1) 
{

    TStopwatch timer1;
    timer1.Start();

    // Число первых бинов для пропуска и число последующих бинов для подсчета подложки
    Int_t iBinPedestalSkip = 2;
    Int_t iBinPedestalCount = 10;
    //Для поиска пика
    Int_t iBinStart = 20;
    Int_t iBinStop  = 50;
    //Интеграл вокруг пика
    Int_t iBinLeft = 4;
    Int_t iBinRight = 16;

    InitChannelMap("basket_channel_map_phiZ.csv");

    MpdDataConverter converter;
    
    converter.OpenInputFile(inputData);
    converter.OpenOutRootFile(outputData);
    //converter.WriteChannelSamplesVector();
    converter.SetTypeWaveForms(TypeWaveForms::invert);
    converter.InitTree();

    Int_t nBinsSamples = 60;
    Float_t SamplesMin = 0-0.5;
    Float_t SamplesMax = 60-0.5;
    
    TH1F *h1_wf = new TH1F("hist", "", nBinsSamples, SamplesMin, SamplesMax);

    Int_t MaxNChannels = 64*12;
    std::vector<TH1F*> eventHistograms;
    eventHistograms.reserve(MaxNChannels);
    for(int i=0; i<MaxNChannels;i++){
        eventHistograms.push_back(new TH1F(Form("histo_%i",i), "", nBinsSamples, SamplesMin, SamplesMax));
    }
    TH2D *h2_int_ch = new TH2D(Form("ch_integral"), ";channel Num;Integral", MaxNChannels,0.5,MaxNChannels+0.5,700,0,700000);


    std::vector<ChannelData> ChannelDataInEvent;
    ChannelDataInEvent.reserve(MaxNChannels);

    TH2D *h2_integral_z_phi=new TH2D("h2",";Z;Phi",64,0.5,64.5,12,0.5,12.5);
    h2_integral_z_phi->GetXaxis()->SetNdivisions(13, 5, kTRUE);
    h2_integral_z_phi->GetYaxis()->SetNdivisions(12, 0, kTRUE);
    h2_integral_z_phi->SetBit(TH1::kNoStats);

    std::vector<TLine*> UserGridXY;
    for (int j = 1; j <= 12; j++) {
        UserGridXY.push_back(new TLine(0.5, j + 0.5, 64.5, j + 0.5));
        
    }
    for (int i = 1; i <= 64; i++) {
        UserGridXY.push_back(new TLine(i + 0.5, 0.5, i + 0.5, 12.5));
    }

    TDirectory* WfDir = converter.outFile->mkdir("channel_wf");
    WfDir->cd();
    
    while (converter.ReadEvent() && converter.EventNumber<1000){
        
        if (targetEvent > 0 && converter.EventNumber != (uint32_t)targetEvent) {
            continue;
        }

        ChannelDataInEvent.clear();
        std::vector<int> _PhiCount(12);

        while (converter.ReadADC()) {
            while(converter.ReadChannel()){
                                
                //Считаем подложку и вычитем ее
                converter.channelEvent.pedestal=SetPedestal(converter.channelEvent.adcValues,iBinPedestalSkip,iBinPedestalCount);
                PedestalSubtraction(converter.channelEvent.adcValues,converter.channelEvent.pedestal);
                
                //Далее гистограмм для WF и подсчета интеграла
                int evNum=converter.channelEvent.eventNum;
                int chNum=converter.channelEvent.channelNums;
                int chPhi=converter.channelEvent.channel_Phi;
                int chZ  =converter.channelEvent.channel_Z;
                std::string h_name = Form("h_ev_%i_ch_%i_phi_%i_z_%i",evNum,chNum,chPhi,chZ);
                std::string h_title = Form("Basket 38, Event %i, channel %i, Phi %i, Z %i",evNum,chNum,chPhi,chZ);

                SetWaveformHistogram(converter.channelEvent.adcValues,h1_wf,h_name,h_title);                
                SetChannelIntegralAmp(converter.channelEvent, h1_wf,iBinStart,iBinStop,iBinLeft,iBinRight);
                
                if(converter.channelEvent.amplitude<100)continue;
                _PhiCount.at(chPhi-1)++;

                ChannelDataInEvent.push_back(converter.channelEvent);                
                
                //converter.outFile->cd();
                //converter.FillTreeData();
            } //end ReadChannel
        } //end ReadADC

        WfDir->cd();
        
        if(ChannelDataInEvent.size() < 5) continue;
        if(ChannelDataInEvent.size() > 128) continue;

        int max_PhiCount = *std::max_element(_PhiCount.begin(), _PhiCount.end());
        if( max_PhiCount < 5) continue;

        for(int i=0; i<ChannelDataInEvent.size();i++){
            converter.channelEvent = ChannelDataInEvent[i];
            converter.FillTreeData();
            h2_int_ch->Fill(converter.channelEvent.channelNums,converter.channelEvent.integral);
        }
        DrawAllWaveformsOnOneCanvas(ChannelDataInEvent, eventHistograms,h2_integral_z_phi, UserGridXY, WfDir);

        if (targetEvent > 0) {
            break;
        }
    }

    converter.outFile->cd();
    h2_int_ch->Write();

    //ProgressBar(converter.inFile.tellg(), converter.inFile.tellg());
    std::cout << "\n\nOutput saved to: " << outputData <<"\n\n"<< std::endl;

    converter.WriteTreeAndClose();

    timer1.Stop();
    timer1.Print();

}