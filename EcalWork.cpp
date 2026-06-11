#include "MpdDataConverter.h"
#include <fstream>
#include <ostream>
#include <iostream>
#include <ctime>
#include <string>
#include <cstdlib>
#include <vector>

#include "RtypesCore.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "TSystem.h"

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(./build/libMPDDataConverter.so)
#endif

template<typename T>
TH1F* CreateWaveformHistogram(
    const std::vector<T>& data, 
    const std::string& name = "waveform",
    const std::string& title = "",
    const std::string& xTitle = "Sample Number",
    const std::string& yTitle = "ADC")
{
    if (data.empty()) return nullptr;
    
    int nSamples = data.size();

    TH1F* hist = new TH1F(name.c_str(), title.c_str(), nSamples, 0, nSamples);
    
    for (int i = 0; i < nSamples; i++) {
        hist->SetBinContent(i + 1, static_cast<Float_t>(data[i]));
    }
    
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());
    
    hist->SetLineColor(kBlue);
    hist->SetLineWidth(2);
    //hist->SetFillColor(kBlue - 9);  // Светло-синий для заливки
    hist->SetDrawOption("L");
    hist->SetBit(TH1::kNoStats);

    return hist;
}

void DrawAllWaveformsOnOneCanvas(MpdDataConverter& converter, 
                                  int eventNumber,
                                  std::vector<TH1F*>& histograms,
                                  TDirectory* eventDir) {
    
    if (histograms.empty()) return;
    
    // Создаем canvas для этого события (один pad)
    TCanvas* canvas = new TCanvas(Form("event_%d_all_waveforms", eventNumber),
                                   Form("Event %d - All Waveforms Overlay", eventNumber),
                                   1000, 600);
    
    // Настройка canvas
    canvas->SetGridx(0);
    canvas->SetGridy(0);
    
    // Цвета для разных каналов
    std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, 
                               kTeal, kPink, kViolet, kSpring, kAzure, kYellow,
                               kGray, kBlack, kOrange+1, kGreen+2, kRed+2};
    
    

    // Рисуем первую гистограмму (она создает оси)
    std::vector<Float_t> axis_y;
    switch (converter.g_TypeWF) {
        case TypeWaveForms::raw: 
            axis_y = {-35000, +40000};
            break;
        case TypeWaveForms::baseline: 
            axis_y = {-70000, 5000};
            break;
        case TypeWaveForms::invert: 
            axis_y = {-5000, 70000};
            break;
    }

    TH2F* first;

    if (!histograms.empty()) {
        first = new TH2F(Form("h_axis_%i",eventNumber),"",2,0,60,2,axis_y[0],axis_y[1]);
        first->SetLineColor(colors[0 % colors.size()]);
        first->SetBit(TH1::kNoStats);
        first->SetLineWidth(2);
        first->Draw();
        
        // Настройка осей
        first->GetXaxis()->SetTitle("Sample Number");
        first->GetYaxis()->SetTitle("ADC Value");
        first->SetTitle(Form("Event %d - All %d Channels", eventNumber, (int)histograms.size()));
    }
    
    // Рисуем остальные гистограммы поверх (same)
    for (size_t i = 0; i < histograms.size(); i++) {
        histograms[i]->SetLineColor(colors[i % colors.size()]);
        histograms[i]->SetLineWidth(2);
        histograms[i]->Draw("L same");
    }
    
    // Добавляем легенду
    TLegend* legend = new TLegend(0.85, 0.1, 0.98, 0.9);
    legend->SetBorderSize(1);
    legend->SetFillStyle(1001);
    
    for (size_t i = 0; i < histograms.size() && i < 20; i++) {  // максимум 20 в легенде
        TH1F* hist = histograms[i];
        std::string legendEntry = Form("Ch %s", histograms[i]->GetName());
        legend->AddEntry(hist, legendEntry.c_str(), "l");
    }
    
    if (histograms.size() > 20) {
        legend->AddEntry((TObject*)0, Form("+ %d more", (int)histograms.size() - 20), "");
    }
    
    legend->Draw();
    
    // Сохраняем canvas
    eventDir->cd();
    canvas->Write();
    
    // Очищаем память
    delete canvas;
    delete legend;
}

Int_t SetChannelIntegral(TH1F* hist, int binStart, int binEnd, int nBinLeft, int nBinRight) {
    if (!hist) return 0;
    // Пик
    Float_t maxValue = 0;
    Int_t peakBin=0;

    for (int i = binStart + 1; i <= binEnd; i++) {
        Float_t currentValue = TMath::Abs(hist->GetBinContent(i));
        if (currentValue > maxValue) {
            maxValue = currentValue;
            peakBin = i;
        }
    }

    double peakX = hist->GetBinCenter(peakBin);
    double peakY = hist->GetBinContent(peakBin);
    
    Int_t int_left = (peakBin-nBinLeft)<=0 ? 10 : peakBin-nBinLeft;
    Int_t int_right = (peakBin+nBinRight)>=60 ? 50 : peakBin+nBinRight;


    //std::cout<<"Integral: "<<int_left<<"\t"<<peakBin<<"\t"<<peakY<<"\t"<<int_right<<"\t"<<hist->Integral(int_left, int_right)<<std::endl;

    // Интеграл
    return hist->Integral(int_left, int_right);
    
    return 0;
}










void EcalWork(std::string inputData = "../run_rc-hs1_088.data",
              std::string outputData = "test2.root", 
              int targetEvent = -1) 
{

    TStopwatch timer1;
    timer1.Start();

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

    TH1F *h1_wf;
    
    TDirectory* WfDir = converter.outFile->mkdir("channel_wf");
    WfDir->cd();

    //converter.FullConvertData2TreeRoot();
    
    while (converter.ReadEvent() && converter.EventNumber<100) {
        
        if (targetEvent > 0 && converter.EventNumber != (uint32_t)targetEvent) {
            continue;  // Ищем дальше
        }

        // std::string eventDirName = Form("ev_%d", converter.EventNumber);
        // TDirectory* EventDir = gDirectory->GetDirectory(eventDirName.c_str());
        // if (!EventDir) {
        //     EventDir = gDirectory->mkdir(eventDirName.c_str());
        // }

        TDirectory* savedDir = gDirectory;
        std::vector<TH1F*> eventHistograms;

        while (converter.ReadADC()) {
            while(converter.ReadChannel()){
                                
                //EventDir->cd();

                int evNum=converter.channelEvent.eventNum;
                int chNum=converter.channelEvent.channelNums;
                int chPhi=converter.channelEvent.channel_Phi;
                int chZ  =converter.channelEvent.channel_Z;
                std::string h_name = Form("h_ev_%i_ch_%i_phi_%i_z_%i",evNum,chNum,chPhi,chZ);
                std::string h_title = Form("Basket 38, Event %i, channel %i, Phi %i, Z %i",evNum,chNum,chPhi,chZ);

                h1_wf=CreateWaveformHistogram(converter.channelEvent.adcValues,h_name,h_title);
                converter.channelEvent.integral=SetChannelIntegral(h1_wf,iBinStart,iBinStop,iBinLeft,iBinRight);

                if (h1_wf != nullptr) {
                    //h1_wf->Write();
                    eventHistograms.push_back(h1_wf);
                }
                
                converter.outFile->cd();
                if(targetEvent>0)converter.channelEvent.Print();
                converter.FillTreeData();
            }
        }

        if (!eventHistograms.empty()) {
            WfDir->cd();
            DrawAllWaveformsOnOneCanvas(converter, converter.EventNumber, eventHistograms, WfDir);
        }

        // for (auto hist : eventHistograms) {
        //     delete hist;
        // }
        savedDir->cd();
        
        if (targetEvent > 0) {
            break;
        }
    }

    ProgressBar(converter.inFile.tellg(), converter.inFile.tellg());
    std::cout << "Output saved to: " << outputData << std::endl;

    timer1.Stop();
    timer1.Print();

    //gSystem->Exit(0);

}