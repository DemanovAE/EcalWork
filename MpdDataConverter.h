#ifndef MPD_DATA_CONVERTER_H
#define MPD_DATA_CONVERTER_H

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>
#include "RtypesCore.h"
#include "TH1.h"


class TFile;
class TTree;

enum class TypeReadByte { unknown, event, adc, channel, chSamples };
enum class TypeWaveForms { raw, baseline, invert};

// Структура 1го event в дереве = 1 канал
struct ChannelData {
    // read
    UInt_t eventNum = 0;
    UInt_t timestamp_sec = 0;
    UInt_t timestamp_ns = 0;
    ULong64_t timestamp_f = 0;
    UInt_t deviceID = 0;
    Short_t channel_Phi = 0;
    Short_t channel_Z = 0;
    Short_t channelNums = 0;
    Short_t numSamples = 0;
    std::vector<Int_t> adcValues = {};

    // calculated
    Int_t integral = 0;
    Int_t peak = 0;
    Int_t amplitude = 0;
    Int_t pedestal = 0;
    Float_t baseline = 0.;

    void Clear();
    void Print(bool printADC = false, int maxADC = 20) const;
};

// Класс для чтения бинарного файла и записи в root tree
class MpdDataConverter {

    private:

        const uint32_t SYNC_WORD = 0x2A50D5AF;
        bool print4ByteFormat = true; // выводить счетчик в байтах или счетчик по 4 байт

        uint64_t g_CountReadFileByte = 0;
        uint64_t g_CountReadEventByte = 0;
        uint64_t g_CountReadAdcByte = 0;
        uint64_t g_CountReadChannelByte = 0;
        uint64_t g_CountReadChSamplesByte = 0;

        uint32_t g_LenFile = 0;
        uint32_t g_LenEvent = 0;
        uint32_t g_LenAdc = 0;
        uint32_t g_LenChannel = 0;
        uint32_t g_LenChannelSamples = 0;
        
        //Функция чтения по 4 байта
        uint32_t readWord(TypeReadByte type);

    public:

        bool wChSamplesVector=false;

        TypeWaveForms g_TypeWF = TypeWaveForms::raw;
        Int_t wf_add_val = 30000;
        Int_t wf_scale = -1;

        TFile* outFile = nullptr;
        std::ifstream inFile;
        TTree* tree = nullptr;
        ChannelData channelEvent;
        uint32_t ChannelNumAndLen;
        uint32_t EventNumber;
        uint32_t deviceSerial;
        uint32_t deviceTimeStampSec;
        uint32_t deviceTimeStampSecNs;
        uint32_t unclearUsage;

        uint32_t devIdAndLen;
        uint32_t CountEventNumber = 0;
        
        void WriteChannelSamplesVector(){this->wChSamplesVector=true;};
        void SetTypeWaveForms(TypeWaveForms _type){ this->g_TypeWF=_type; };

        void OpenInputFile(std::string inputData);
        void OpenOutRootFile(std::string outputData);
        void InitTree();
        void FillTreeData();
        void WriteTreeAndClose();
        //Для проверки
        void PrintCurrentCounByte();
        void PrintWord(uint32_t word, const std::string& label = "");
        //Функция пропуска указанного кол-ва байтов байт
        void SkipBytes(uint64_t bytesToSkip);
        //Это если просто надо конвертнуть данные в root дерево без обработки.
        void FullConvertData2TreeRoot();
        void showProgress(uint64_t current, uint64_t total, int barWidth = 45);

        bool ReadEvent();
        bool ReadADC();
        bool ReadChannel();

        ULong64_t SetTimeStampNs(uint32_t _sec, uint32_t _ns);
        Short_t SetChannelNum(uint32_t _AdcId, uint32_t _chNum);
        Int_t AdcSampleValueFirst(uint32_t _adcValue);
        Int_t AdcSampleValueSecond(uint32_t _adcValue);
        Int_t ChangeWfValue(Int_t _value);

        ~MpdDataConverter();

};

// Вспомогательные функции для маппы каналов
void InitChannelMap(const std::string& filename = "basket_channel_map_phiZ.csv");
Short_t GetPhi(short channel);
Short_t GetZ(short channel);
Short_t getAddressForBasket_38(uint32_t channel);

// Прогресс-бар
void ProgressBar(uint64_t current, uint64_t total);

#endif