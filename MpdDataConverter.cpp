#include "MpdDataConverter.h"
#include "RtypesCore.h"
#include <TFile.h>
#include <TTree.h>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>

using namespace std;

// ==================== ChannelData ====================

void ChannelData::Clear() {
    eventNum = 0;
    timestamp_sec = 0;
    timestamp_ns = 0;
    timestamp_f = 0;
    deviceID = 0;
    channel_Phi = 0;
    channel_Z = 0;
    channelNums = 0;
    numSamples = 0;
    adcValues.clear();
 
    integral = 0;
    peak = 0;
    baseline = 0;
    amplitude = 0;
    pedestal = 0;
    
}

void ChannelData::Print(bool printADC, int maxADC) const {
    std::cout << "\n========== CHANNEL DATA ==========" << std::endl;
    std::cout << std::setw(25) << "Event Number:  " << eventNum << std::endl;
    std::cout << std::setw(25) << "Device ID:  " << std::hex << deviceID << std::dec << std::endl;
    std::cout << std::setw(25) << "Channel Number:  " << channelNums << std::endl;
    std::cout << std::setw(25) << "Channel Phi:  " << channel_Phi << std::endl;
    std::cout << std::setw(25) << "Channel Z:  " << channel_Z << std::endl;
    std::cout << std::setw(25) << "Timestamp (sec):  " << timestamp_sec << " s" << std::endl;
    std::cout << std::setw(25) << "Timestamp (ns):  " << timestamp_ns << " ns" << std::endl;
    std::cout << std::setw(25) << "Timestamp (sec+ns):  " << timestamp_f << " ns" << std::endl;
    std::cout << std::setw(25) << "Number of Samples:  " << numSamples << std::endl;
    std::cout << std::setw(25) << "Integral:  " << integral << std::endl;
    std::cout << std::setw(25) << "Peak:  " << peak << std::endl;
    std::cout << std::setw(25) << "Amplitude:  " << amplitude << std::endl;
    std::cout << std::setw(25) << "Pedestal:  " << pedestal << std::endl;
    std::cout << std::setw(25) << "Baseline:  " << std::fixed << std::setprecision(2) << baseline << std::endl;
    
    if (printADC && !adcValues.empty()) {
        std::cout << std::setw(25) << "ADC Samples:" << std::endl;
        int count = std::min(maxADC, (int)adcValues.size());
        for (int i = 0; i < count; i++) {
            std::cout << "  [" << std::setw(4) << i << "] " << adcValues[i];
            if ((i + 1) % 10 == 0 && i + 1 < count) std::cout << std::endl;
        }
        if (adcValues.size() > maxADC) {
            std::cout << "  ... and " << (adcValues.size() - maxADC) << " more" << std::endl;
        }
    }
    std::cout << "==================================" << std::endl;
}

Int_t ChannelData::GetMaximumAdcValue(Int_t startIndex, Int_t endIndex) {
    if (adcValues.empty()) return 0;
    if (startIndex >= adcValues.size()) return 0;
    if (endIndex > adcValues.size()) endIndex = adcValues.size();
    if (startIndex >= endIndex) return 0;
    return *std::max_element(adcValues.begin() + startIndex, adcValues.begin() + endIndex);
}

Int_t ChannelData::GetMinimumAdcValue(Int_t startIndex, Int_t endIndex) {
    if (adcValues.empty()) return 0;
    if (startIndex >= adcValues.size()) return 0;
    if (endIndex > adcValues.size()) endIndex = adcValues.size();
    if (startIndex >= endIndex) return 0;
    return *std::min_element(adcValues.begin() + startIndex, adcValues.begin() + endIndex);
}

// ==================== Channel Map ====================

//маппа и функции для возврата координаты x и phi
static std::map<short, std::pair<short, short>> g_channelMap;

void InitChannelMap(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout<<"File "<<filename<<" not found!"<<std::endl;
        return;
    };
    
    std::string line;
    std::getline(file, line);  // Пропускаем заголовок
    
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string ch, phi, z;
        std::getline(ss, ch, ',');
        std::getline(ss, phi, ',');
        std::getline(ss, z, ',');
        
        g_channelMap[std::stoi(ch)] = {std::stoi(phi), std::stoi(z)};
    }
}

Short_t GetPhi(short channel) {
    auto it = g_channelMap.find(channel);
    return it != g_channelMap.end() ? it->second.first : -1;
}

Short_t GetZ(short channel) {
    auto it = g_channelMap.find(channel);
    return it != g_channelMap.end() ? it->second.second : -1;
}

// Прогресс-бар
void ProgressBar(uint64_t current, uint64_t total) {
    static int lastPercent = -1;
    
    if (total == 0) return;
    
    int percent = (current * 100) / total;
    
    // Обновляем только при изменении процента
    if (percent != lastPercent) {
        int barWidth = 45;
        int pos = barWidth * percent / 100;
        
        std::cout << "\r[";
        for (int i = 0; i < barWidth; i++) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        
        std::cout << "] " << std::setw(3) << percent << "% "
                  << "(" << current / 1024 / 1024 << "/" << total / 1024 / 1024 << " MB)";
        std::cout.flush();
        
        lastPercent = percent;
        
        if (current == total) {
            std::cout << std::endl;
        }
    }
}

Short_t getAddressForBasket_38(uint32_t channel)
{
    switch (channel)
    {
    case 0x097fec5f:
        return 1;
    case 0x097f5e39:
        return 2;
    case 0x0d6c59c3:
        return 3;
    case 0x097f65eb:
        return 4;
    case 0x09803a84:
        return 5;
    case 0x097fd82a:
        return 6;
    case 0x0980ae1a:
        return 7;
    case 0x097f5243:
        return 8;
    case 0x0cd13128:
        return 9;
    case 0x0e74f06a:
        return 10;
    case 0x0cd130ff:
        return 11;
    case 0x0cd0a3ac:
        return 12;
    default:
        return 0;
    }
}




// ==================== MpdDataConverter ====================

void MpdDataConverter::WriteTreeAndClose() {
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Warning: Output file is not open" << std::endl;
        return;
    }
    
    if (!tree) {
        std::cerr << "Warning: Tree is null" << std::endl;
        return;
    }
    
    Long64_t nEntries = tree->GetEntries();
    if (nEntries == 0) {
        std::cerr << "Warning: Tree is empty, nothing to write" << std::endl;
        return;
    }
    
    outFile->cd();
    
    std::cout << "Writing tree (" << nEntries << " entries) to file..." << std::endl;
    tree->Write("", TObject::kOverwrite);
    std::cout << "Tree written successfully: " << nEntries << " entries" << std::endl;
}

MpdDataConverter::~MpdDataConverter() {
    if (inFile.is_open()) {
        inFile.close();
        std::cout << "Input file closed" << std::endl;
    }
    if (outFile) {
        try {
            if (!outFile->IsZombie()) {
                outFile->Close();
                std::cout << "ROOT file closed successfully" << std::endl;
            } else {
                std::cerr << "Warning: ROOT file is zombie, skipping..." << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error closing ROOT file: " << e.what() << std::endl;
        }
        
        // Удаляем объект файла
        delete outFile;
        outFile = nullptr;
        tree = nullptr;  // tree удаляется вместе с outFile
    }

}

// Открытие бинарного файла
void MpdDataConverter::OpenInputFile(string inputData){
    inFile.open(inputData, std::ios::binary);
    if (!inFile.is_open()) {
        std::cerr << "ERROR: Could not open binary file: " << inputData << std::endl;
        return;
    }
    // Сохраняем текущую позицию
    std::streampos currentPos = inFile.tellg();
    // Перемещаемся в конец
    inFile.seekg(0, std::ios::end);
    // Получаем размер
    g_LenFile = inFile.tellg();
    // Возвращаемся обратно
    inFile.seekg(currentPos);
    std::cout << "File opened: " << inputData << " (" << g_LenFile << " bytes)" << std::endl;
}

// Открытие ROOT файла
void MpdDataConverter::OpenOutRootFile(string outputData){
    outFile = TFile::Open(outputData.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "ERROR: Could not create ROOT file: " << outputData << std::endl;
        inFile.close();
    }
}

// Создаем дерево
void MpdDataConverter::InitTree(){
    tree = new TTree("events", "Ecal Data Tree");
    tree->Branch("eventNum", &channelEvent.eventNum, "eventNum/i");
    tree->Branch("deviceID", &channelEvent.deviceID, "deviceID/i");
    tree->Branch("timestamp_sec", &channelEvent.timestamp_sec, "timestamp_sec/i");
    tree->Branch("timestamp_ns", &channelEvent.timestamp_ns, "timestamp_ns/i");
    tree->Branch("timestamp", &channelEvent.timestamp_f, "timestamp/l");
    tree->Branch("channelNums", &channelEvent.channelNums,"channelNums/s");
    tree->Branch("channel_Phi", &channelEvent.channel_Phi,"channel_Phi/s");
    tree->Branch("channel_Z", &channelEvent.channel_Z,"channel_Z/s");
    tree->Branch("integral", &channelEvent.integral, "integral/I");
    tree->Branch("peak", &channelEvent.peak, "peak/I");
    tree->Branch("amplitude", &channelEvent.amplitude, "amplitude/I");
    tree->Branch("pedestal", &channelEvent.pedestal, "pedestal/I");
    if(wChSamplesVector)tree->Branch("adcValues", &channelEvent.adcValues);
}

void MpdDataConverter::FillTreeData(){
    tree->Fill();
}

//Это если просто надо конвертнуть данные в root дерево без обраотки.
void MpdDataConverter::FullConvertData2TreeRoot(){
    while (ReadEvent()) {
        while (ReadADC()) {
            while(ReadChannel()){
                FillTreeData(); 
            }
        }
    }
}

// Для проверки
void MpdDataConverter::PrintCurrentCounByte(){       
    uint32_t word;
    uint64_t scale = print4ByteFormat ? 1 : 4;
    
    std::cout << "\n\n==================== READ STATISTICS ===================" << std::endl;
    if(CountEventNumber!=0)std::cout<<std::setw(24)<<"Current event number:  "<<std::setw(12)<<CountEventNumber<< std::endl;
    std::cout<<std::setw(24)<<"Total bytes:[ "<<std::setw(12)<<g_CountReadFileByte/scale<<" / "<<std::setw(12)<<g_LenFile/scale << " ]"<<std::endl;
    std::cout<<std::setw(24)<<"Event bytes:[ "<<std::setw(12)<<g_CountReadEventByte/scale<<" / "<<std::setw(12)<<g_LenEvent/scale << " ]"<<std::endl;
    std::cout<<std::setw(24)<<"ADC bytes:[ "<<std::setw(12)<<g_CountReadAdcByte/scale<<" / "<<std::setw(12)<<g_LenAdc/scale << " ]"<<std::endl;
    std::cout<<std::setw(24)<<"Channel bytes:[ "<<std::setw(12)<<g_CountReadChannelByte/scale<<" / "<<std::setw(12)<<g_LenChannel/scale << " ]"<<std::endl;
    std::cout<<std::setw(24)<<"Ch samples bytes:[ "<<std::setw(12)<<g_CountReadChSamplesByte/scale<<" / "<<std::setw(12)<<g_LenChannelSamples/scale << " ]"<<std::endl;
    inFile.seekg(-4, std::ios::cur);
    inFile.read(reinterpret_cast<char*>(&word), sizeof(word));
    PrintWord(word,"Last 4 bytes read:  ");
    std::cout << "========================================================" << std::endl;
}

void MpdDataConverter::PrintWord(uint32_t word, const std::string& label) {
    if (!label.empty()) {
        std::cout<<std::setw(24)<<label;
    }
    std::cout << "0x" << std::hex<<std::setw(8)<<std::setfill('0')<<word<< std::dec;
    // Дополнительно выводим как число
    std::cout << " (" << word << ")" << std::endl;
    std::cout << std::setfill(' ');
}

//прогресс чтения исходного файла
void MpdDataConverter::showProgress(uint64_t current, uint64_t total, int barWidth) {
    if (total == 0) return;
    
    float percent = (float)current / total;
    int pos = barWidth * percent;
    
    std::cout << "\r[";
    for (int i = 0; i < barWidth; i++) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    
    std::cout << "] " << std::setw(3) << int(percent * 100) << "% "
                << "(" << current / 1024 / 1024 << "/" << total / 1024 / 1024 << " MB)";
    std::cout.flush();
    
    if (current == total) {
        std::cout << std::endl;
    }
}

//Функция пропуска указанного кол-ва байтов байт
void MpdDataConverter::SkipBytes(uint64_t bytesToSkip) {
    if (bytesToSkip == 0) return;
    if(bytesToSkip>g_LenEvent){
        std::cout<<"Error! After the skip "<<bytesToSkip<<" byte we go beyond the event!"<<std::endl;
        return;
    }
    inFile.seekg(bytesToSkip, std::ios::cur);
    g_CountReadFileByte += bytesToSkip;
    //Это сделано, если идет небольшой скип внутри ивента, не выходя за его границу. (или канала, или adc)
    if(g_CountReadEventByte+bytesToSkip<=g_LenEvent)g_CountReadEventByte += bytesToSkip;
    if(g_CountReadAdcByte+bytesToSkip<=g_LenAdc)g_CountReadAdcByte += bytesToSkip;
    if(g_CountReadChannelByte+bytesToSkip<=g_LenChannel)g_CountReadChannelByte += bytesToSkip;
    if(g_CountReadChSamplesByte+bytesToSkip<=g_LenChannelSamples)g_CountReadChSamplesByte += bytesToSkip;

    if (inFile.fail()) {
        throw std::runtime_error("Failed to skip " + std::to_string(bytesToSkip) + " bytes");
    }
}

//Функция чтения по 4 байта
uint32_t MpdDataConverter::readWord(TypeReadByte type) {
    uint32_t _word = 0;
    if (!inFile.is_open()) throw std::runtime_error("File is not open");
    if (g_CountReadFileByte + 4 > g_LenFile) throw std::runtime_error("Not enough bytes in file");
    inFile.read(reinterpret_cast<char*>(&_word), sizeof(_word));
    if (inFile.gcount() != sizeof(_word)) throw std::runtime_error("Failed to read 4 bytes from file");
    
    g_CountReadFileByte+=4;
    switch (type) {
        case TypeReadByte::event: 
            g_CountReadEventByte += 4; 
            break;
        case TypeReadByte::adc :
            g_CountReadEventByte += 4; 
            g_CountReadAdcByte += 4;
            break;
        case TypeReadByte::channel: 
            g_CountReadEventByte += 4; 
            g_CountReadAdcByte += 4; 
            g_CountReadChannelByte += 4; 
            break;
        case TypeReadByte::chSamples:
            g_CountReadEventByte += 4; 
            g_CountReadAdcByte += 4; 
            g_CountReadChannelByte += 4; 
            g_CountReadChSamplesByte += 4; 
            break;
        case TypeReadByte::unknown:
            break;
    }

    if (g_CountReadFileByte % (1024 * 1024) < 4) {
        ProgressBar(g_CountReadFileByte, g_LenFile);
    }

    return _word;
}

bool MpdDataConverter::ReadEvent() {
    if (g_LenFile == 0) return false;
    if (g_CountReadFileByte >= g_LenFile) return false;
    
    uint32_t _word = 0;
    g_CountReadEventByte=0;

    // Поиск синхрослова
    while (g_CountReadFileByte < g_LenFile) {
        _word = readWord(TypeReadByte::unknown);
        if (_word == SYNC_WORD) {
            CountEventNumber+=1;
            break;
        }
    }

    if (_word != SYNC_WORD) {
        return false;
    }

    g_LenEvent = readWord(TypeReadByte::unknown);
    EventNumber = readWord(TypeReadByte::event);
    //SkipBytes(g_LenEvent-g_CountReadEventByte);

    return true;
}

bool MpdDataConverter::ReadADC(){
    if (g_LenFile == 0) return false;
    if (g_CountReadFileByte >= g_LenFile) return false;
    if (g_CountReadEventByte>=g_LenEvent)return false;
    
    uint32_t _word = 0;
    g_CountReadAdcByte=0;

    // [i+3] Device Serial
    deviceSerial=readWord(TypeReadByte::event);
    // [i+4] Device ID and Len
    devIdAndLen=readWord(TypeReadByte::event);
    g_LenAdc = devIdAndLen & 0xFFFFFF;
    
    // [i+5] (subtype and words)
    unclearUsage=readWord(TypeReadByte::adc);
    // [i+6]
    deviceTimeStampSec=readWord(TypeReadByte::adc);
    // [i+7]
    deviceTimeStampSecNs=readWord(TypeReadByte::adc);
    // [i+8]
    unclearUsage=readWord(TypeReadByte::adc);
    // [i+9]
    unclearUsage=readWord(TypeReadByte::adc);
    
    //SkipBytes(g_LenAdc-g_CountReadAdcByte);
    
    return true;
};

bool MpdDataConverter::ReadChannel(){
    if (g_LenFile == 0) return false;
    if (g_CountReadFileByte >= g_LenFile) return false;
    if (g_CountReadAdcByte>=g_LenAdc) return false;

    uint32_t _word=0;
    g_CountReadChannelByte=0;
    channelEvent.Clear();
    channelEvent.adcValues.reserve(60);

    // [i+10] Ch_{num}-1 and Len byte channel
    ChannelNumAndLen=readWord(TypeReadByte::adc);
    g_LenChannel = (ChannelNumAndLen & 0xFF) & 0xFC;
    uint32_t channel_raw = ((ChannelNumAndLen >> 24) & 0xFF) + 1;
    // [i+11]
    channelEvent.timestamp_sec=readWord(TypeReadByte::channel);
    // [i+12]
    channelEvent.timestamp_ns=readWord(TypeReadByte::channel);
    // [i+13;i+42] channel samples
    g_CountReadChSamplesByte=0;
    g_LenChannelSamples = g_LenChannel - 8;
    while(g_CountReadChSamplesByte<g_LenChannelSamples){
        _word=readWord(TypeReadByte::chSamples);
        //PrintCurrentCounByte();
        channelEvent.adcValues.push_back( ChangeWfValue ( AdcSampleValueFirst(_word) ) );
        channelEvent.adcValues.push_back( ChangeWfValue ( AdcSampleValueSecond(_word) ) );
    }

    channelEvent.eventNum = EventNumber;
    channelEvent.deviceID = deviceSerial;
    channelEvent.channelNums = SetChannelNum(getAddressForBasket_38(deviceSerial), channel_raw);
    channelEvent.channel_Phi = GetPhi(channelEvent.channelNums);
    channelEvent.channel_Z = GetZ(channelEvent.channelNums);
    channelEvent.timestamp_f = SetTimeStampNs(channelEvent.timestamp_sec,channelEvent.timestamp_ns);    
    if(channelEvent.adcValues.size()>0){
        channelEvent.numSamples=channelEvent.adcValues.size();
    }

    //SkipBytes(g_LenChannel-g_CountReadChannelByte);

    return true;
};

Int_t MpdDataConverter::AdcSampleValueFirst(uint32_t _adcValue){
    uint16_t raw = (_adcValue >> 16) & 0xFFFF;
    return static_cast<Int_t>(static_cast<int16_t>(raw));
}

Int_t MpdDataConverter::AdcSampleValueSecond(uint32_t _adcValue){
    uint16_t raw = _adcValue & 0xFFFF;
    return static_cast<Int_t>(static_cast<int16_t>(raw));
}

Int_t MpdDataConverter::ChangeWfValue(Int_t _value){
    switch (g_TypeWF) {
        case TypeWaveForms::raw : return _value;
        case TypeWaveForms::baseline : return (_value - wf_add_val);
        case TypeWaveForms::invert : return (wf_scale * _value) + wf_add_val;
    }
    return _value;
}

ULong64_t MpdDataConverter::SetTimeStampNs(uint32_t _sec, uint32_t _ns){
    uint32_t nanoseconds = (_ns & 0xFFFFFFFC) >> 2;
    return static_cast<ULong64_t>(_sec) * 1000000000ULL + nanoseconds;
}

Short_t MpdDataConverter::SetChannelNum(uint32_t _AdcId, uint32_t _chNum){    
    if(_AdcId==0)return -1;
    return (_AdcId-1)*64 + _chNum;
}