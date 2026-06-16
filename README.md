# MPD Data Converter

Библиотека для чтения бинарных `.data` файлов с калориметра MPD и конвертации в ROOT деревья.

## Сборка библиотеки

```bash
# Клонирование репозитория
git clone https://github.com/DemanovAE/EcalWork.git
cd EcalWork

# Сборка
mkdir build && cd build
cmake ..
make

# Библиотека будет создана: build/libMPDDataConverter.so

# Интерактивный режим
root -l EcalWork.cpp

# Batch режим
root -l -b -q EcalWork.cpp

# С параметрами
root -l -b -q 'EcalWork.cpp("input.data", "output.root", -1)'
