# MPD Data Converter / ECalWork

Библиотека и утилиты для:
1. конвертации бинарных файлов калориметра MPD (`.data`) в ROOT-файлы с TTree `events`;
2. анализа этих ROOT-файлов (продольный / поперечный анализ ECAL) локально и на кластере (SLURM job array).

Library and tools for:
1. converting MPD calorimeter binary `.data` files to ROOT files with an `events` TTree;
2. running longitudinal / transverse ECAL analysis on these ROOT files, locally or on a SLURM cluster (job arrays).

---

## 1. Сборка / Build

```bash
# Клонирование репозитория / Clone repository
git clone https://github.com/DemanovAE/EcalWork.git
cd EcalWork

# Сборка CMake / Build with CMake
mkdir build && cd build
cmake ..
make
```

После сборки будут созданы (имена могут отличаться, но логика такая):
- основная библиотека преобразования данных: `build/libMPDDataConverter.so` (или `.dylib` на macOS);
- исполнительный файл анализа: `build/EcalWorkExec`.

After building you get (names may differ slightly, but conceptually):
- core data conversion library: `build/libMPDDataConverter.so` (or `.dylib` on macOS);
- analysis executable: `build/EcalWorkExec`.

---

## 2. Шаг 1: Конвертация `.data` → ROOT / Step 1: `.data` → ROOT

Конвертация выполняется функцией:

```cpp
void ConvertToRoot(int targetEvent = -1,
                   std::vector<int> targetEvents = {},
                   std::string inputData  = "../run_rc-hs1_088.data",
                   std::string outputData = "out_all.root");
```

Типичный запуск в ROOT (интерактивный узел, без SLURM):

```bash
# Интерактивный ROOT
root -l EcalWork.cpp
```

Внутри ROOT:

```cpp
.L EcalWork.cpp+

// Простой пример: конвертация всего файла в один ROOT TTree "events"
ConvertToRoot(-1, {}, "/path/to/run_rc-hs1_088.data", "run_rc-hs1_088.root");
```

Что делает `ConvertToRoot`:
- читает бинарный `.data` файл, используя `MpdDataConverter`;
- применяет настройки pedestal и обработку пиков/интегралов (см. код в `EcalWork.cpp`);
- создаёт ROOT-файл `outputData` с TTree `events` и нужными полями;
- при необходимости может сохранять waveform-гистограммы и картинки (каталог `pict/` рядом с ROOT-файлом).

What `ConvertToRoot` does:
- reads the raw `.data` file via `MpdDataConverter`;
- applies pedestal and ADC waveform processing (see `EcalWork.cpp`);
- produces a ROOT file `outputData` with an `events` TTree;
- optionally saves waveform histograms and plots (under `pict/` next to the ROOT file).

---

## 3. Шаг 2: Анализ ROOT TTree / Step 2: ROOT TTree analysis

Основная функция анализа:

```cpp
void EcalWork(std::string inputDataTree =
                  "/nica/mpd1/demanov/ecal_mpd/run_rc1-hs4_133_basket5_1616162.root",
              std::string outputData = "hs4_133_basket5_1616162.root",
              Long64_t firstEntry    = 0,
              Long64_t lastEntry     = 20.e6);
```

Функция:
- открывает входной ROOT-файл `inputDataTree`;
- читает TTree `"events"` (обязательно наличие дерева с таким именем);
- использует набор ветвей:
  - `eventNum`
  - `integral`
  - `amplitude`
  - `channelNums`
  - `channel_Z`
  - `channel_Phi`
- строит события по `eventNum`, выполняет либо **продольный** (`LongAnalysis`) либо **поперечный** (`TransverseAnalysis`) отбор;
- заполняет сводные гистограммы (амплитуды, интегралы, распределения по каналам);
- пишет выходной ROOT-файл `outputData`:
  - основные гистограммы (event, channel, amplitude, integral);
  - директорию `target_channels` с гистограммами по выбранным «целевым» каналам;
  - директорию `channel_wf` при необходимости (waveforms).

The analysis function:
- opens the input ROOT file `inputDataTree`;
- reads the `events` TTree (must exist in the file);
- uses branches `eventNum`, `integral`, `amplitude`, `channelNums`, `channel_Z`, `channel_Phi`;
- builds events by `eventNum`, applies **longitudinal** or **transverse** selection;
- fills summary histograms and per-channel target histograms;
- writes an output ROOT file `outputData` with analysis results.

---

## 4. Локальный запуск `EcalWorkExec` / Local run of `EcalWorkExec`

Исполнительный файл `EcalWorkExec` — это обёртка над `EcalWork(...)`, которая принимает аргументы командной строки:

```bash
./EcalWorkExec <input.root> <output.root> <firstEntry> <lastEntry>
```

Пример локального запуска на интерактивном узле:

```bash
cd build

./EcalWorkExec \
  /nica/mpd1/demanov/ecal_mpd/run_rc1-hs4_133_basket5_1616162.root \
  /scratch3/dflusova/ECalWork/root/test_small.root \
  0 \
  10000000
```

Это запустит `EcalWork` на событиях с индексами от `firstEntry` до `lastEntry` (верхняя граница автоматически ограничивается числом записей в TTree `events` внутри кода).

This runs `EcalWork` on entries in the range `[firstEntry, lastEntry)` (upper bound is clamped to `tree->GetEntries()` inside the code).

---

## 5. Запуск на кластере (SLURM job array) / Cluster run (SLURM job array)

Для массового анализа большого ROOT-файла используется скрипт `run_script/run_EcalWork.sh`, который запускает `EcalWorkExec` как SLURM job array.

### 5.1. Скрипт `run_EcalWork.sh`

Основные элементы скрипта:

- SLURM директивы:

```bash
#SBATCH -p nica
#SBATCH -J ECAL
#SBATCH -a 0-9                 # 10 задач массива: TASK_ID = 0..9
#SBATCH -N 1
#SBATCH -o /scratch3/dflusova/ECalWork/log/qa_%A_%a.out
#SBATCH -e /scratch3/dflusova/ECalWork/log/qa_%A_%a.err
```

- Инициализация окружения MPDROOT:

```bash
source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev/v24.09.24-1
source /lhep/users/dflusova/mpdroot_polarization/mpdroot/install/config/env.sh
```

- Пути:

```bash
export BUILD_DIR=/lhep/users/dflusova/EcalWork/EcalWork/build
export EXEC=${BUILD_DIR}/EcalWorkExec

export INPUT_DIR=/nica/mpd1/demanov/ecal_mpd
export INPUT_FILE=${INPUT_DIR}/run_rc1-hs4_133_basket5_1616162.root

export OUT_BASE=/scratch3/dflusova/ECalWork
export OUT_DIR=${OUT_BASE}/root
export LOG_DIR=${OUT_BASE}/log

export OUT_FILE=${OUT_DIR}/run_rc1-hs4_133_basket5_1616162_part${TASK_ID}.root
export LOG=${LOG_DIR}/log_qa_run8_${JOB_ID}_${TASK_ID}.log
```

- Разбиение по диапазонам записей (chunking):

```bash
# total entries in TTree "events" (задаётся вручную)
TOTAL_ENTRIES=584947682

N_CHUNKS=$(( SLURM_ARRAY_TASK_MAX - SLURM_ARRAY_TASK_MIN + 1 ))
CHUNK_SIZE=$(( TOTAL_ENTRIES / N_CHUNKS ))

FIRST_ENTRY=$(( TASK_ID * CHUNK_SIZE ))
if [ "${TASK_ID}" -eq "$((N_CHUNKS-1))" ]; then
  LAST_ENTRY=${TOTAL_ENTRIES}
else
  LAST_ENTRY=$(( (TASK_ID + 1) * CHUNK_SIZE ))
fi
```

Каждая задача массива (`TASK_ID = 0..9`) получает свой диапазон `[FIRST_ENTRY, LAST_ENTRY)` и пишет отдельный выходной ROOT-файл.

Each array task (`TASK_ID = 0..9`) processes its own entry range `[FIRST_ENTRY, LAST_ENTRY)` and writes a separate output ROOT file.

- Запуск:

```bash
cd "${BUILD_DIR}"

"${EXEC}" "${INPUT_FILE}" "${OUT_FILE}" ${FIRST_ENTRY} ${LAST_ENTRY} &>> "${LOG}"

echo "Job is finished" &>> "${LOG}"
```

### 5.2. Подготовка и отправка / Submission

1. Убедитесь, что:
   - `EcalWorkExec` собран и доступен в `${BUILD_DIR}`;
   - путь к входному файлу `INPUT_FILE` правильный;
   - значение `TOTAL_ENTRIES` соответствует `events->GetEntries()` для данного входного файла (его нужно обновлять вручную при смене файла).

2. Отправка job array:

```bash
cd run_script
sbatch run_EcalWork.sh
```

SLURM создаст 10 задач, каждая обработает свою часть TTree и положит результат в:
`$OUT_DIR/run_rc1-hs4_133_basket5_1616162_part<TASK_ID>.root`.

---

## 6. Объединение выходных файлов / Merging output ROOT files

После завершения job array вы получите набор файлов:

```text
run_rc1-hs4_133_basket5_1616162_part0.root
run_rc1-hs4_133_basket5_1616162_part1.root
...
run_rc1-hs4_133_basket5_1616162_part9.root
```

Возможные стратегии:
- анализ каждого файла отдельно (например, по задачам/чанкам);
- объединение гистограмм через ROOT (`hadd` или ручной merge).

Пример с `hadd`:

```bash
hadd -f run_rc1-hs4_133_basket5_1616162_all.root \
     run_rc1-hs4_133_basket5_1616162_part*.root
```

Или более аккуратный merge через собственный ROOT-скрипт, если структура выходных файлов сложнее (несколько директорий, дополнительные TTree и т.п.).

---

## 7. Режимы анализа / Analysis modes

Переключение между **продольным** и **поперечным** анализом контролируется глобальной логикой в `EcalWork.cpp` (флаг `TransverAnalysis` и соответствующие функции `LongAnalysis` / `TransverseAnalysis`).

Основные принципы:
- Longitudinal mode:
  - поиск самой «горячей» ячейки и её ближайшего соседа (3×3 окно);
  - расчёт энергии в 5×5 окне вокруг hottest cell;
  - различие между «core» и «diffuse» энергией (ratio);
  - строгие cuts по интегралу, положениям по Z/Phi, отношению энергий.

- Transverse mode:
  - поиск стрипов по Z или Phi с последовательными ячейками;
  - отбор событий с достаточной длиной стрипа и малой контаминацией;
  - выбор целевой ячейки в стрипе и заполнение гистограмм.

Все детали реализованы в `LongAnalysis(...)` и `TransverseAnalysis(...)` в `EcalWork.cpp`. При модификации критериев отбора или порогов (например, `gl_long_max_1st_Integral`, `gl_long_max_3x3_ratio`) эти значения нужно менять в исходнике.

---

## 8. Требования / Requirements

- ROOT 6.x
- C++17
- MPDROOT environment (модули и env-файл, как показано в `run_EcalWork.sh`)
- SLURM для запуска на кластере (если используете job array)

---

## 9. Краткое резюме / Short summary

1. Сначала конвертируете бинарный `.data` → ROOT TTree `events` с `ConvertToRoot(...)`.
2. Затем анализируете TTree с `EcalWork(...)`:
   - локально: `./EcalWorkExec input.root output.root firstEntry lastEntry`;
   - на кластере: `sbatch run_script/run_EcalWork.sh` с job array и разбиением на чанки.
3. Объединяете или анализируете выходные ROOT-файлы по необходимости.

First, convert `.data` to a ROOT `events` TTree with `ConvertToRoot(...)`, then run `EcalWork(...)` (via `EcalWorkExec`) locally or as a SLURM job array, and finally inspect or merge the per-chunk output ROOT files as needed.