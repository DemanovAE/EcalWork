#include "EcalDrawClass.h"
#include "MpdDataConverter.h"
#include "Rtypes.h"

// Forward declare EcalWork from EcalWork.cpp
void EcalWork(std::string inputDataTree,
              std::string outputData,
              Long64_t firstEntry,
              Long64_t lastEntry);

int main(int argc, char** argv)
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
              << " input.root output.root [firstEntry] [lastEntry]\n";
    return 1;
  }

  std::string input  = argv[1];
  std::string output = argv[2];

  Long64_t firstEntry = 0;
  Long64_t lastEntry  = -1;

  if (argc >= 4) {
    firstEntry = std::stoll(argv[3]);
  }
  if (argc >= 5) {
    lastEntry = std::stoll(argv[4]);
  }

  EcalWork(input, output, firstEntry, lastEntry);
  return 0;
}