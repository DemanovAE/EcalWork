cmake_minimum_required(VERSION 3.10)
project(MPDDataConverter)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set nlohmann_json path BEFORE finding ROOT
set(nlohmann_json_DIR /cvmfs/nica.jinr.ru/sw/202309/slc7_x86-64/nlohmann_json/v3.11.3-2/share/cmake/nlohmann_json)
# Also add to CMAKE_PREFIX_PATH as a backup
list(APPEND CMAKE_PREFIX_PATH /cvmfs/nica.jinr.ru/sw/202309/slc7_x86-64/nlohmann_json/v3.11.3-2)

# Найти ROOT
find_package(ROOT REQUIRED)

# Создаем динамическую библиотеку
add_library(MPDDataConverter SHARED
    MpdDataConverter.cpp
    EcalDrawClass.cpp
)

# Подключаем ROOT
target_include_directories(MPDDataConverter PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(MPDDataConverter ${ROOT_LIBRARIES})

# Устанавливаем библиотеку (опционально)
install(TARGETS MPDDataConverter DESTINATION lib)
install(FILES MpdDataConverter.h DESTINATION include)

message(STATUS "ROOT found: ${ROOT_VERSION}")
message(STATUS "Library will be built as: libMPDDataConverter.so")