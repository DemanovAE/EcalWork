cmake_minimum_required(VERSION 3.10)
project(MPDDataConverter)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Включаем позиционно-независимый код (обязательно для .so)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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

# Выводим информацию
message(STATUS "ROOT found: ${ROOT_VERSION}")
message(STATUS "Library will be built as: libMPDDataConverter.so")