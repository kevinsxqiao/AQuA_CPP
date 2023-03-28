cmake_minimum_required(VERSION 3.24)
project(AQuA_CPP)

set(CMAKE_CXX_STANDARD 11)


# Include FFTW library
link_directories(${CMAKE_SOURCE_DIR}/fftw3.3.5-64)
include_directories(${CMAKE_SOURCE_DIR}/fftw3.3.5-64)

add_executable(${PROJECT_NAME}  src/draft3.cpp
        src/preProcess/preProcessRun.h src/phaseRun.h src/getappdata.h src/feaRun.h src/evtRun.h src/actRun.h
        src/flow/flow.h src/data/data.h src/preProcess/regCrossCorrelation.cpp)

target_link_libraries(${PROJECT_NAME} libfftw3-3 libfftw3f-3 libfftw3l-3)
