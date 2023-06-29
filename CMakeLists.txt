cmake_minimum_required(VERSION 3.24)
project(AQuA_CPP)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-fopenmp")


# include opencv library
set(OpenCV_DIR "C:/Program Files (x86)/opencv/mingw_build/install")
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})


# include FFTW library
set(FFTW_DIR "C:/Program Files (x86)/fftw3.3.5-64")
set(FFTW_LIBS libfftw3-3 libfftw3f-3 libfftw3l-3)
include_directories(${FFTW_DIR})
link_directories(${FFTW_DIR})


#include MATLAB
set(MATLAB_ROOT "C:/Program Files/MATLAB/R2022b")
include_directories("${MATLAB_ROOT}/extern/include")
find_library(MAT_LIBRARY mat "${MATLAB_ROOT}/extern/lib/win64/mingw64")
find_library(MX_LIBRARY mx "${MATLAB_ROOT}/extern/lib/win64/mingw64")


#include Boost
set(BOOST_ROOT "C:/Program Files (x86)/boost_1_82_0/boost_1_82_0")
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


add_executable(${PROJECT_NAME}
                    src/data/data.h src/data/data.cpp
                    src/phaseRun.h src/getappdata.h src/feaRun.h src/evtRun.h src/actRun.h
                    src/flow/flow.h src/preProcess/preProcessRun.h
                    src/draft1.cpp src/draft.h
                    src/preProcess/baselineRemoveAndNoiseEstimation.cpp src/preProcess/regCrossCorrelation.cpp src/preProcess/preProcessRun.cpp)

target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS} ${FFTW_LIBS} ${MAT_LIBRARY} ${MX_LIBRARY} ${Boost_LIBRARIES})

