cmake_minimum_required(VERSION 3.24)
project(AQuA_CPP)

set(CMAKE_CXX_STANDARD 11)


# include opencv library
set(OpenCV_DIR "C:/Program Files (x86)/opencv/mingw_build/install")
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})


# include FFTW library
set(FFTW_LIBS libfftw3-3 libfftw3f-3 libfftw3l-3)
link_directories(${CMAKE_SOURCE_DIR}/fftw3.3.5-64 )
include_directories(${CMAKE_SOURCE_DIR}/fftw3.3.5-64)


add_executable(${PROJECT_NAME}
                    src/data/data.h src/data/data.cpp
                    src/phaseRun.h src/getappdata.h src/feaRun.h src/evtRun.h src/actRun.h
                    src/flow/flow.h  src/preProcess/preProcessRun.h
                    src/draft1.cpp src/draft.h )

target_link_libraries(${PROJECT_NAME} ${FFTW_LIBS} ${OpenCV_LIBS} )

