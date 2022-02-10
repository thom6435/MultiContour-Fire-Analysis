cmake_minimum_required(VERSION 3.17)
project(MultiContour)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(MultiContour
        contourregion.cpp
        contourregion.h
        fireline.cpp
        fireline.h
        imageanalysis.cpp
        imageanalysis.h
        mca_config.cpp
        mca_config.h
        MultiContour.cpp
        MultiContour.hpp
        radiometer.cpp
        radiometer.h
        normal.cpp
        normal.h thermalArea.cpp thermalArea.h)

target_link_libraries(MultiContour opencv_highgui opencv_core opencv_xphoto opencv_imgcodecs opencv_imgproc)