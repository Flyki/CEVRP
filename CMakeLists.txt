cmake_minimum_required(VERSION 3.15)
project(CEVRP_Yahui)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3") # -O3 optimization argument

include_directories(.)

add_executable(CEVRP_Yahui
        ant.h
        BACO2.cpp
        BACO2.h
        CACO.cpp
        CACO.h
        case.cpp
        case.h
        struct.h
        utilities.cpp
        utilities.h
        main.cpp)
