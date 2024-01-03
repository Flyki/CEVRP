cmake_minimum_required(VERSION 3.15)
project(CEVRP-Yahui)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3") # -O3 optimization argument

include_directories(.)

add_executable(CEVRP-Yahui
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
        main.cpp
        stats.h
        stats.cpp)

# Link with pthread
find_package(Threads REQUIRED)
target_link_libraries(CEVRP-Yahui PRIVATE Threads::Threads)