cmake_minimum_required(VERSION 3.15)

project(MoleculesDynamics VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets DataVisualization)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

add_executable(MoleculesDynamics
    WIN32
    main.cpp
    Source/MoleculesDynamics.cpp
    Source/MoleculesDynamics.h
    Source/mathFunc.hpp
)

target_link_libraries(MoleculesDynamics PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Qt6::DataVisualization
)

target_include_directories(MoleculesDynamics PRIVATE
    Source/
)