cmake_minimum_required (VERSION 3.16)
project(final_project_2024)
set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
include_directories(include)

macro(macro1 arg)
  add_executable(dspSimMain src/dspSimMain.cpp src/dspSolver_${arg}.cpp src/dspLinkConfig.cpp src/dspJoint.cpp src/dspLink.cpp)
  target_include_directories(dspSimMain PUBLIC include)
  target_link_libraries (dspSimMain Eigen3::Eigen)
endmacro()

macro1(20190837)