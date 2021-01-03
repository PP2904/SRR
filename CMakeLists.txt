cmake_minimum_required(VERSION 3.16)
project(SRR)

set(CMAKE_CXX_STANDARD 14)

#3.1.21: mainPR_D.cpp hinzugefügt => Fehler: ld: 3 duplicate symbols for architecture x86_64
add_executable(SRR main.cpp)