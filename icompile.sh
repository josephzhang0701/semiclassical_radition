#!/usr/bin/env bash
# use clang++ (Apple Clang) 或者 brew 的 g++
# -std=c++11 确保 C++11 支持
clang++ -std=c++11 *.cpp \
  $(pkg-config --cflags --libs gsl) \
  -lm \
  -o irun.x

# from old authors
# icc -Wall *.cpp -lgsl -lgslcblas -o ./irun.x