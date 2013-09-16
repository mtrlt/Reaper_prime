#!/bin/sh
mkdir build
cp *.cl build/
cp *.conf build/
cd build
cmake -G "Unix Makefiles" -D CMAKE_BUILD_TYPE=Release ..
make
