#!/bin/sh
#mkdir build
#cp *.cl build/
#cp *.conf build/
cd build
cmake -G "MSYS Makefiles" -D CMAKE_BUILD_TYPE=Debug ..
make
