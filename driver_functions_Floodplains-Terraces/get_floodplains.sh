#!/bin/sh
mkdir -p build_floodplains/
cd build_floodplains/
cp ../CMakeLists_get_floodplains.txt CMakeLists.txt
cmake .
make
mv get_floodplains.out ../
cd ..
