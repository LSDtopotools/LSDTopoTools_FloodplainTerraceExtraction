#!/bin/sh
mkdir -p build/
cd build/
cp ../CMakeLists_get_terraces.txt CMakeLists.txt
cmake .
make
mv get_terraces.out ../
cd ..
