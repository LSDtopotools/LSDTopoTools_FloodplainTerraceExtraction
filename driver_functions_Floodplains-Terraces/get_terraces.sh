#!/bin/sh
cd build/
cmake .
make
mv get_terraces.out ../
cd ..
