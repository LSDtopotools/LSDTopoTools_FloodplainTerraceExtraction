#!/bin/sh
cd build/
cmake .
make
mv terraces_swath_driver.out ../
cd ..
