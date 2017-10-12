#!/bin/sh
cd build_shapefile/
cmake .
make
mv get_terraces_from_shapefile.out ../
cd ..
