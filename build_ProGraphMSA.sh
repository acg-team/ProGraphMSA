#!/bin/bash
# Example for how to build ProGraphMSA

CMAKE_OPTIONS=(-G "Unix Makefiles")

# Some customizations for homebrew
if [ -x /usr/local/include/eigen3 ]; then
    CMAKE_OPTIONS+=(-D EIGEN_INCLUDE_DIR=/usr/local/include/eigen3)
fi

if [ -x /usr/local/bin/c++-7 ]; then
    CMAKE_OPTIONS+=(-D CMAKE_CXX_COMPILER=/usr/local/bin/c++-7)
elif [ -x /usr/local/bin/c++-8 ]; then
    CMAKE_OPTIONS+=(-D CMAKE_CXX_COMPILER=/usr/local/bin/c++-8)
fi

if [ -x /usr/local/bin/gcc-7 ]; then
    CMAKE_OPTIONS+=(-D CMAKE_C_COMPILER=/usr/local/bin/gcc-7)
elif [ -x /usr/local/bin/gcc-8 ]; then
    CMAKE_OPTIONS+=(-D CMAKE_C_COMPILER=/usr/local/bin/gcc-8)
fi

# After changing cmake options, clean the BUILD directory
#rm -rf BUILD

mkdir -p BUILD &&
cd BUILD &&
cmake "${CMAKE_OPTIONS[@]}" .. &&
make ProGraphMSA &&
cp src/ProGraphMSA ../