#!/bin/sh
cd BUILD
cmake -G "Unix Makefiles" ..
make ProGraphMSA
cd ..
cp BUILD/src/ProGraphMSA .
