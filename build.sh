#!/bin/sh

#  build.sh
#  
#
#  Created by Ollie Steiner on 30/07/17.
#

mkdir build
cd build
cmake ..
make
./UnscentedKF
