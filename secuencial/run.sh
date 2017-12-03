#!/bin/bash
let "N  = $1 + $1"
g++ ecuacionCalor2D.cpp -o t -O2 -larmadillo
./t < in > ./out/out
head -n $1 ./out/out > ./out/out0
head -n $N ./out/out | tail -n $1 > ./out/out1
tail -n $1 ./out/out > ./out/out2
gnuplot plot