#!/bin/bash
g++ ./secuencial/ecuacionCalor2D.cpp -o EC_Sec -O2 -larmadillo

#cmake .. dentro de paralelo/build
#luego make para crear el ejecutable

for i in {0..4}
do
	for j in {1..5}
	do
		./EC_Sec < ./in/in$i >> ./times/secuencial/times$i.txt
        #./paralelo/build/EC_Par < ./in/in$i >> ./times/paralelo/times$i.txt
	done
	echo "Tiempos finalizados con "in$i
done