#!/bin/bash
#No esta funcionando {0..9}
for i in 0 1 2 3 4 5 6 7 8 9
do
	for j in 1 2 3 4 5
	do
		./secuencial/EC_2D < ./in/in$i >> ./times/secuencial/times$i.txt
        #./paralelo/build/EC_Par < ./in/in$i >> ./times/paralelo/times$i.txt
	done
	echo "Tiempos finalizados con "in$i
done
