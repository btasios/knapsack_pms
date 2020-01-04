#!/bin/bash
echo "Running the Generator"
for i in 1 2 3 4 5; do
	for n in {10,50,100,500}; do
		for r in 50 100 500 1000; do
			for t in 1 2 3 4; do
				./generator $n $r $t $i 5
				mv test.in problem_${n}_${r}_${t}_${i}_5.txt
			done
		done
	done
done
