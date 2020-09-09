#!/bin/bash
for i in $(seq 8 2 200)
do
   ./greedy.out -f preprocessing/small/min-$i-1_D1.txt -t 10 -greedy 3 -l result/result_D1_small_greedy3_w-type1.txt
done
