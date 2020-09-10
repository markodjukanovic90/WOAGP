#!/bin/bash
for k in $(seq 8 1 8)
do
   for j in $(seq 3 1 3)
   do
      for i in $(seq 8 2 200)
      do
         ./greedy.out -f preprocessing/small/min-$i-1_D1.txt -t 10 -greedy $j -w_type $k -l result/result_D1_small-W-type$k-Greedy$j.txt
      done
   done
done
