#!/bin/bash
for k in $(seq 8 1 8)
do
   for j in $(seq 0 1 3)
   do
      for i in $(seq 8 2 200)
      do
         ./greedy.out -f ../milan/preprocessing/large/fat-$i-1_D1.txt -t 10 -greedy $j -w_type $k -l ../milan/result/result_D1_large-W-type$k-Greedy$j.txt
      done
   done
done
