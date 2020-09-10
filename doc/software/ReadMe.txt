Parametri programa: 

-f : sraza ka instanci 
-t: vrijeme izvrsavana CPLEX-a 
-alg: tip alg. koji pozivamo
      0: CPLEX
      1: CP

Primjer poziva:
./cplex_woagp -f milan/preprocessing/small/min-10-1_D1.txt -t 1000000 -alg 1 -w_type 0
