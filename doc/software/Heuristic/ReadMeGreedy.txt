Kompilacija: g++ greedy.cpp -o greedy 
Poziv: ./greedy -f ../milan/preprocessing/small/min-100-1_D1.txt -alg 1 -greedy 3 -turn_ls 1 -w_type 3

-f: path to instance file;
-t: amount of time to allow execution; 
-alg: algoritam koji pozovamo:
      0: pure Greedy
      1: Greedy+LS;

-greedy: type of the greedy method
         0: Lovasz
         3: Dragan's greedy 

-turn_ls: ukljuciti / iskljuciti LS
          0: Greedy + remove Vertices
          1: Greedy + LS + remove vertices

w_type: type of weight involved into benchmarks 
        0: tezina proporcionalna sa velicinom skupa S[i]
        1: srednja vrijednost duzine ivica koje idu iz tjemena 
        2: avg visibility tezina
        3: random 
        4: non-weighted verzija
