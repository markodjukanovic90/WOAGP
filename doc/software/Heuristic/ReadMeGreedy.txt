Kompilacija: g++ greedy.cpp -o greedy 
Poziv: ./greedy -f ../milan/preprocessing/small/min-100-1_D1.txt -alg 1 -greedy 3 -turn_ls 1 -w_type 3 -partial

-f: path to instance file;
-t: amount of time to allow execution;

-alg: algoritam koji pozivamo:
      0: pure Greedy
      1: Greedy+LS;
      2: Greedy + CPLEX

-greedy: type of the greedy method
         0: Lovasz
         3: Dragan's greedy 

-turn_ls: ukljuciti / iskljuciti LS
          0: Greedy + remove Vertices
          1: Greedy + LS + remove vertices

-partial: ( partial * n )- solution of a greedy will be 
          included into CPLEX 
          if partial == 0, then run pure CPLEX model 
 
w_type: type of weight involved into benchmarks 
        0: tezina proporcionalna sa velicinom skupa S[i]
        1: srednja vrijednost duzine ivica koje idu iz tjemena 
        2: avg visibility tezina
        3: random 
        4: non-weighted verzija

dodatni parametri za greedy + shaking:

-dist:  int, distanca od cvora i i ostalih cvorova u parcijalnom rjesenju --> za kandidate za izbacivanje se uzmu samo oni cvorovi koji su na udaljenosti vecoj od <dist> od cvora i
-offset: float (= p = offset * Parcijalno_rjesenje): Ako je cvor v_i u parcijalnom rjesenju, 
               za shaking razmatramo samo cvorove {v_{i-p},...., v_i,... v_{i+p}}
