Kompilacija: g++ greedy.cpp -o greedy 
Poziv: ./greedy -f instance/random/rand-84-19.pol -t 1 -greedy 0 -ls 1 -w_type 4

-f: path to instance file;
-t: amount of time to allow execution; 
-greedy: type of the greedy method
         0: Lovasz
         3: Dragan's greedy 
-ls: exclude/include LS into greedy 
     0: false
     1: true

w_type: type of weight involved into benchmarks 
        0: tezina proporcionalna sa velicinom skupa S[i]
        1: srednja vrijednost duzine ivica koje idu iz tjemena 
        2: avg visibility tezina
        3: random 
        4: non-weighted verzija
