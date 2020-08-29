\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 obj: 4 x1 + 8 x2 + 6 x3 + 4 x4 + 7 x5 + x6
Subject To
 c1: x1 + x2 + x3 + x4 >= 1
 c2: x1 + x3 + x4 + x5 >= 1
 c3: x1 + x2 + x4 >= 1
 c4: x1 + x4 >= 1
 c5: x2 + x4 >= 1
Bounds
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
      x6 = 0
Binaries
 x1  x2  x3  x4  x5 
End
