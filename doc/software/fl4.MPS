* ENCODING=ISO-8859-1
NAME          IloCplex
ROWS
 N  obj     
 G  c1      
 G  c2      
 G  c3      
 G  c4      
 G  c5      
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    x1        obj                            87
    x1        c2                              1
    x2        obj                            84
    x2        c3                              1
    x3        obj                            16
    x3        c3                              1
    x4        obj                            78
    x4        c4                              1
    MARK0001  'MARKER'                 'INTEND'
    x5        obj                             1
RHS
    rhs       c1                              1
    rhs       c2                              1
    rhs       c3                              1
    rhs       c4                              1
    rhs       c5                              1
BOUNDS
 UP bnd       x1                              1
 UP bnd       x2                              1
 UP bnd       x3                              1
 UP bnd       x4                              1
 FX bnd       x5                              0
ENDATA
