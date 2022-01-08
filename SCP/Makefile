TARGET = CMSA
CXXFLAGS = -ansi -O3
OBJS = Timer.o #Random.o

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

CPLEXDIR      = /home/marko/Desktop/CPLEX_Studio127/cplex
CONCERTDIR    = /home/marko/Desktop/CPLEX_Studio127/concert
CPDIR    = /home/marko/Desktop/CPLEX_Studio127/cpoptimizer

CCC = g++ -std=c++17  #g++
CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBCP = $(CPDIR)/lib/x86-64_sles10_4.1/static_pic

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -pthread
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CPINDIR = $(CPDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(CPINDIR)

all: ${TARGET}

CMSA: CMSA.o $(OBJS)
	$(CCC) $(CCFLAGS) CMSA.o $(OBJS) -o CMSA $(CCLNFLAGS)

CMSA.o: cplex_woagp.cpp
	$(CCC) -c $(CCFLAGS) CMSA.cpp -o CMSA.o 

clean:
	@rm -f *~ *.o ${TARGET} core
