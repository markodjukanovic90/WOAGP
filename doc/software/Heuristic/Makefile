TARGET = greedy
CXXFLAGS = -ansi -O3
OBJS = #Timer.o #Random.o

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

CPLEXDIR      = /home/djukanovic/Desktop/projects/LCAPS_software/cplex-12.5/cplex
CONCERTDIR    = /home/djukanovic/Desktop/projects/LCAPS_software/cplex-12.5/concert
CPDIR    = /home/djukanovic/Desktop/projects/LCAPS_software/cplex-12.5/cpoptimizer

CCC = g++ -std=c++17  #g++
CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBCP = $(CPDIR)/lib/x86-64_sles10_4.1/static_pic

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBCP) -lcp -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -pthread
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CPINDIR = $(CPDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(CPINDIR)

all: ${TARGET}

greedy: greedy.o $(OBJS)
	$(CCC) $(CCFLAGS) greedy.o $(OBJS) -o greedy $(CCLNFLAGS)

greedy.o: greedy.cpp
	$(CCC) -c $(CCFLAGS) greedy.cpp -o greedy.o 

clean:
	@rm -f *~ *.o ${TARGET} core
