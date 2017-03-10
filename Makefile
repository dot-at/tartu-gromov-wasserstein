INCLUDE_GUROBI=-I /opt/gurobi/linux64/include/
LINK_GUROBI=-L/opt/gurobi/linux64/lib/ -lgurobi_c++ -lgurobi70

CCC=g++
#CCCDEBUGFLAGS= -g -O0 -fno-inline -fno-eliminate-unused-debug-types
CCCDEBUGFLAGS= -O
CCCFLAGS= $(INCLUDE_GUROBI) $(CCCDEBUGFLAGS) -std=gnu++11  -Wall #-DGromov_Wasserstein_DEBUG

HEADERS        = gromov_wasserstein.hh memoli.hh
TEMPLATE_FILES = memoli.t.cc
SOURCES        = $(TEMPLATE_FILES) gromov_wasserstein.cc
OBJECTS        = $(SOURCES:.cc=.o)
EXECS          = estimate-gromov-wasserstein

all: $(OBJECTS) $(EXECS) # all_transitions

# lib: $(LIBRARY)
#
# obj: $(OBJECTS)
#
# pslistings: $(PSLISTINGS)

clean:
	rm $(OBJECTS) # $(PSLISTINGS) $(LIBRARY)

depend:
	makedepend -o.o -f Makefile.depend -- $(CCCFLAGS) -- $(SOURCES)


# E x e c u t a b l e s

estimate-gromov-wasserstein: $(OBJECTS) estimate-gromov-wasserstein.cc
	$(CCC) $(CCCFLAGS) estimate-gromov-wasserstein.cc $(OBJECTS) $(LINK_GUROBI) -o estimate-gromov-wasserstein


# ------------------------------------------------------------

%.o: %.cc
	$(CCC) $(CCCFLAGS) -c $<



#
# %.a: $(OBJECTS)
# 	ar -r $@ $(OBJECTS)



include Makefile.depend
# DO NOT DELETE
