
## University of Tartu

# Gromov-Wasserstein Estimator


## Installation

Basically, it works only on Unix systems.  You need g++, in a version which supports C++14.

1. Install the Gurobi optimizer http://www.gurobi.com/
2. Edit the Makefile, to provide the paths to the Gurobi library and include files
3. Create an empty `Makefile.depend`
4. `make depend`
5. `make'
6. ./estimate-gromov-wasserstein -h

That should do it.


DOT
