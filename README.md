
## University of Tartu

# Gromov-Wasserstein Estimator


## Installation

Basically, it works only on Unix systems.  You need g++, in a version which supports C++14.

1. Install the Gurobi optimizer http://www.gurobi.com/
2. Edit the `Makefile`, to provide the paths to the Gurobi library and include files
3. Create an empty `Makefile.depend`
4. `make depend`
5. `make`

To use the program, type something like
* `./estimate-gromov-wasserstein -d -Oresult.p facebook.mmspace twitter.mmspace`
which would compute the distance between Facebook and Twitter, both scaled to have diameter 1 (-d). The distance is manifested through a pair of random variables (X,Y), whose joint probability distribution is written to `result.p`. (Note that there's no space between -O and the file name.)

Try `./estimate-gromov-wasserstein -h` for full help text describing all options and the file formats everything.

That should do it.


DOT
