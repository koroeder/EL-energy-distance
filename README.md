# EL-energy-distance

Script to obtain the energy distance between all pairs of minima for a energy landscape represented as database of minima and transition states.
The energy distance for two minima is the sum of the uphill energy barriers on the shortest discrete path between the minima.
The program uses Dijkstra's shortest path algorithm.

When graphs are created form the database, we ignore multiple edges and use the one with the lowest transition state energy.
Transition states connecting to a minimum with itself are also ignored.
To create directional edges, each TS is added as two edges with a weight given by the energy barrier.

**Compilation and example**
For compilation, run the standard compilation for fortran on your chosen compiler 
(e.g. f95 getEdist.f90 -o CalEdist). Compilation has been tested on gfortran and ifort, though individual compiler version will not work (generally due to compiler issues).

There are two required arguments and one optional argument. The required arguments are the number of minima and number of transition states.
There is a parser for those files implemented, but it seg faults with the inquire function on some clusters, so it is commented out for now.

The optional argument is the starting point for the distance calculations. The program uses Dijkstra's search from a single source to all other nodes.
It iterates over all nodes to complete this search for every minimum. The distance from src to all minima is written out after the search is complete,
and so the calculatiosn can be restarted at any point. The optional argument gives the index of the first minimum for the iteration.

There is a test for the connectivity of the database, which is currently not used.

The example can be run with CalcEdist 11 11 > out.log and should yield the output provided in expected_out.
