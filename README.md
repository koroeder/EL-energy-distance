# EL-energy-distance

Scripts to obtain the energy distance between all pairs of minima for a energy landscape represented as database of minima and transition states.

## Fortran script

### Performance
This program is tuned for fast performance. The scripts have not been benchmarked, but some data from our testing:
Fortran script on around 40,000 minima (1.6 billion distances) took around 24 h on a single CPU.
The python scripts on the same CPU managed around 800 out of 7000 minima in 24 h (5.6 million distances).
The speed up is around 285 times , though there are algorithmic differences.

The Fortran script also uses less memory. The python scripts require roughly 10Mb per 250,000 minima and scaling as the square of the number of minima.
This translates to around 4GB for 10,000 minima and 64GB for 40,000 minima. The Fortran version stores a couple of MB for 10,000 minima.
This is mainly related to arithmetic differences (we use memorisation in the python scripts and a dictionary of dictionary to have a hashable data set to minimise the number of necessary distance calculations).

### Basic algorithm
The energy distance for two minima is the sum of the uphill energy barriers on the shortest discrete path between the minima.
The program uses Dijkstra's shortest path algorithm.

When graphs are created form the database, we ignore multiple edges and use the one with the lowest transition state energy.
Transition states connecting to a minimum with itself are also ignored.
To create directional edges, each TS is added as two edges with a weight given by the energy barrier.

### Compilation and example
For compilation, run the standard compilation for fortran on your chosen compiler 
(e.g. f95 getEdist.f90 -o CalEdist). Compilation has been tested on gfortran and ifort, though individual compiler version will not work (generally due to compiler issues).

There are two required arguments and one optional argument. The required arguments are the number of minima and number of transition states.
There is a parser for those files implemented, but it seg faults with the inquire function on some clusters, so it is commented out for now.

The optional argument is the starting point for the distance calculations. The program uses Dijkstra's search from a single source to all other nodes.
It iterates over all nodes to complete this search for every minimum. The distance from src to all minima is written out after the search is complete,
and so the calculatiosn can be restarted at any point. The optional argument gives the index of the first minimum for the iteration.

There is a test for the connectivity of the database, which is currently not used.

The example can be run with 

> `CalcEdist 11 11 > out.log'

and should yield the output provided in expected_out.


## Python scripts
The python scripts provided are mainly for conceptual purposes. 
There are two sets of scripts:
- get_Edist.py
- get_ELgraph.py and calc_Edist.py

The first one parses the database, creates the graph and then runs the calculations.
The second set separates setup and calculations. get_ELgraph.py sets up the database and saves them (via pickling).
Distances are then calculated with calc_Edist.py. Backup files are created here, so the scripts can be restarted.

### Basic algorithm

The python scripts use networkx and numpy. The use of networkx leads to poor performance (igraph or networkit might give better performance),
but it seems networkx is (currently) easier to set up.

Graphs are created in networkx and multi edges are avoided the same way they are in the Fortran script.
We obtain the paths from dijkstra searches for pairs of nodes (the runtime for allpaths from single source was too slow to be effective for testing, but could be used). All paths contained in the obtained shortest paths must be shortest paths themselves, so every calaculation yields a large number of shortest paths, which are stored in dictionaries.

Two distances are calculated - one using the sum of barriers, and another using the energy difference between the source minimum and the highest energy transition state on the path.

### Examples and their output
Expected_out1 uses the first script as 
> `python get_Edist.py > out.log'

Expected_out2 uses 
> `python get_ELgraph.py > out.graph.log'

followed by 

> `python ../../calc_Edist.py > out.dist.log'

The output files are the two distance matrix as csv files, and the pickeld and backup data for the second set of scripts.


