Genetic Algorithm for Weighted MaxIS
====================================

This program uses a parallel genetic algorithm to approximate solutions
to the Weighted Maximum Independent Set problem on large graphs. It
implements both single-threaded and multi-threaded versions of the
algorithm found in Beasley and Chu's [A Genetic Algorithm for the Set
Covering
Problem](http://www.sciencedirect.com/science/article/pii/037722179500159X).
I used the
[BHOSLIB](http://www.nlsde.buaa.edu.cn/~kexu/benchmarks/graph-benchmarks.htm)
MIS benchmark graphs for testing, which are non-weighted. The algorithm
generates optimal solutions for most of the 450 vertex graphs. Finding
optimal search parameters for the larger graphs is a work in progress,
but I suspect that 450 vertices is near the upper limit if an optimal
solution is required in a reasonable amount of time. Like all genetic
algorithms, the solver provides no guarantee that a given solution is
optimal.

This solver is designed mostly to explore the efficacy of genetic algorithms on
highly-constrained NP optimization problems. Emphasis was put on flexibility
over performance. It also allows the user to select from several genetic
operators or implement their own by inheriting from a virtual base class.

Compilation and Usage
---------------------

The solver must be built by a compiler with c++14 support and depends on boost.
This means at least clang-3.5 or gcc-5.

The solver will log the state of the population every few thousand iterations.
Object files will be put in the build directory, and the executable will appear
in the bin directory. You may need to edit the Makefile to point your compiler
to the boost header files on your system.

Run `bin/graph --help` for usage information. The `--file` option can be
passed positionally to the executable.

Algorithm
---------

Familiarity with genetic algorithms is assumed for this section.

The fitness of a cantidate solution is the sum of the vertex weights
contained in that solution. For the unweighted version, every vertex is
given a weight of 1. The solver uses the incremental replacement method
and always replaces the weakest memeber of the population. In doing so,
we ensure that the average fitness of the population increases
monotonically at the cost of a higher likelihood of premature
convergence. This is functionally equivalent to a very high level of
elitism.

The minimum set covering problem is closely related to the minimum
vertex covering problem on a graph, whose solution is the complement to
that graph's minimum independent set. It is therefore reasonable to
expect that approaches which are effective on minimum set covering
problems will also work well for finding the maximum independent set
(MIS). Beasely and Chu describe (as far as I can tell) two novel
improvements on a typical genetic algorithm specifically for the minimum
set covering. Both of these are implemented in this solver.

The first is the heuristic feasibility operator. Given a randomly
generated set, such as one that is created by recombinating or mutating
other sets, it is highly likely that some vertices are now over- or
under-covered. A traditional genetic algorithm would solve this by
implementing a penalty function, but such a function is difficult to get
right, and results in yet another opaque parameter to the algorithm.
Instead, Beasely and Chu's heuristic feasibility operator is applied to
each new cantidate solution after recombination and mutation. The
operator is basically a greedy solver for creating maximal independent
sets from a random set. It first eliminates vertices using a heuristic
to determine which are most likely to appear in the MIS. Then it adds
back vertices where it can without invalidating the independent set
constraint. In this case, the heuristic is the weight of the vertex
divided by the number of neighbors it has. The algorithm runs in
O(|V|^2) time, where |V| is the number of vertices in the graph.

This operator ensures that all solutions created by the genetic process
are independent sets. It also increases the chance of generating
duplicate members of the population, thereby reducing genetic diversity.
To mitigate this, we ensure each newly created solution is unique,
rerunning the genetic operators if it is not. This poses some
difficulties for the parallel version of the algorithm (see below).

The other major improvement is the use of a "generalised fitness-based
crossover operator" for recombination. It uses a parent's fitness to
determine the probability that each bit appears in the child. Although
theoretically we could use any number of parents, the algorithm
currently only selects two for efficiency reasons. According to Beasely
and Chu, this operator creates more diverse solutions "when the two
parent solutions are similar in structure than the one-point or
two-point crossover operators." The heuristic feasibility operator makes
genetically similar parents very likely, as most solutions will contain
vertices that the heuristic considers good ones.

We also use Beasely and Chu's variable rate mutator, although it is more
of a minor optimization than a major one.

Results
-------

These tests were performed on a 2008 Dell laptop with a 64-bit dual-core
2.10GHz processor. The algorithm should scale almost linearly to
additional cores, so using a more modern computer would greatly improve
performance. I have not yet tested the algorithm with random seeds,
instead using a fixed seed to determine optimal parameters for the
mutator and recombinator. Until testing with random seeds is performed,
these results should be viewed with skepticism. Tuning has focused on
the 450 vertex graphs with a MIS of 30, and the best invocation seemed
to be:

    bin/graph -w 30 -p 600 -s 8 -m 8 {input_file}

This uses a tournament size of 8, a mutation rate of 8, and a population
size of 600. Using these parameters, optimal solutions were found for
BHOSLIB instances 1, 2, and 4 in 7.49s, 7.73s, and 1.63s respectively.
Instance 3 converges to 29 vertices in 3.48s but makes no further
progress through 1e6 iterations. Instance 5 takes even longer to
converge, requiring 175.86s and just over 1e6 iterations to find an
independent set of 29 vertices. This behavior can probably be
attributed to the heuristic feasibility operator, as the MISs for the
first three graphs contain many of the vertices which are positively
evaluated by the heuristic. This shows that while the heuristic
feasibility operator finds good solutions rapidly, it can converge too
quickly to find optimal solutions on some graphs.

For fun, I also ran the genetic algorithm on the frb100 instance which
contains 4000 vertices and has a MIS of 100. The algorithm finds an
independent set of 90 in 140s, but stalls at 92 through 5e5 iterations
or roughly 1600s. According to the BHOSLIB website, the best solutions
are found with stochastic local search solvers, which discover solutions
of over 95 vertices in around 1000s. The current record of 99 vertices
took just under 24 hours to find. I am hopeful that using adaptive
genetic operators could prevent premature convergence of the algorithm
on large graphs such as these.

Implementation Details
----------------------

The solver is written in C++ and makes heavy use of features introduced
in C++11. Currently, it depends on boost's program options, dynamic
bitset, and graph libraries. Most of the implementation
details for the single-threaded case are fairly obvious. We reorder the
input graph to give "good" vertices lower indices which simplifies the
heuristic feasibility operator.

In the parallel version, each thread maintains its own population, and performs
a number of iterations of the single threaded genetic algorithm without
interacting with another thread. At this point, each thread copies some subset
of its population onto a SPSC channel which is read by the next thread spawned.
Each thread then continues iterating, but reads from a second SPSC channel,
expecting the other threads in the cycle to perform roughly the same operation.
There is very little synchronization here, meaning the multi-threaded approach
incurs very little overhead, but it does not guarantee that each thread gets
equal time on the processor. Uniqueness is guaranteed only within a thread's
population, not amongst all threads, as that would require a synchronized
global hash table.

Limitations and Possible Improvements
-------------------------------------

The genetic operators are implemented with virtual base
classes, meaning they incur the overhead of a virtual function call at
every iteration of the algorithm. If performance is critical,  could
templatize the genetic solver and remove the virtual base classes, but I
like the clarity of interfaces.

The main performance bottleneck is the heuristic feasibility operator,
which could be optimized further by rearranging the graph to put highly-rated
vertices first. This would mean we no longer have to check every bit in a
genome, only the set ones.

