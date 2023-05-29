# UAIC_Twin_Width - An Exact and Heuristic solver for the Twin-Width Problem
 Twin-width is a graph-theoretic invariant defined in terms of contractions on trigraphs.  
 This repository provides an exact and heuristic solver for the twin-width problem.
 
Requirements
-----------

  - A 64-bit Linux operating system.
  - A 9.4.0 or higher version of `G++` compiler.

Run Application
-----------

Both `heuristic` and `exact` solvers are built within a single `C++` file that reads a twin-width instance from stdin and print the solution to stdout.
For the input and output format, please refer to the [PACE challenge web page](https://pacechallenge.org/2023/io/).

To compile and run, use the following commands:

    g++ heuristic.cpp -o heuristic
    ./heuristic
    
    g++ exact.cpp -o exact
    ./exact

