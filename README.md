# UAIC_Twin_Width - An Exact and Heuristic solver for the Twin-Width Problem
 Twin-width is a graph-theoretic invariant defined in terms of contractions on trigraphs.  
 This repository provides an exact and heuristic solver for the twin-width problem.
 
 - Exact Solver ([PDF][exact_description])
 - Heuristic Solver ([PDF][heuristic_description])

Requirements
-----------

  - A 64-bit Linux operating system.
  - A 9.4.0 or higher version of `G++` compiler.

Run Application
-----------

Both `heuristic` and `exact` solvers are built within a single `C++` file that reads a twin-width instance from stdin and print the solution to stdout.
For the input and output format, please refer to the [PACE challenge web page](https://pacechallenge.org/2023/io/).

To compile and run, use the following commands:

    g++ -o heuristic -static heuristic.cpp 
    ./heuristic
    
    g++ -o exact -static exact.cpp 
    ./exact

[heuristic_description]: https://andrei-arhire.web.app/assets/UAIC_Twin_Width__A_Heuristic_Twin_Width_Algorithm.pdf "Heuristic Solver Description"
[exact_description]: https://andrei-arhire.web.app/assets/UAIC_Twin_Width__An_Exact_Twin_Width_Algorithm.pdf "Exact Solver Description"
