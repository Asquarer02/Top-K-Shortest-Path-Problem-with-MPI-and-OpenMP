# Top-K-Shortest-Path-Problem-with-MPI-and-OpenMP

## Overview
This project implements a parallel solution for the Top K Shortest Path Problem using a combination of MPI (Message Passing Interface) for distributed computing and OpenMP for shared memory parallelization within MPI processes. The goal is to efficiently find the Kth shortest path among all nodes of a given graph.

## Features
- Parallel implementation of the Kth Shortest Path algorithm
- Utilizes MPI for distributed computing across multiple processes
- Employs OpenMP for shared memory parallelization within each MPI process
- Handles large-scale graphs efficiently
- Supports weighted graphs
- Includes preprocessing capabilities for adapting various graph representations

## Technologies Used
- C++
- MPI (Message Passing Interface)
- OpenMP

## Installation and Setup
1. Ensure you have an MPI implementation (e.g., OpenMPI, MPICH) installed on your system
2. Install OpenMP (typically included with modern C++ compilers)
3. Clone the repository:
   ```
   git clone https://github.com/Asquarer02/Top-K-Shortest-Path-Problem-with-MPI-and-OpenMP.git
   cd top-k-shortest-path
   ```
4. Compile the project:
   ```
   mpic++ -fopenmp -o ksp main.cpp
   ```

## Usage
1. Prepare your input graph file in the required format
2. Run the program:
   ```
   mpirun -np <num_processes> ./ksp <input_file> <K> <num_pairs>
   ```
   Where:
   - `<num_processes>` is the number of MPI processes to use
   - `<input_file>` is the path to your graph file
   - `<K>` is the Kth shortest path you want to find
   - `<num_pairs>` is the number of random node pairs to test

## Datasets Tested
- Doctor Who Network
- Enron Email Network
- EU Email Network

## Contributing
Contributions to improve the Top K Shortest Path Solver are welcome. Please feel free to submit a pull request or open an issue for discussion.

## Contact
Ahmed Abd-ur-Rahman - i210404@nu.edu.pk

Project Link: [https://github.com/Asquarer02/Top-K-Shortest-Path-Problem-with-MPI-and-OpenMP](https://github.com/Asquarer02/Top-K-Shortest-Path-Problem-with-MPI-and-OpenMP)
