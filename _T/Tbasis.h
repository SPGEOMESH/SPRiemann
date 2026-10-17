
#pragma once

#include <iostream>
using std::cout; using std::endl; using std::cerr;
#include <iomanip>
#include <fstream>
#include <filesystem>

#include <kiwi/mpi/sim.h>

// testing: gate for a serial processing step (e.g., post-processing) in a test.
bool test_serial_processing() {
   MPI_Barrier(MPI_COMM_WORLD);
   int node_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &node_rank);
   return (node_rank == 0);
}

void test_end_serial_processing() {
   MPI_Barrier(MPI_COMM_WORLD);
}


void savetxt(const std::string& filename, const std::string& content) {
   std::ofstream s;
   s.open(filename.c_str(), std::ofstream::out);
   s << content;
   s.close();
}

// delete a file
bool delete_file(const std::string& filename) {
   return std::filesystem::remove(filename);
}

// compare doubles
bool isclose(double x, double y) {
   static const double epsilon = 1e-10;
   return (x == y) || (x < y && x + epsilon > y) || (y < x && y + epsilon > x);
}

// todo work in progress.
std::vector<double> loadtxt(const std::string& filename) {
   std::fstream s;
   s.open(filename.c_str(), std::fstream::in);
   std::vector<double> out;
   double x;
   while (s >> x) out.push_back(x);
   s.close();
   return out;
}


//#define value_precision 6
//#define column_size 15
//#define value(x) std::setw(column_size) << std::setprecision(value_precision) << (x)

