
#pragma once

#include <iostream>
using std::cout; using std::endl; using std::cerr;
#include <iomanip>
#include <fstream>

#include <kiwi/mpi/sim.h>

void savetxt(const std::string& filename, const std::string& content) {
   std::ofstream s;
   s.open(filename.c_str(), std::ofstream::out);
   s << content;
   s.close();
}


//#define value_precision 6
//#define column_size 15
//#define value(x) std::setw(column_size) << std::setprecision(value_precision) << (x)

