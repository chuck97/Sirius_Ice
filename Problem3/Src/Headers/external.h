#pragma once

#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <iomanip>
#include <cmath>
#include <map>
#include <algorithm>
#include <cctype>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <iomanip>
#include <functional>
#include <limits>
#include <initializer_list>
#include "inmost.h"
#include <nlohmann/json.hpp>

#define USE_MPI
#define USE_PARTITIONER
#define USE_PARTITIONER_PARMETIS
#define USE_SOLVER


#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

#define INMOST_ICE_ERR(message) {std::cerr << "Error: " << message  << std::endl; MPI_Finalize(); exit(1);}
#define VAREPS 1e-15
