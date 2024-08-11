
#include <iostream>
#include <vector>
#include <random>

#include "../include/constants.h"
#include "../include/galoisfield.h"
#include "../include/galoisfieldpolynomial.h"
#include "../include/galoisfieldelement.h"
#include "../include/dualtrellis.h"
#include "../include/feedbacktrellis.h"
#include "../include/feedforwardtrellis.h"
#include "../include/listdecoder.h"
#include "../include/minheap.h"

int MAXERRORS = 100;
int NUMTRIALS = 1e4;
int LISTSIZE  = 1000;

namespace {

std::vector<double>
ComputeSquaredDifferences(const std::vector<int> &vector1,
                          const std::vector<double> &vector2) {
  // Check if the vectors have the same size
  if (vector1.size() != vector2.size()) {
    // You can handle this error in your preferred way, e.g., throw an exception
    throw std::invalid_argument("Vectors must have the same size");
  }

  // Calculate squared differences
  std::vector<double> squaredDifferences;
  squaredDifferences.reserve(vector1.size()); // Reserve space for efficiency

  for (std::size_t i = 0; i < vector1.size(); ++i) {
    double diff = vector1[i] - vector2[i];
    squaredDifferences.push_back(diff * diff);
  }

  return squaredDifferences;
}

} // namespace

static std::default_random_engine generator;

void elf_turbo_simulation(codeInformation code);

int main() {
  codeInformation code;
  code.k = 1;      // numerator of the rate
  code.n = 2;      // denominator of the rate
  code.v = 6;      // number of memory elements
  code.crcDeg = 8; // m+1, degree of CRC, # bits of CRC polynomial
  code.crc = 215;
  code.numInfoBits = 64;

  // optimal code numerators and denominator are known, and are given in octal
  code.numerators = {133};
  code.denominator = 171;
  code.hMatrix = {{1, 0, 0, 1, 1, 1, 1}, {1, 1, 0, 1, 1, 0, 1}};

  elf_turbo_simulation(code);

  return 0;
}
