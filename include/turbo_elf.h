#pragma once
#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include "feedbacktrellis.h"
#include "listdecoder.h"

namespace awgn {

std::vector<double> addNoise(std::vector<int> encodedMsg, double SNR);

} // namespace awgn

namespace crc {

// binary sum, used in crc_check
int bin_sum(int i, int j);

// checks the decoded message against the crc
bool crc_check(std::vector<int> input_data, int crc_bits_num, int crc_dec);

void crc_calculation(std::vector<int>& input_data, int crc_bits_num, int crc_dec);

} // namespace crc

namespace turbo_elf_utils {

// prints a vector of doubles, with commas seperating elements
void print_double_vector(std::vector<double> vector);

// prints a vector of ints, with commas seperating elements
void print_int_vector(std::vector<int> vector);

// outputs a vector of ints to a file
void output_int_vector(std::vector<int> vector, std::ofstream& file);

} // namespace turbo_elf_utils


int make_file_interleaver(char interleaver_file[],
                          unsigned short int interleaver[], int n);

void elf_turbo_simulation(codeInformation code);