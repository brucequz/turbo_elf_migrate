#include "../include/feedforwardtrellis.h"
#include <iostream>
#include <string>

static const int V = 6;

namespace fft_utils {

// converts decimal input to binary output, with a given number of bits
// since we need to keep track of leading zeros
static void dec_to_binary(int input, std::vector<int> &output, int bit_number) {
  output.assign(bit_number, -1);
  for (int i = bit_number - 1; i >= 0; i--) {
    int k = input >> i;
    if (k & 1)
      output[bit_number - 1 - i] = 1;
    else
      output[bit_number - 1 - i] = 0;
  }
}

// converts decimal output to n-bit BPSK
static std::vector<int> get_point(int output, int n) {
  std::vector<int> bin_output;
  dec_to_binary(output, bin_output, n);
  for (int i = 0; i < n; i++) {
    bin_output[i] = -2 * bin_output[i] + 1;
  }
  return bin_output;
}

} // namespace fft_utils

FeedforwardTrellis::FeedforwardTrellis(int k, int n, int v,
                                       std::vector<int> numerators) {
  this->k_ = k;
  this->n_ = n;
  this->v_ = v;
  for (int i = 0; i < numerators.size(); i++) {
    this->numerators_.push_back(numerators[i]);
  }
  this->numInputSymbols_ = pow(2.0, k);
  this->numOutputSymbols_ = pow(2.0, n);
  this->numStates_ = pow(2.0, v);
  this->nextStates_ = std::vector<std::vector<int>>(
      numStates_, std::vector<int>(numInputSymbols_));
  this->outputs_ = std::vector<std::vector<int>>(
      numStates_, std::vector<int>(numInputSymbols_));

  if (v != V) {
    std::cout
        << "MAJOR ISSUE: CONST V DOES NOT MATCH v. EDIT IN FEEDBACKTRELLIS"
        << std::endl;
    exit(1);
  }
  computeNextStates();
}

void FeedforwardTrellis::computeNextStates() {
  // convert to binary numerators
  std::vector<std::vector<int>> bin_numerators(n_, std::vector<int>(v_ + 1));
  for (int i = 0; i < n_; i++) {
    int tempNum = numerators_[i]; // octal number in numerator(135)
    int decIn = 0;
    std::string in = std::to_string(tempNum);
    for (int p = (in.length() - 1); p >= 0; p--)
      decIn += (int)(in[p] - '0') * pow(8, (in.length() - p - 1));
    for (int j = v_; j >= 0; j--) {
      if (decIn % 2 == 0)
        bin_numerators[i][j] = 0;
      else
        bin_numerators[i][j] = 1;
      decIn = decIn / 2;
    }
  }
  // calculate next states and outputs
  for (int currentState = 0; currentState < numStates_; currentState++) {
    std::vector<int> mem_elements = dec2Bin(currentState, v_ + 1);
    for (int input = 0; input < numInputSymbols_; input++) {
      mem_elements[0] = input;
      std::vector<int> output(n_);
      for (int i = 0; i < n_; i++) {
        output[i] = 0;
      }
      for (int x_bit = 0; x_bit < n_; x_bit++) {
        for (int m_bit = 0; m_bit < V + 1; m_bit++) {
          if (bin_numerators[x_bit][m_bit] == 1) {
            output[x_bit] ^= mem_elements[m_bit];
          }
        }
      }
      outputs_[currentState][input] = bin2Dec(output);
      std::vector<int> temp(V);
      for (int i = 0; i < V; i++) {
        temp[i] = mem_elements[i];
        // std::cout << temp[i] << std::endl;
      }
      nextStates_[currentState][input] = bin2Dec(temp);
    }
  }
}

std::vector<int> FeedforwardTrellis::encoder(std::vector<int> originalMessage) {
  // brute force approach, there is a better way to do this assuming
  // invertibility that allows us to precompute starting / ending states,
  // reducing complexity in each encoding from O(numStates) to O(2). revisit
  // when available

  for (int m = 0; m < numStates_; m++) {
    std::vector<int> output;
    int State = m;
    for (int i = 0; i < originalMessage.size(); i += k_) {
      int decimal = 0;
      for (int j = 0; j < k_; j++) {
        decimal += (originalMessage[i + j] * pow(2, k_ - j - 1));
      }
      std::vector<int> outputBinary =
          fft_utils::get_point(outputs_[State][decimal], n_);
      State = nextStates_[State][decimal];
      for (int j = 0; j < n_; j++) {
        output.push_back(outputBinary[j]);
      }
    }
    if (m == State) {
      return output;
    }
  }
  return originalMessage;
}

std::vector<int> FeedforwardTrellis::dec2Bin(int decimal, int length) {
  std::vector<int> binary(length);
  for (int j = (length - 1); j >= 0; j--) {
    if (decimal % 2 == 0)
      binary[j] = 0;
    else
      binary[j] = 1;
    decimal = decimal / 2;
  }
  return binary;
}

int FeedforwardTrellis::bin2Dec(std::vector<int> binary) {
  int decimal = 0;
  for (int i = (binary.size() - 1); i >= 0; i--) {
    decimal += (binary[i] * pow(2, (binary.size() - i - 1)));
  }
  return decimal;
}