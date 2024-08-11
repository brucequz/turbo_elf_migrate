#ifndef FEEDBACKTRELLIS_H
#define FEEDBACKTRELLIS_H
#include <vector>

namespace fbt_utils {

// converts decimal input to binary output, with a given number of bits
// since we need to keep track of leading zeros
static void dec_to_binary(int input, std::vector<int>& output, int bit_number) {
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
	for (int i=0; i<n; i++){
		bin_output[i] = -2 * bin_output[i] + 1;
	}
	return bin_output;
}

} 

class FeedbackTrellis {
public:
	FeedbackTrellis(int k, int n, int v, std::vector<int> numerators, int denominator);
	std::vector<int> ztencoder(std::vector<int> originalMessage);
	std::vector<int> terminateMsg(std::vector<int> originalMessage);
	std::vector<int> encoder(std::vector<int> originalMessage);
	std::vector<std::vector<int> > getNextStates() {return nextStates_;};
	std::vector<std::vector<int> > getOutputs() {return outputs_;};
	std::vector<std::vector<int> > getTerminations() {return terminations_;};
	int getNumInputSymbols() {return numInputSymbols_;};
	int getNumOutputSymbols() {return numOutputSymbols_;};
	int getNumStates() {return numStates_;};
private:
	int k_;
	int n_;
	int v_;
	int numInputSymbols_;
	int numOutputSymbols_;
	int numStates_;
	std::vector<int> numerators_;
	int denominator_;
	std::vector<std::vector<int> > nextStates_;
	std::vector<std::vector<int> > outputs_;
	std::vector<std::vector<int> > terminations_;

	void computeTerminations();
	void computeNextStates();

	std::vector<int> dec2Bin(int decimal, int length);
};
#endif
