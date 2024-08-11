#ifndef FEEDFORWARDTRELLIS_H
#define FEEDFORWARDTRELLIS_H
#include <vector>
#include <string>

class FeedforwardTrellis {
public:
	FeedforwardTrellis(int k, int n, int v, std::vector<int> numerators);
	std::vector<int> encoder(std::vector<int> originalMessage);
	std::vector<std::vector<int>> getNextStates() {return nextStates_;};
	std::vector<std::vector<int>> getOutputs() {return outputs_;};
	int getNumInputSymbols() {return numInputSymbols_;};
	int getNumOutputSymbols() {return numOutputSymbols_;};
	int getNumStates() {return numStates_;};
	int getV() {return v_;};
	int getN() {return n_;};
private:
	int k_;
	int n_;
	int v_;
	int numInputSymbols_;
	int numOutputSymbols_;
	int numStates_;
	std::vector<int> numerators_;
	std::vector<std::vector<int>> nextStates_;
	std::vector<std::vector<int>> outputs_;

	void computeNextStates();

	std::vector<int> dec2Bin(int decimal, int length);
    int bin2Dec(std::vector<int> binary);
};
#endif