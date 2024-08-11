#include <iostream>
#include <queue>

#include "../include/feedbacktrellis.h"
#include "../include/galoisfield.h"
#include "../include/galoisfieldpolynomial.h"

static const int V = 6;
// table of primitive polynomials for v=2-15
std::vector<std::string> primitive_polys = {"111",
                                            "1011",
                                            "11001",
                                            "111101",
                                            "1100001",
                                            "10110111",
                                            "110110001",
                                            "1011010111",
                                            "11000001011",
                                            "110111100111",
                                            "1110111110101",
                                            "11111011110011",
                                            "100000101000011",
                                            "1000010000100011"};

FeedbackTrellis::FeedbackTrellis(int k, int n, int v,
                                 std::vector<int> numerators, int denominator) {
  this->k_ = k;
  this->n_ = n;
  this->v_ = v;
  this->denominator_ = denominator;
  for (int i = 0; i < numerators.size(); i++) {
    this->numerators_.push_back(numerators[i]);
  }
  this->numInputSymbols_ = pow(2.0, k_);
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
  // computeTerminations();
}

void FeedbackTrellis::computeNextStates() {
  std::vector<std::vector<int>> gs(k_, std::vector<int>(v_ + 1)); // rename plz
  unsigned int b_arr[V + 1];

  for (int i = 0; i < k_; i++) {
    int tempNum = numerators_[i]; // octal number in numerator(135)
    int decIn = 0;
    std::string in = std::to_string(tempNum);
    for (int p = (in.length() - 1); p >= 0; p--)
      decIn += (int)(in[p] - '0') * pow(8, (in.length() - p - 1));
    for (int j = v_; j >= 0; j--) {
      if (decIn % 2 == 0)
        gs[k_ - i - 1][j] = 0;
      else
        gs[k_ - i - 1][j] = 1;
      decIn = decIn / 2;
    }
  }

  // first, construct GF outside the loop with a primitive poly
  unsigned int base_poly[V + 1];
  for (int i = 0; i < V + 1; i++) {
    base_poly[i] = primitive_polys[V - 2][i] - '0';
  }
  galois::GaloisField gf1 = galois::GaloisField(V, base_poly);

  // convert denominator to binary
  std::string b;
  for (int j = 0; j <= v_; j++) {
    b += '1';
  }
  int decIn = 0;
  std::string in = std::to_string(denominator_);
  for (int p = (in.length() - 1); p >= 0; p--)
    decIn += (int)(in[p] - '0') * pow(8, (in.length() - p - 1));
  for (int j = v_; j >= 0; j--) {
    if (decIn % 2 == 0)
      b[j] = '0';
    decIn = decIn / 2;
  }
  for (int i = 0; i < V + 1; i++) {
    b_arr[i] = b[i] - '0';
  }

  // computing nextStates
  for (int currentState = 0; currentState < numStates_; currentState++) {
    std::vector<int> cur_sigmas = dec2Bin(currentState, v_);
    for (int input = 0; input < numInputSymbols_; input++) {
      std::vector<int> us = dec2Bin(input, k_);
      unsigned int total_numerator[V + 1];
      for (int i = 0; i < (V + 1); i++) {
        total_numerator[i] = 0;
      }
      for (int ii = 0; ii < us.size(); ii++) {
        std::vector<int> temp(v_ + 1);
        for (int i = 0; i < (v_ + 1); i++) {
          temp[i] = gs[ii][i] * us[ii]; // GFCONV
        }
        for (int i = 0; i < (v_ + 1); i++) {
          total_numerator[i] = temp[i] ^ total_numerator[i];
        }
      }
      for (int i = 1; i < (v_ + 1); i++) {
        total_numerator[i] = cur_sigmas[i - 1] ^ total_numerator[i];
      }

      // //	DEBUGGING: USED TO DISPLAY NUMERATOR AND DENOMINATOR FOR GFDIV
      // std::cout << "_____________________________________________________" <<
      // std::endl; std::cout << "total numerator: " << std::endl; for(int i =
      // 0; i < V + 1; i++){ 	std::cout << total_numerator[i] << "\t";
      // }
      // std::cout << std::endl;
      // std::cout << "b_arr: " << std::endl;
      // for(int i = 0; i < V + 1; i++){
      // 	std::cout << b_arr[i] << "\t";
      // }
      // std::cout << "decIn: " << decIn << std::endl;

      galois::GaloisFieldElement gfe1[V + 1];
      galois::GaloisFieldElement gfe2[V + 1];

      for (int i = 0; i < (V + 1); i++) {
        gfe1[i] = galois::GaloisFieldElement(&gf1, total_numerator[i]);
      }
      for (int i = 0; i < (V + 1); i++) {
        gfe2[i] = galois::GaloisFieldElement(&gf1, b_arr[i]);
      }
      galois::GaloisFieldPolynomial poly1 =
          galois::GaloisFieldPolynomial(&gf1, V, gfe1);
      galois::GaloisFieldPolynomial poly2 =
          galois::GaloisFieldPolynomial(&gf1, V, gfe2);

      galois::GaloisFieldPolynomial polyremd = poly1 % poly2;
      galois::GaloisFieldPolynomial polyq = poly1 / poly2;

      std::vector<galois::GaloisFieldElement> remdTemp = polyremd.getPoly();
      std::vector<galois::GaloisFieldElement> qTemp = polyq.getPoly();
      std::vector<int> remd;
      std::vector<int> q;
      for (int i = 0; i < remdTemp.size(); i++) {
        remd.push_back(remdTemp[i].poly());
      }
      for (int i = 0; i < qTemp.size(); i++) {
        q.push_back(qTemp[i].poly());
      }
      if (remd.size() == 0) {
        remd.push_back(0);
      }
      if (q.size() == 0) {
        q.push_back(0);
      }

      std::vector<int> testRemd;
      std::vector<int> testQ;

      if (total_numerator[v_] == 1) {
        testQ.push_back(1);
        for (int i = 0; i < v_ + 1; i++) {
          // for(int i = 0; i < v; i++){
          testRemd.push_back(total_numerator[i] ^ b_arr[i]);
        }
      } else {
        testQ.push_back(0);
        for (int i = 0; i < v_ + 1; i++) {
          testRemd.push_back(total_numerator[i]);
        }
      }
      while (testRemd.size() > 1 && testRemd[testRemd.size() - 1] == 0) {
        testRemd.pop_back();
      }

      // VERBOSE DEBUGGING: USED TO CHECK EQUALITY BETWEEN MY GFDIV AND PROPER
      // GFDIV std::cout << std::endl << "remd: " << std::endl; for(int i = 0; i
      // < remd.size(); i++){ 	std::cout << remd[i] << "\t";
      // }
      // std::cout << std::endl << "testRemd: " << std::endl;
      // for(int i = 0; i < testRemd.size(); i++){
      // 	std::cout << testRemd[i] << "\t";
      // }
      // std::cout << std::endl << "quotient: " << std::endl;
      // for(int i = 0; i < q.size(); i++){
      // 	std::cout << q[i] << std::endl;
      // }
      // std::cout << std::endl << "test quotient: " << std::endl;
      // for(int i = 0; i < testQ.size(); i++){
      // 	std::cout << testQ[i] << std::endl;
      // }

      // CHECKS EQUALITY BETWEEN GALOIS FIELD LIBRARY AND MY IMPLEMENTATION-
      // SHOULD BE INCLUDED UNTIL EQS ARE VALIDATED ACROSS ALL RELEVENT CODES
      if (V != v_) {
        std::cout << "constant V does not match- check feedbacktrellis"
                  << std::endl;
      }
      if (testRemd != remd) {
        std::cout << "gfdiv machine broke" << std::endl;
        return;
      }
      if (testQ != q) {
        std::cout << "gfdiv machine broke 2" << std::endl;
        return;
      }

      std::vector<int> next_sigmas(v_, 0);
      for (int i = 0; i < v_ && i < remd.size(); i++) {
        next_sigmas[i] = remd[i];
      }

      int nextState = 0;
      for (int i = (next_sigmas.size() - 1); i >= 0; i--) {
        nextState += next_sigmas[i] *
                     pow(2, (next_sigmas.size() - 1 -
                             i)); // convert next_sigmas back into decimal value
      }

      nextStates_[currentState][input] = nextState;
      std::vector<int> nextOutput(n_); // k_+1
      nextOutput[0] = q[0];            // q has to be size 0 or 1
      for (int i = 1; i < n_; i++) {
        nextOutput[i] = us[i - 1];
      }
      int output = 0;
      for (int i = (nextOutput.size() - 1); i >= 0; i--) {
        output += nextOutput[i] * pow(2, nextOutput.size() - 1 - i);
      }
      outputs_[currentState][input] = output;
      // std::cout << output << std::endl;
    }
  }
}
std::vector<int> FeedbackTrellis::encoder(std::vector<int> originalMessage) {
  // brute force approach, there is a better way to do this assuming
  // invertibility that allows us to precompute starting / ending states,
  // reducing complexity in each encoding from O(numStates_) to O(2). revisit
  // when available

  for (int m = 0; m < numStates_; m++) {
    std::vector<int> output;
    int State = m;
    for (int i = 0; i < originalMessage.size(); i += k_) {
      int decimal = 0;
      for (int j = 0; j < k_; j++) {
        decimal += (originalMessage[i + j] * pow(2, k_ - j - 1));
      }
      std::vector<int> outputBinary = dec2Bin(outputs_[State][decimal], n_);
      State = nextStates_[State][decimal];
      for (int j = 0; j < n_; j++) {
        output.push_back(-1 * ((outputBinary[j] * 2) - 1));
      }
    }
    if (m == State) {
      return output;
    }
  }
  return originalMessage;
}

std::vector<int> FeedbackTrellis::dec2Bin(int decimal, int length) {
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

// termination bits for primal trellis
void FeedbackTrellis::computeTerminations() {

  int num_transitions = std::ceil((double)v_ / (n_ - 1));
  std::vector<int> tree(numStates_);
  std::vector<int> vis(numStates_);
  vis[0] = 1;
  std::queue<int> queue;
  queue.push(0);
  while (!queue.empty()) {
    int target_state = queue.front();
    queue.pop();
    for (int i = 0; i < numStates_; i++) {
      if (std::count(nextStates_[i].begin(), nextStates_[i].end(),
                     target_state) > 0 &&
          vis[i] == 0) {
        vis[i] = 1;
        queue.push(i);
        tree[i] = target_state;
      }
    }
  }

  terminations_ = std::vector<std::vector<int>>(
      numStates_, std::vector<int>(num_transitions));
  for (int i = 0; i < num_transitions; i++) {
    terminations_[0][i] = 0;
  }
  int cur = 0;
  for (int i = 1; i < numStates_; i++) {
    cur = i;
    for (int j = 0; j < num_transitions; j++) {
      int fa_state = tree[cur];
      int index = std::distance(nextStates_[cur].begin(),
                                std::find(nextStates_[cur].begin(),
                                          nextStates_[cur].end(), fa_state));
      terminations_[i][j] = index;
      cur = fa_state;
    }
  }
  std::cout << "made it to the end of the termination computations"
            << std::endl;
}

std::vector<int>
FeedbackTrellis::terminateMsg(std::vector<int> originalMessage) {
  // first encode to find the final state
  std::vector<int> output;
  int finalState = 0;
  for (int i = 0; i < originalMessage.size(); i += k_) {
    int decimal = 0;
    for (int j = 0; j < k_; j++) {
      decimal += (originalMessage[i + j] * pow(2, k_ - j - 1));
    }
    std::vector<int> outputBinary = dec2Bin(outputs_[finalState][decimal], n_);
    finalState = nextStates_[finalState][decimal];
  }

  for (int i = 0; i < originalMessage.size(); i++) {
    output.push_back(originalMessage[i]);
  }

  // appending the termination bits that force the zero state
  std::vector<int> localTerminations = terminations_[finalState];
  for (int i = 0; i < localTerminations.size(); i++) {
    std::vector<int> outputBinary = dec2Bin(localTerminations[i], k_);
    for (int j = 0; j < k_; j++) {
      output.push_back(outputBinary[j]);
    }
  }
  return output;
}

std::vector<int>
FeedbackTrellis::ztencoder(std::vector<int> terminatedMessage) {
  // unlike for tbcc, we always start at the zero state, add the crc, then add
  // bits to the message to guarentee we return to the zero state when we encode

  std::vector<int> output;
  int State = 0;
  for (int i = 0; i < terminatedMessage.size(); i += k_) {
    int decimal = 0;
    for (int j = 0; j < k_; j++) {
      decimal += (terminatedMessage[i + j] * pow(2, k_ - j - 1));
    }
    std::vector<int> outputBinary = dec2Bin(outputs_[State][decimal], n_);
    State = nextStates_[State][decimal];
    for (int j = 0; j < n_; j++) {
      output.push_back(-1 * ((outputBinary[j] * 2) - 1));
    }
  }
  if (State != 0) {
    std::cout << "error: not zero terminated" << std::endl;
  }
  return output;
}