#include <algorithm>

#include "../include/turbo_elf.h"
#include "../include/constants.h"
#include "../include/feedbacktrellis.h"
#include "../include/listdecoder.h"

namespace awgn {

std::default_random_engine generator;

std::vector<double> addNoise(std::vector<int> encodedMsg, double SNR) {
  std::vector<double> noisyMsg;

  double variance = pow(10.0, -SNR / 10.0);
  double sigma = sqrt(variance);
  std::normal_distribution<double> distribution(0.0, sigma);

  for (int i = 0; i < encodedMsg.size(); i++) {
    noisyMsg.push_back(encodedMsg[i] + distribution(generator));
  }
  return noisyMsg;
}

} // namespace awgn

namespace crc {

// binary sum, used in crc_check
int bin_sum(int i, int j) {
	return (i + j) % 2;
}

// checks the decoded message against the crc
bool crc_check(std::vector<int> input_data, int crc_bits_num, int crc_dec) {
	std::vector<int> CRC;
	fbt_utils::dec_to_binary(crc_dec, CRC, crc_bits_num);

	for (int ii = 0; ii <= (int)input_data.size() - crc_bits_num; ii++) {
		if (input_data[ii] == 1) {
			// Note: transform doesn't include .end
			std::transform(input_data.begin() + ii, input_data.begin() + (ii + crc_bits_num), CRC.begin(), input_data.begin() + ii, bin_sum);
		}
	}
	bool zeros = std::all_of(input_data.begin(), input_data.end(), [](int i) { return i == 0; });
	return zeros;
}

void crc_calculation(std::vector<int>& input_data, int crc_bits_num, int crc_dec){
	// crc_bits_num: the number of CRC bits, redundancy bits number is 1 less.
	int length = (int)input_data.size();
	std::vector<int> CRC;
	fbt_utils::dec_to_binary(crc_dec, CRC, crc_bits_num);
	input_data.resize(length + crc_bits_num - 1, 0);

	std::vector<int> output_data = input_data;
	for (int ii = 0; ii <= length - 1; ii++)
	{
		if (output_data[ii] == 1)
		{
			std::transform(output_data.begin() + ii, output_data.begin() + (ii + crc_bits_num), CRC.begin(), output_data.begin() + ii, bin_sum);
		}
	}

	for (int ii = length; ii < (int)output_data.size(); ii++){ input_data[ii] = output_data[ii];}
}

} // namespace crc

namespace turbo_elf_utils {

// prints a vector of doubles, with commas seperating elements
void print_double_vector(std::vector<double> vector){
	if(vector.size() == 0)
		return;
	for(int i = 0; i < vector.size() - 1; i++){
		std::cout << vector[i] << ", ";
	}
	std::cout << vector[vector.size() - 1] << std::endl;
}

// prints a vector of ints, with commas seperating elements
void print_int_vector(std::vector<int> vector){
	if(vector.size() == 0)
		return;
	for(int i = 0; i < vector.size() - 1; i++){
		std::cout << vector[i] << ", ";
	}
	std::cout << vector[vector.size() - 1] << std::endl;
}

// outputs a vector of ints to a file
void output_int_vector(std::vector<int> vector, std::ofstream& file){
	if(vector.size() == 0)
		return;
	for(int i = 0; i < vector.size() - 1; i++){
		file << vector[i] << ", ";
	}
	file << vector[vector.size() - 1] << std::endl;
}

} // namespace turbo_elf_utils

int make_file_interleaver(char interleaver_file[],
                          unsigned short int interleaver[], int n) {
  FILE *fp;
  if ((fp = fopen(interleaver_file, "r")) == NULL) {
    printf("Error opening interleaver file for reading.\n");
    printf("%s\n", interleaver_file);
    exit(1); /* exit program with error condition  */
  }
  printf("Reading interleaver from file.\n");
  for (int i = 0; i < n; i++) {
    if (fscanf(fp, "%hu", &interleaver[i])) {
    };
    /* Translate from 1..K to 0..(K-1) */
    //       interleaver[i]--;
  }

  printf("%s\n", interleaver_file);
  fclose(fp);
  return 0;
}

void elf_turbo_simulation(codeInformation code, std::vector<double> SNR) {

  std::ofstream outFile;
  outFile.open("../output/output.txt");

  int Km = code.numInfoBits + code.crcDeg - 1; // 71

  // puncture 43 bits to achieve overall rate-1/2 when m=7
  std::vector<int> punc_idx = {0,  1,  2,  3,  4,  6,  9,  13, 14, 15, 18,
                               19, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31,
                               32, 33, 34, 35, 36, 37, 38, 39, 40, 42, 45,
                               46, 52, 53, 54, 55, 59, 60, 61, 62, 63};

  FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators,
                                  code.denominator);

  // interleaver
  char interleaver_file[1024]; // make sure the size is large enough
  // only K+m because half of the bits are systematic
  snprintf(interleaver_file, sizeof(interleaver_file),
           "../data/interleaver_%d_S5_T1.dat", Km);
  unsigned short int interleaver[Km]; // Adjust size as needed
  make_file_interleaver(interleaver_file, interleaver, Km);
  unsigned short int deinterleaver[Km]; // Adjust size as needed
  for (int i = 0; i < Km; i++) {
    deinterleaver[interleaver[i]] = i;
  }

  // outer loop: SNR
  for (int s = 0; s < SNR.size(); s++) {
    double cur_SNR = SNR[s];

    int numerror = 0;
    int numtrial = 0;
    int default_cnt = 0;
    int agreed_cnt = 0;
    int best_available_cnt = 0;
    std::vector<double> metric_error_diff = {};

    std::tuple<int, int, int> error_counts{0, 0, 0};
    // inner loop: MC trials
    while (numerror < MAXERRORS && numtrial < NUMTRIALS) {

      if (numtrial % 1000 == 0) {
        std::cout << "trial number: " << numtrial << std::endl;
      }

      std::vector<int> original_message;
      for (int i = 0; i < code.numInfoBits; i++) {
        // original_message.push_back(rand() % 2);
        original_message.push_back(0);
      }

      crc::crc_calculation(original_message, code.crcDeg, code.crc);

      std::vector<int> encodedMessage =
          encodingTrellis.encoder(original_message);

      // get parity bits
      std::vector<int> X_R1;
      std::vector<int> X_sys;
      for (int j = 0; j < encodedMessage.size(); j += 2) {
        X_sys.push_back(encodedMessage[j + 1]);
        X_R1.push_back(encodedMessage[j]);
      }

      // interleaved original message
      std::vector<int> pi_original_message;
      for (int ii = 0; ii < Km; ii++) {
        pi_original_message.push_back(original_message[interleaver[ii]]);
      }

      // encode interleaved message
      std::vector<int> encodedMessage_inter =
          encodingTrellis.encoder(pi_original_message);

      // get parity bits
      std::vector<int> X_R2;
      for (int j = 0; j < encodedMessage_inter.size(); j += 2) {
        X_R2.push_back(encodedMessage_inter[j]);
      }

      // add noise

      // std::vector<double> Y_R1;
      // std::vector<double> Y_R2;
      // std::vector<double> Y_sys;

      // for (int i = 0; i < X_R1.size(); i++) {
      //   Y_R1.push_back(static_cast<double>(X_R1[i]));
      // }

      // for (int i = 0; i < X_R2.size(); i++) {
      //   Y_R2.push_back(static_cast<double>(X_R2[i]));
      // }

      // for (int i = 0; i < X_sys.size(); i++) {
      //   Y_sys.push_back(static_cast<double>(X_sys[i]));
      // }
      std::vector<double> Y_R1 = awgn::addNoise(X_R1, cur_SNR);
      std::vector<double> Y_R2 = awgn::addNoise(X_R2, cur_SNR);
      std::vector<double> Y_sys = awgn::addNoise(X_sys, cur_SNR);

      // puncture both parity sequences
      for (int p = 0; p < Y_R1.size(); p++) {
        for (int q = 0; q < punc_idx.size(); q++) {
          if (p == punc_idx[q]) {
            Y_R1[p] = 0;
            Y_R2[p] = 0;
          }
        }
      }
      // interleave Y_sys for bottom list decoder
      std::vector<double> pi_Y_sys;
      for (int ii = 0; ii < Km; ii++) {
        pi_Y_sys.push_back(Y_sys[interleaver[ii]]);
      }

      std::vector<double> squaredDifference_R1(X_R1.size());
      std::vector<double> squaredDifference_R2(X_R2.size());
      std::vector<double> squaredDifference_sys(X_sys.size());

      for (int i = 0; i < X_R1.size(); i++) {
        squaredDifference_R1[i] =
            (std::find(punc_idx.begin(), punc_idx.end(), i) == punc_idx.end())
                ? pow(X_R1[i] - Y_R1[i], 2)
                : 0;
      }
      double R1_squraed_diff = std::accumulate(squaredDifference_R1.begin(),
                                               squaredDifference_R1.end(), 0.0);

      for (int i = 0; i < X_R2.size(); i++) {
        squaredDifference_R2[i] =
            (std::find(punc_idx.begin(), punc_idx.end(), i) == punc_idx.end())
                ? pow(X_R2[i] - Y_R2[i], 2)
                : 0;
      }
      double R2_squraed_diff = std::accumulate(squaredDifference_R2.begin(),
                                               squaredDifference_R2.end(), 0.0);

      for (int i = 0; i < X_sys.size(); i++) {
        squaredDifference_sys[i] = pow(X_sys[i] - Y_sys[i], 2);
      }
      double sys_squraed_diff = std::accumulate(
          squaredDifference_sys.begin(), squaredDifference_sys.end(), 0.0);

      double sum_squared_diff =
          R1_squraed_diff + R2_squraed_diff + sys_squraed_diff;

      std::vector<double> DLD_R1;
      std::vector<double> DLD_R2;

      for (int i = 0; i < Y_R1.size(); i++) {
        DLD_R1.push_back(Y_R1[i]);
        DLD_R1.push_back(Y_sys[i]);
      }
      for (int i = 0; i < Y_R2.size(); i++) {
        DLD_R2.push_back(Y_R2[i]);
        DLD_R2.push_back(pi_Y_sys[i]);
      }

      std::vector<codeInformation> codeList = {code, code};
      DualListDecoder DLD(codeList, LISTSIZE, punc_idx);
      DLDInfo result = DLD.DualListDecoding_TurboELF_BAM_distance_spectrum(
          DLD_R1, DLD_R2, interleaver, deinterleaver);

      // process return type
      if (result.return_type == "default") {
        default_cnt++;
      } else if (result.return_type == "agreed") {
        agreed_cnt++;
      } else if (result.return_type == "best_available") {
        best_available_cnt++;
      }

      numtrial++;
      assert(default_cnt + agreed_cnt + best_available_cnt == numtrial);

      if (result.message.size() == original_message.size() &&
          result.message == original_message) {
        assert(std::fabs(result.combined_metric - sum_squared_diff) < 0.0001);
      } else if (result.message.size() != original_message.size() ||
                 result.message != original_message) {
        numerror++;
        assert(result.combined_metric - sum_squared_diff > 0.0001);
        metric_error_diff.push_back(result.combined_metric - sum_squared_diff);
        if (result.return_type == "default") {
          std::get<0>(error_counts)++;
        } else if (result.return_type == "agreed") {
          std::get<1>(error_counts)++;
        } else if (result.return_type == "best_available") {
          std::get<2>(error_counts)++;
        }

        if (outFile.is_open()) {
          outFile << "//////// ERROR EVENT ///////// " << std::endl;
          outFile << "SNR: " << cur_SNR << " , max list size: " << LISTSIZE
                  << std::endl;
          outFile << "trial number: " << numtrial << std::endl;
          outFile << "error type: " << result.return_type << std::endl;
          outFile << "return metric: " << result.combined_metric
                  << ", found at: ";
          if (result.return_type == "best_available") {
            outFile << "discovered decoder: " << result.discovered_decoder_idx
                    << ", partial metric: " << result.discovered_partial_metric
                    << std::endl;
            outFile << "undiscovered metric: "
                    << result.undiscovered_partial_metric << std::endl;
          }
          outFile << "true metric: " << R1_squraed_diff << " + "
                  << sys_squraed_diff << " + " << R2_squraed_diff << " -> "
                  << sum_squared_diff << std::endl;
          outFile << std::endl;
        } else {
          std::cerr << "Unable to open file" << std::endl;
        }

      } else {
        std::cout << "error in error counting" << std::endl;
      }
    } // while (numerror < MAXERRORS && numtrial < NUMTRIALS)

    std::cout << "For SNR = " << cur_SNR << ", number of trials: " << numtrial
              << ", number of errors: " << numerror << std::endl;
    std::cout << "FER: " << (double)numerror / numtrial << std::endl;
    std::cout << "default: " << default_cnt << ", agreed: " << agreed_cnt
              << ", best_available: " << best_available_cnt << std::endl;
    std::cout << "default error: " << std::get<0>(error_counts)
              << ", agreed error: " << std::get<1>(error_counts)
              << ", best_available error: " << std::get<2>(error_counts)
              << std::endl;

    std::cout << "printing metric error difference: ";
    turbo_elf_utils::print_double_vector(metric_error_diff);
    std::cout << std::endl;
  } // for (int s = 0; s < SNR.size(); s++)

  outFile.close();
}