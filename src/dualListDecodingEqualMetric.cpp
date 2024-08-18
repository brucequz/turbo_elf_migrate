#include "../include/listdecoder.h"

#include <vector>


// DLDInfo DualListDecoding_TurboELF_BAM_Equal_Metric(std::vector<double> txSig_0, std::vector<double> txSig_1, unsigned short int* interleaver_ptr, unsigned short int* deinterleaver_ptr) {
  /*
  Input: 
    - txSig_0: input vector of R1 + Sys symbols
    - txSig_1: input vector of R2 + Sys symbols
    - interleaver_ptr: an array of interleaver indices
    - deinterleaver_ptr: an array of deinterleaver indices

  Output:
    - DLDInfo (struct):
        -  discovered_decoder_idx (int)         : if best available is returned, which decoder discovered it. Default = -1
        -  discovered_partial_metric (double)   : if best available is returned, what is the partial metric. Default = -1.0
        -  discovered_systematic_metric (double): if best available is returned, what is the systematic metric = partial metric - R1/R2 metric. Default = -1.0
        -  undiscovered_partial_metric (double) : if best available is returned, we re-encode, and then compare with the undiscovered codeword, accumulate only the R1/R2 metrics. Default = 
        -  combined_metric (double)             : full metric for all three types of symbols. Default: 0.0
        -  return_type (std::string)            : one of the following: "best_available", "default", "agreed". Default: "default"
        -  message (std::vector<int>)           : decoded message. Default: empty vector
        -  list_ranks (std::vector<int>)        : pair of list ranks. Default: {-1, -1}

  Pseudo Code:
    Initialize minHeap for both Decoders
    Build trellises
    Either a match has been found or 
  */
// }