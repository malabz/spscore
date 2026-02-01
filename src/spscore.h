#ifndef SPSCORE_H
#define SPSCORE_H

#include <cstdint>
#include <vector>
#include <string>

// Normalize to {'A','C','G','T','N','-'} with U->T and invalid->N.
char normalize_base(char c);

// Calculate binomial coefficient: n choose 2
uint64_t choose2_u64(uint64_t x);

// Calculate SP score for a column given base counts
// cnt: [A,C,G,T,N,'-'] in indices 0..5
double score_of_counts(const uint32_t cnt[6],
                       double matchS, double mismatchS,
                       double gap1S, double gap2S);

// Calculate total SP score for multiple sequence alignment
// Returns: {total_sp, avg_sp, scaled_sp}
struct SPScoreResult {
    double total_sp;
    double avg_sp;
    double scaled_sp;
};

SPScoreResult calculate_sp_score(const std::vector<std::string>& seqs,
                                  double matchS = 1.0,
                                  double mismatchS = -1.0,
                                  double gap1S = -2.0,
                                  double gap2S = 0.0);

// Load FASTA file (supports gzip)
std::vector<std::string> load_fasta(const std::string& path);

#endif // SPSCORE_H

