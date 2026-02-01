// sp_omp_kseq.cpp
// Requires: kseq.h (in include path) + zlib
//
// Build:
//   g++ -O3 -march=native -fopenmp sp_omp_kseq.cpp -lz -o sp
//
// Run:
//   ./sp -i alignment.fa[.gz] [--match v] [--mismatch v] [--gap1 v] [--gap2 v]

#include "spscore.h"
#include <zlib.h>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#ifdef _OPENMP
  #include <omp.h>
#endif

// Your local kseq.h
#include "kseq.h"

// kseq needs this macro instantiation for gzFile:
KSEQ_INIT(gzFile, gzread)

// Normalize to {'A','C','G','T','N','-'} with U->T and invalid->N.
char normalize_base(char c) {
    switch (c) {
        case 'A': case 'a': return 'A';
        case 'C': case 'c': return 'C';
        case 'G': case 'g': return 'G';
        case 'T': case 't': return 'T';
        case 'U': case 'u': return 'T';
        case 'N': case 'n': return 'N';
        case '-': return '-';
        default: return 'N';
    }
}

uint64_t choose2_u64(uint64_t x) {
    return (x < 2) ? 0ULL : (x * (x - 1)) / 2ULL;
}

// cnt: [A,C,G,T,N,'-'] in indices 0..5
double score_of_counts(const uint32_t cnt[6],
                       double matchS, double mismatchS,
                       double gap1S,  double gap2S) {
    const uint64_t A = cnt[0], C = cnt[1], G = cnt[2], T = cnt[3], N = cnt[4], D = cnt[5];

    // match: AA + CC + GG + TT
    const uint64_t match =
        choose2_u64(A) + choose2_u64(C) + choose2_u64(G) + choose2_u64(T);

    // mismatch among {A,C,G,T} with different letters:
    // A*C + A*G + A*T + C*G + C*T + G*T
    const uint64_t mismatch =
        A*C + A*G + A*T + C*G + C*T + G*T;

    // gap1: '-' with (A/C/G/T/N)
    const uint64_t gap1 = (A + C + G + T + N) * D;

    // gap2: gap-gap + N-N + N-(A/C/G/T)
    const uint64_t gap2 =
        choose2_u64(D) + choose2_u64(N) + (A + C + G + T) * N;

    return match * matchS + mismatch * mismatchS + gap1 * gap1S + gap2 * gap2S;
}

std::vector<std::string> load_fasta(const std::string& path) {
    gzFile fp = gzopen(path.c_str(), "rb");
    if (!fp) {
        throw std::runtime_error("Cannot open: " + path);
    }

    kseq_t* ks = kseq_init(fp);
    std::vector<std::string> seqs;
    seqs.reserve(1024);

    while (kseq_read(ks) >= 0) {
        std::string s;
        s.resize(ks->seq.l);
        for (size_t i = 0; i < ks->seq.l; ++i) {
            s[i] = normalize_base(ks->seq.s[i]);
        }
        seqs.emplace_back(std::move(s));
    }

    kseq_destroy(ks);
    gzclose(fp);

    return seqs;
}

SPScoreResult calculate_sp_score(const std::vector<std::string>& seqs,
                                  double matchS,
                                  double mismatchS,
                                  double gap1S,
                                  double gap2S) {
    if (seqs.size() < 2) {
        throw std::invalid_argument("Error: need at least 2 sequences.");
    }

    const size_t M = seqs.size();
    const size_t L = seqs[0].size();
    if (L == 0) {
        throw std::invalid_argument("Error: empty sequences.");
    }

    for (size_t i = 1; i < M; ++i) {
        if (seqs[i].size() != L) {
            throw std::invalid_argument(
                "Error: sequences are not the same length (MSA required). "
                "seq[0] length=" + std::to_string(L) +
                ", seq[" + std::to_string(i) + "] length=" + std::to_string(seqs[i].size()));
        }
    }

    // Parallel over columns
    double total = 0.0;

    #pragma omp parallel for reduction(+:total) schedule(static)
    for (long long j = 0; j < (long long)L; ++j) {
        uint32_t cnt[6] = {0,0,0,0,0,0};

        // count column j
        for (size_t i = 0; i < M; ++i) {
            const char c = seqs[i][(size_t)j];
            // c is in {'A','C','G','T','N','-'} after normalization
            switch (c) {
                case 'A': cnt[0]++; break;
                case 'C': cnt[1]++; break;
                case 'G': cnt[2]++; break;
                case 'T': cnt[3]++; break;
                case 'N': cnt[4]++; break;
                case '-': cnt[5]++; break;
                default:  cnt[4]++; break; // should not happen
            }
        }

        total += score_of_counts(cnt, matchS, mismatchS, gap1S, gap2S);
    }

    // Pair count = M choose 2
    const long double pair_num_ld = (long double)M * (long double)(M - 1) / 2.0L;
    const double pair_num = (double)pair_num_ld;

    SPScoreResult result;
    result.total_sp = total;
    result.avg_sp = total / pair_num;
    result.scaled_sp = result.avg_sp / (double)L;

    return result;
}

static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " -i <msa.fasta[.gz]> [--match v] [--mismatch v] [--gap1 v] [--gap2 v]\n"
        << "Defaults: match=1 mismatch=-1 gap1=-2 gap2=0\n";
}

static bool parse_double_arg(int& idx, int argc, char** argv, double& out) {
    if (idx + 1 >= argc) return false;
    out = std::strtod(argv[idx + 1], nullptr);
    idx += 1;
    return true;
}

#ifndef NO_MAIN
int main(int argc, char** argv) {
    std::string input_path;
    double matchS = 1.0, mismatchS = -1.0, gap1S = -2.0, gap2S = 0.0;

    // CLI parsing
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (std::strcmp(a, "-i") == 0 || std::strcmp(a, "--input") == 0) {
            if (i + 1 >= argc) { usage(argv[0]); return 2; }
            input_path = argv[++i];
        } else if (std::strcmp(a, "--match") == 0) {
            if (!parse_double_arg(i, argc, argv, matchS)) { usage(argv[0]); return 2; }
        } else if (std::strcmp(a, "--mismatch") == 0) {
            if (!parse_double_arg(i, argc, argv, mismatchS)) { usage(argv[0]); return 2; }
        } else if (std::strcmp(a, "--gap1") == 0) {
            if (!parse_double_arg(i, argc, argv, gap1S)) { usage(argv[0]); return 2; }
        } else if (std::strcmp(a, "--gap2") == 0) {
            if (!parse_double_arg(i, argc, argv, gap2S)) { usage(argv[0]); return 2; }
        } else if (std::strcmp(a, "-h") == 0 || std::strcmp(a, "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown arg: " << a << "\n";
            usage(argv[0]);
            return 2;
        }
    }

    if (input_path.empty()) {
        usage(argv[0]);
        return 2;
    }

    try {
        // Load sequences
        std::vector<std::string> seqs = load_fasta(input_path);

        // Calculate SP score
        SPScoreResult result = calculate_sp_score(seqs, matchS, mismatchS, gap1S, gap2S);

        // Output results
        std::cout.setf(std::ios::fixed);
        std::cout << std::setprecision(6);
        std::cout << "SP score: " << result.total_sp << "\n";
        std::cout << "Avg SP: " << result.avg_sp << "\n";
        std::cout << "Scaled SP: " << result.scaled_sp << "\n";

    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    return 0;
}
#endif // NO_MAIN
