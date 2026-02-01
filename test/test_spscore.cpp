#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../src/spscore.h"
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>

// Helper function for floating point comparison
bool approx_equal(double a, double b, double epsilon = 1e-9) {
    return std::abs(a - b) < epsilon;
}

// ======================== Unit Tests ========================

TEST_SUITE("normalize_base") {
    TEST_CASE("uppercase nucleotides") {
        CHECK(normalize_base('A') == 'A');
        CHECK(normalize_base('C') == 'C');
        CHECK(normalize_base('G') == 'G');
        CHECK(normalize_base('T') == 'T');
        CHECK(normalize_base('N') == 'N');
    }

    TEST_CASE("lowercase nucleotides") {
        CHECK(normalize_base('a') == 'A');
        CHECK(normalize_base('c') == 'C');
        CHECK(normalize_base('g') == 'G');
        CHECK(normalize_base('t') == 'T');
        CHECK(normalize_base('n') == 'N');
    }

    TEST_CASE("uracil to thymine") {
        CHECK(normalize_base('U') == 'T');
        CHECK(normalize_base('u') == 'T');
    }

    TEST_CASE("gap character") {
        CHECK(normalize_base('-') == '-');
    }

    TEST_CASE("invalid characters become N") {
        CHECK(normalize_base('X') == 'N');
        CHECK(normalize_base('Y') == 'N');
        CHECK(normalize_base('1') == 'N');
        CHECK(normalize_base('*') == 'N');
    }
}

TEST_SUITE("choose2_u64") {
    TEST_CASE("edge cases") {
        CHECK(choose2_u64(0) == 0);
        CHECK(choose2_u64(1) == 0);
        CHECK(choose2_u64(2) == 1);
    }

    TEST_CASE("small values") {
        CHECK(choose2_u64(3) == 3);  // 3*2/2 = 3
        CHECK(choose2_u64(4) == 6);  // 4*3/2 = 6
        CHECK(choose2_u64(5) == 10); // 5*4/2 = 10
        CHECK(choose2_u64(10) == 45); // 10*9/2 = 45
    }

    TEST_CASE("large values") {
        CHECK(choose2_u64(100) == 4950);    // 100*99/2
        CHECK(choose2_u64(1000) == 499500); // 1000*999/2
    }
}

TEST_SUITE("score_of_counts") {
    TEST_CASE("all matches - identical bases") {
        uint32_t cnt[6];

        // All A's: 4 sequences, all A
        cnt[0] = 4; cnt[1] = 0; cnt[2] = 0; cnt[3] = 0; cnt[4] = 0; cnt[5] = 0;
        // 4 choose 2 = 6 matches, score = 6 * 1 = 6
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), 6.0));
    }

    TEST_CASE("all mismatches") {
        uint32_t cnt[6];

        // 2 A's, 2 C's: 2*2 = 4 mismatches
        cnt[0] = 2; cnt[1] = 2; cnt[2] = 0; cnt[3] = 0; cnt[4] = 0; cnt[5] = 0;
        // score = 4 * (-1) = -4
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), -4.0));
    }

    TEST_CASE("gap penalties - gap1") {
        uint32_t cnt[6];

        // 2 A's, 2 gaps: (2+0+0+0+0) * 2 = 4 gap1 interactions
        cnt[0] = 2; cnt[1] = 0; cnt[2] = 0; cnt[3] = 0; cnt[4] = 0; cnt[5] = 2;
        // 1 match (A-A) + 4 gap1 = 1*1 + 4*(-2) = 1 - 8 = -7
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), -7.0));
    }

    TEST_CASE("gap-gap interactions - gap2") {
        uint32_t cnt[6];

        // 4 gaps: choose2(4) = 6 gap-gap interactions
        cnt[0] = 0; cnt[1] = 0; cnt[2] = 0; cnt[3] = 0; cnt[4] = 0; cnt[5] = 4;
        // score = 6 * 0 = 0 (gap2 score is 0 by default)
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), 0.0));
    }

    TEST_CASE("N character penalties - gap2") {
        uint32_t cnt[6];

        // 2 N's: choose2(2) = 1 N-N interaction
        cnt[0] = 0; cnt[1] = 0; cnt[2] = 0; cnt[3] = 0; cnt[4] = 2; cnt[5] = 0;
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), 0.0));

        // 2 A's, 2 N's: (2+0+0+0)*2 = 4 A-N interactions (gap2)
        cnt[0] = 2; cnt[1] = 0; cnt[2] = 0; cnt[3] = 0; cnt[4] = 2; cnt[5] = 0;
        // 1 match (A-A) + 4 A-N (gap2) + 1 N-N = 1 + 0 + 0 = 1
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), 1.0));
    }

    TEST_CASE("mixed case") {
        uint32_t cnt[6];

        // 1A, 1C, 1G, 1T
        cnt[0] = 1; cnt[1] = 1; cnt[2] = 1; cnt[3] = 1; cnt[4] = 0; cnt[5] = 0;
        // All pairs are mismatches: A-C, A-G, A-T, C-G, C-T, G-T = 6 mismatches
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), -6.0));
    }

    TEST_CASE("custom scoring scheme") {
        uint32_t cnt[6];

        // 4 A's
        cnt[0] = 4; cnt[1] = 0; cnt[2] = 0; cnt[3] = 0; cnt[4] = 0; cnt[5] = 0;
        // 6 matches with score 2.0
        CHECK(approx_equal(score_of_counts(cnt, 2.0, -1.0, -2.0, 0.0), 12.0));
    }
}

TEST_SUITE("calculate_sp_score") {
    TEST_CASE("minimal case - 2 identical sequences") {
        std::vector<std::string> seqs = {"ACGT", "ACGT"};
        auto result = calculate_sp_score(seqs);

        // 1 pair, 4 columns, all matches
        // Total = 4 * 1 (match) = 4
        // Avg = 4 / 1 = 4
        // Scaled = 4 / 4 = 1
        CHECK(approx_equal(result.total_sp, 4.0));
        CHECK(approx_equal(result.avg_sp, 4.0));
        CHECK(approx_equal(result.scaled_sp, 1.0));
    }

    TEST_CASE("2 completely different sequences") {
        std::vector<std::string> seqs = {"AAAA", "CCCC"};
        auto result = calculate_sp_score(seqs);

        // 1 pair, 4 columns, all mismatches
        // Total = 4 * (-1) = -4
        CHECK(approx_equal(result.total_sp, -4.0));
        CHECK(approx_equal(result.avg_sp, -4.0));
        CHECK(approx_equal(result.scaled_sp, -1.0));
    }

    TEST_CASE("3 identical sequences") {
        std::vector<std::string> seqs = {"AC", "AC", "AC"};
        auto result = calculate_sp_score(seqs);

        // 3 pairs (3 choose 2), 2 columns, all matches
        // Column 1: 3 A's -> 3 matches, score = 3
        // Column 2: 3 C's -> 3 matches, score = 3
        // Total = 6
        // Avg = 6 / 3 = 2
        // Scaled = 2 / 2 = 1
        CHECK(approx_equal(result.total_sp, 6.0));
        CHECK(approx_equal(result.avg_sp, 2.0));
        CHECK(approx_equal(result.scaled_sp, 1.0));
    }

    TEST_CASE("sequences with gaps") {
        std::vector<std::string> seqs = {"A-", "A-"};
        auto result = calculate_sp_score(seqs);

        // 1 pair, 2 columns
        // Column 1: 2 A's -> 1 match = 1
        // Column 2: 2 gaps -> 1 gap-gap (gap2) = 0
        // Total = 1
        CHECK(approx_equal(result.total_sp, 1.0));
    }

    TEST_CASE("sequences with gap penalties") {
        std::vector<std::string> seqs = {"AA", "A-"};
        auto result = calculate_sp_score(seqs);

        // 1 pair, 2 columns
        // Column 1: 2 A's -> 1 match = 1
        // Column 2: 1 A, 1 gap -> 1 gap1 = -2
        // Total = -1
        CHECK(approx_equal(result.total_sp, -1.0));
    }

    TEST_CASE("larger alignment") {
        std::vector<std::string> seqs = {
            "ACGT",
            "ACGT",
            "ACGT",
            "ACGT"
        };
        auto result = calculate_sp_score(seqs);

        // 4 sequences: 6 pairs (4 choose 2)
        // 4 columns, each column has 4 identical bases
        // Each column: 4 choose 2 = 6 matches
        // Total = 4 * 6 = 24
        // Avg = 24 / 6 = 4
        // Scaled = 4 / 4 = 1
        CHECK(approx_equal(result.total_sp, 24.0));
        CHECK(approx_equal(result.avg_sp, 4.0));
        CHECK(approx_equal(result.scaled_sp, 1.0));
    }

    TEST_CASE("mixed quality alignment") {
        std::vector<std::string> seqs = {
            "ACGT",
            "ACCT",
            "AGGT"
        };
        auto result = calculate_sp_score(seqs);

        // 3 sequences: 3 pairs
        // Column 0: AAA -> 3 matches = 3
        // Column 1: CCC -> 3 matches = 3
        // Column 2: GCG -> 1 match (GG) + 2 mismatches (G-C, G-C) = 1 - 2 = -1
        // Column 3: TTT -> 3 matches = 3
        // Total = 3 + 3 - 1 + 3 = 8
        CHECK(approx_equal(result.total_sp, 8.0));
    }

    TEST_CASE("custom scoring parameters") {
        std::vector<std::string> seqs = {"AA", "CC"};
        auto result = calculate_sp_score(seqs, 2.0, -3.0, -4.0, -1.0);

        // 1 pair, 2 columns, all mismatches
        // Total = 2 * (-3) = -6
        CHECK(approx_equal(result.total_sp, -6.0));
    }
}

TEST_SUITE("edge cases and error handling") {
    TEST_CASE("single sequence - should throw") {
        std::vector<std::string> seqs = {"ACGT"};
        CHECK_THROWS_AS(calculate_sp_score(seqs), std::invalid_argument);
    }

    TEST_CASE("empty sequence vector - should throw") {
        std::vector<std::string> seqs;
        CHECK_THROWS_AS(calculate_sp_score(seqs), std::invalid_argument);
    }

    TEST_CASE("empty sequences - should throw") {
        std::vector<std::string> seqs = {"", ""};
        CHECK_THROWS_AS(calculate_sp_score(seqs), std::invalid_argument);
    }

    TEST_CASE("unequal length sequences - should throw") {
        std::vector<std::string> seqs = {"ACG", "ACGT"};
        CHECK_THROWS_AS(calculate_sp_score(seqs), std::invalid_argument);
    }

    TEST_CASE("very long sequences") {
        std::string seq1(10000, 'A');
        std::string seq2(10000, 'A');
        std::vector<std::string> seqs = {seq1, seq2};

        auto result = calculate_sp_score(seqs);
        // 1 pair, 10000 columns, all matches
        CHECK(approx_equal(result.total_sp, 10000.0));
    }

    TEST_CASE("many sequences") {
        std::vector<std::string> seqs;
        for (int i = 0; i < 100; ++i) {
            seqs.push_back("ACGT");
        }

        auto result = calculate_sp_score(seqs);
        // 100 sequences: 100*99/2 = 4950 pairs
        // 4 columns, each with 100 identical bases
        // Each column: 100 choose 2 = 4950 matches
        // Total = 4 * 4950 = 19800
        CHECK(approx_equal(result.total_sp, 19800.0));
    }

    TEST_CASE("all gap sequences") {
        std::vector<std::string> seqs = {"---", "---"};
        auto result = calculate_sp_score(seqs);

        // 1 pair, 3 columns of gaps
        // Each column: 2 gaps -> 1 gap-gap (gap2 = 0)
        CHECK(approx_equal(result.total_sp, 0.0));
    }

    TEST_CASE("all N sequences") {
        std::vector<std::string> seqs = {"NNN", "NNN"};
        auto result = calculate_sp_score(seqs);

        // 1 pair, 3 columns of N's
        // Each column: 2 N's -> 1 N-N (gap2 = 0)
        CHECK(approx_equal(result.total_sp, 0.0));
    }

    TEST_CASE("mixed valid/invalid characters") {
        // Invalid characters should be converted to N
        std::vector<std::string> seqs = {"AXY", "AXY"};
        // After normalization: ANN, ANN
        auto result = calculate_sp_score(seqs);

        // Column 0: AA -> 1 match = 1
        // Column 1: NN -> 1 N-N (gap2 = 0) = 0
        // Column 2: NN -> 1 N-N (gap2 = 0) = 0
        // Total = 1
        CHECK(approx_equal(result.total_sp, 1.0));
    }
}

// ======================== Performance Tests ========================

TEST_SUITE("performance tests") {
    TEST_CASE("small alignment (10 sequences x 100 bp)") {
        std::vector<std::string> seqs;
        for (int i = 0; i < 10; ++i) {
            std::string seq;
            for (int j = 0; j < 100; ++j) {
                seq += "ACGT"[j % 4];
            }
            seqs.push_back(seq);
        }

        auto start = std::chrono::high_resolution_clock::now();
        auto result = calculate_sp_score(seqs);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "10 sequences x 100 bp: " << duration.count() << " μs\n";

        CHECK(duration.count() < 100000); // Should be < 100ms
    }

    TEST_CASE("medium alignment (50 sequences x 500 bp)") {
        std::vector<std::string> seqs;
        for (int i = 0; i < 50; ++i) {
            std::string seq;
            for (int j = 0; j < 500; ++j) {
                seq += "ACGT"[j % 4];
            }
            seqs.push_back(seq);
        }

        auto start = std::chrono::high_resolution_clock::now();
        auto result = calculate_sp_score(seqs);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "50 sequences x 500 bp: " << duration.count() << " ms\n";

        CHECK(duration.count() < 1000); // Should be < 1s
    }

    TEST_CASE("large alignment (100 sequences x 1000 bp)") {
        std::vector<std::string> seqs;
        for (int i = 0; i < 100; ++i) {
            std::string seq;
            for (int j = 0; j < 1000; ++j) {
                seq += "ACGT"[j % 4];
            }
            seqs.push_back(seq);
        }

        auto start = std::chrono::high_resolution_clock::now();
        auto result = calculate_sp_score(seqs);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "100 sequences x 1000 bp: " << duration.count() << " ms\n";

        CHECK(duration.count() < 5000); // Should be < 5s
    }

    TEST_CASE("very large alignment (200 sequences x 2000 bp)") {
        std::vector<std::string> seqs;
        for (int i = 0; i < 200; ++i) {
            std::string seq;
            for (int j = 0; j < 2000; ++j) {
                seq += "ACGT"[j % 4];
            }
            seqs.push_back(seq);
        }

        auto start = std::chrono::high_resolution_clock::now();
        auto result = calculate_sp_score(seqs);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "200 sequences x 2000 bp: " << duration.count() << " ms\n";

        CHECK(duration.count() < 20000); // Should be < 20s
    }

    TEST_CASE("long sequences (10 sequences x 10000 bp)") {
        std::vector<std::string> seqs;
        for (int i = 0; i < 10; ++i) {
            std::string seq;
            for (int j = 0; j < 10000; ++j) {
                seq += "ACGT"[j % 4];
            }
            seqs.push_back(seq);
        }

        auto start = std::chrono::high_resolution_clock::now();
        auto result = calculate_sp_score(seqs);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "10 sequences x 10000 bp: " << duration.count() << " ms\n";

        CHECK(duration.count() < 1000); // Should be < 1s
    }

    TEST_CASE("many sequences (500 sequences x 100 bp)") {
        std::vector<std::string> seqs;
        for (int i = 0; i < 500; ++i) {
            std::string seq;
            for (int j = 0; j < 100; ++j) {
                seq += "ACGT"[j % 4];
            }
            seqs.push_back(seq);
        }

        auto start = std::chrono::high_resolution_clock::now();
        auto result = calculate_sp_score(seqs);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "500 sequences x 100 bp: " << duration.count() << " ms\n";

        CHECK(duration.count() < 5000); // Should be < 5s
    }
}

// ======================== Realistic Biological Scenarios ========================

TEST_SUITE("biological scenarios") {
    TEST_CASE("conserved region") {
        std::vector<std::string> seqs = {
            "ATGCATGC",
            "ATGCATGC",
            "ATGCATGC",
            "ATGCATGC"
        };
        auto result = calculate_sp_score(seqs);

        // Perfect conservation should give maximum score
        CHECK(result.scaled_sp > 0.99); // Very close to 1.0
    }

    TEST_CASE("variable region") {
        std::vector<std::string> seqs = {
            "AAAACCCC",
            "CCCCAAAA",
            "GGGGTTTT",
            "TTTTGGGG"
        };
        auto result = calculate_sp_score(seqs);

        // High variability should give negative score
        CHECK(result.scaled_sp < 0);
    }

    TEST_CASE("insertion/deletion variation") {
        std::vector<std::string> seqs = {
            "ATGC--GC",
            "ATGCAT--",
            "ATGCATGC",
            "AT--ATGC"
        };
        auto result = calculate_sp_score(seqs);

        // Mixed gaps should give intermediate score
        CHECK(result.total_sp < 0); // Due to gap penalties
    }

    TEST_CASE("SNP variation") {
        std::vector<std::string> seqs = {
            "ATGCATGC",
            "ATGCATGC",
            "ATGCCTGC",  // Single nucleotide change
            "ATGCATGC"
        };
        auto result = calculate_sp_score(seqs);

        // Mostly conserved with one SNP
        CHECK(result.scaled_sp > 0.5); // Still mostly positive
    }

    TEST_CASE("RNA sequences with uracil") {
        // U should be converted to T
        std::vector<std::string> seqs = {
            "AUGC",  // RNA
            "ATGC"   // DNA
        };
        auto result = calculate_sp_score(seqs);

        // Should be treated as identical after normalization
        CHECK(approx_equal(result.total_sp, 4.0));
    }

    TEST_CASE("ambiguous bases") {
        std::vector<std::string> seqs = {
            "ANNGC",
            "ACNGC",
            "ACGNC"
        };
        auto result = calculate_sp_score(seqs);

        // N's should be penalized
        CHECK(result.total_sp < 9.0); // Less than if all were perfect matches
    }
}

