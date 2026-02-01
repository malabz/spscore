#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../src/spscore.h"
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>

// ============================================================================
// test_spscore.cpp（中文注释版）
// ----------------------------------------------------------------------------
// 这份文件是 spscore 的单元测试 + 小型性能测试，使用 doctest 框架。
//
// 主要覆盖点：
// 1) normalize_base：输入字符归一化（大小写、U->T、非法字符->N、gap 字符保留）。
// 2) choose2_u64：组合数 C(n,2)=n*(n-1)/2 的边界与正确性。
// 3) score_of_counts：
//    - 这是“单列（column）”SP score 的 O(1) 计算核心。
//    - 输入是 cnt[6]，分别表示该列中 A/C/G/T/N/'-' 的计数。
//    - 依据实现（见 src/spscore.cpp）：
//        match    = C(A,2)+C(C,2)+C(G,2)+C(T,2)
//        mismatch = A*C + A*G + A*T + C*G + C*T + G*T           （仅 A/C/G/T 之间）
//        gap1     = (A+C+G+T+N) * D                              （D 表示 '-'）
//        gap2     = C(D,2) + C(N,2) + (A+C+G+T)*N                （D-D、N-N、N-ATGC）
//      注意：N 与 '-' 的配对在本实现里计入 gap1（不是 gap2）。
// 4) calculate_sp_score：
//    - 对齐（MSA）要求：至少 2 条序列、长度相同、长度>0。
//    - 逐列统计 cnt[6] 并累加 score_of_counts 得到 total_sp。
//    - avg_sp    = total_sp / C(M,2) （M 为序列条数）
//    - scaled_sp = avg_sp / L        （L 为比对长度/列数）
//
// 运行方式（WSL）：
// - 推荐一键脚本：./run_tests.sh
// - 或者在 test/build 下直接执行：./spscore_test
// ============================================================================

// Helper function for floating point comparison
//
// 中文说明：
// 浮点比较不要直接用 ==，而是判断差值是否落在一个很小的 epsilon 内。
// 这里默认 epsilon=1e-9，对本项目的整数型计数推导（转 double）足够稳。
bool approx_equal(double a, double b, double epsilon = 1e-9) {
    return std::abs(a - b) < epsilon;
}

// ======================== Unit Tests ========================

TEST_SUITE("normalize_base") {
    // 目的：验证归一化逻辑对标准大写碱基的“恒等映射”。
    TEST_CASE("uppercase nucleotides") {
        CHECK(normalize_base('A') == 'A');
        CHECK(normalize_base('C') == 'C');
        CHECK(normalize_base('G') == 'G');
        CHECK(normalize_base('T') == 'T');
        CHECK(normalize_base('N') == 'N');
    }

    // 目的：验证小写输入会被提升到大写。
    TEST_CASE("lowercase nucleotides") {
        CHECK(normalize_base('a') == 'A');
        CHECK(normalize_base('c') == 'C');
        CHECK(normalize_base('g') == 'G');
        CHECK(normalize_base('t') == 'T');
        CHECK(normalize_base('n') == 'N');
    }

    // 目的：RNA 的 U 与 DNA 的 T 视为等价，归一化后都输出 'T'。
    TEST_CASE("uracil to thymine") {
        CHECK(normalize_base('U') == 'T');
        CHECK(normalize_base('u') == 'T');
    }

    // 目的：比对中的缺口符号 '-' 需要保留（不会被转成 N）。
    TEST_CASE("gap character") {
        CHECK(normalize_base('-') == '-');
    }

    // 目的：任何非 A/C/G/T/U/N/'-' 的字符（例如 IUPAC 扩展、数字、符号）
    // 在本实现中都统一视为未知 'N'。
    TEST_CASE("invalid characters become N") {
        CHECK(normalize_base('X') == 'N');
        CHECK(normalize_base('Y') == 'N');
        CHECK(normalize_base('1') == 'N');
        CHECK(normalize_base('*') == 'N');
    }
}

TEST_SUITE("choose2_u64") {
    // 目的：确保边界条件 (x<2) 时返回 0，不会产生负数/溢出。
    TEST_CASE("edge cases") {
        CHECK(choose2_u64(0) == 0);
        CHECK(choose2_u64(1) == 0);
        CHECK(choose2_u64(2) == 1);
    }

    // 目的：若 x 较小，可手算验证 C(x,2)=x*(x-1)/2。
    TEST_CASE("small values") {
        CHECK(choose2_u64(3) == 3);  // 3*2/2 = 3
        CHECK(choose2_u64(4) == 6);  // 4*3/2 = 6
        CHECK(choose2_u64(5) == 10); // 5*4/2 = 10
        CHECK(choose2_u64(10) == 45); // 10*9/2 = 45
    }

    // 目的：相对更大的数值，验证不会被 32 位截断，并且结果仍正确。
    TEST_CASE("large values") {
        CHECK(choose2_u64(100) == 4950);    // 100*99/2
        CHECK(choose2_u64(1000) == 499500); // 1000*999/2
    }
}

TEST_SUITE("score_of_counts") {
    // 本套件专门测试“单列计分”score_of_counts。
    //
    // cnt[6] 的索引约定（与实现保持一致）：
    //   cnt[0]=A, cnt[1]=C, cnt[2]=G, cnt[3]=T, cnt[4]=N, cnt[5]='-'
    //
    // 默认计分参数：match=+1, mismatch=-1, gap1=-2, gap2=0。
    // 其中 gap2=0 表示：gap-gap、N-N、N-ATGC 这三类配对默认不加分不扣分。

    TEST_CASE("all matches - identical bases") {
        uint32_t cnt[6];

        // 场景：一列里 4 条序列全是 A。
        // - A-A 配对数 = C(4,2)=6
        // - mismatch=0, gap1=0, gap2=0
        // 所以列得分 = 6 * match(=1) = 6
        cnt[0] = 4; cnt[1] = 0; cnt[2] = 0; cnt[3] = 0; cnt[4] = 0; cnt[5] = 0;
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), 6.0));
    }


    TEST_CASE("gap penalties - gap1") {
        uint32_t cnt[6];

        // 场景：2 个 A、2 个 '-'。
        // - gap1 的定义：'-' 与 (A/C/G/T/N) 的配对数
        //   => gap1 = (A+C+G+T+N) * D = (2+0+0+0+0) * 2 = 4
        // - match(A-A)=C(2,2)=1
        // - gap2=0
        // 列得分 = 1*match(1) + 4*gap1(-2) = 1 - 8 = -7
        cnt[0] = 2; cnt[1] = 0; cnt[2] = 0; cnt[3] = 0; cnt[4] = 0; cnt[5] = 2;
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), -7.0));
    }

    TEST_CASE("gap-gap interactions - gap2") {
        uint32_t cnt[6];

        // 场景：一列里全是 '-'，共 4 个。
        // - gap2 中包含 gap-gap：C(D,2)=C(4,2)=6
        // - 默认 gap2S=0，因此得分为 6*0 = 0
        cnt[0] = 0; cnt[1] = 0; cnt[2] = 0; cnt[3] = 0; cnt[4] = 0; cnt[5] = 4;
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), 0.0));
    }

    TEST_CASE("N character penalties - gap2") {
        uint32_t cnt[6];

        // 场景 1：一列里只有 2 个 N。
        // - gap2 包含 N-N：C(N,2)=C(2,2)=1
        // - gap2S=0 => 得分 0
        cnt[0] = 0; cnt[1] = 0; cnt[2] = 0; cnt[3] = 0; cnt[4] = 2; cnt[5] = 0;
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), 0.0));

        // 场景 2：2 个 A、2 个 N。
        // - match(A-A)=C(2,2)=1
        // - gap2 包含 N-N：C(2,2)=1 以及 N-ATGC：(A+C+G+T)*N = 2*2 = 4
        // - 但 gap2S=0，所以这些配对都不贡献分数
        // 期望列得分仅来自 match => 1
        cnt[0] = 2; cnt[1] = 0; cnt[2] = 0; cnt[3] = 0; cnt[4] = 2; cnt[5] = 0;
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), 1.0));
    }

    TEST_CASE("mixed case") {
        uint32_t cnt[6];

        // 场景：一列里 A/C/G/T 各 1 个。
        // - match=0（没有同字母）
        // - mismatch = A*C + A*G + A*T + C*G + C*T + G*T = 6
        // 列得分 = 6 * (-1) = -6
        cnt[0] = 1; cnt[1] = 1; cnt[2] = 1; cnt[3] = 1; cnt[4] = 0; cnt[5] = 0;
        CHECK(approx_equal(score_of_counts(cnt, 1.0, -1.0, -2.0, 0.0), -6.0));
    }

    TEST_CASE("custom scoring scheme") {
        uint32_t cnt[6];

        // 目的：验证自定义 match 分值会线性影响最终列得分。
        // 场景：4 个 A -> match=C(4,2)=6；matchS=2.0 => 得分=12
        cnt[0] = 4; cnt[1] = 0; cnt[2] = 0; cnt[3] = 0; cnt[4] = 0; cnt[5] = 0;
        CHECK(approx_equal(score_of_counts(cnt, 2.0, -1.0, -2.0, 0.0), 12.0));
    }
}

TEST_SUITE("calculate_sp_score") {
    // 本套件验证总接口 calculate_sp_score 的整体行为：
    // - total_sp：所有列的 score_of_counts 之和
    // - avg_sp：除以序列对数 C(M,2)
    // - scaled_sp：再除以长度 L

    TEST_CASE("minimal case - 2 identical sequences") {
        std::vector<std::string> seqs = {"ACGT", "ACGT"};
        auto result = calculate_sp_score(seqs);

        // M=2 => C(M,2)=1 对；L=4。
        // 每列都是 match（+1），总分 total_sp=4。
        // avg_sp=4/1=4；scaled_sp=4/4=1（每对每列平均 +1）。
        CHECK(approx_equal(result.total_sp, 4.0));
        CHECK(approx_equal(result.avg_sp, 4.0));
        CHECK(approx_equal(result.scaled_sp, 1.0));
    }

    TEST_CASE("2 completely different sequences") {
        std::vector<std::string> seqs = {"AAAA", "CCCC"};
        auto result = calculate_sp_score(seqs);

        // M=2 => 只有 1 对；L=4。
        // 每列都是 mismatch（-1），total_sp = 4 * (-1) = -4。
        // avg_sp=-4；scaled_sp=-4/4=-1。
        CHECK(approx_equal(result.total_sp, -4.0));
        CHECK(approx_equal(result.avg_sp, -4.0));
        CHECK(approx_equal(result.scaled_sp, -1.0));
    }

    TEST_CASE("3 identical sequences") {
        std::vector<std::string> seqs = {"AC", "AC", "AC"};
        auto result = calculate_sp_score(seqs);

        // M=3 => pair 数 C(3,2)=3；L=2。
        // 列 0：A,A,A => match=C(3,2)=3
        // 列 1：C,C,C => match=C(3,2)=3
        // total_sp=6；avg_sp=6/3=2；scaled_sp=2/2=1。
        CHECK(approx_equal(result.total_sp, 6.0));
        CHECK(approx_equal(result.avg_sp, 2.0));
        CHECK(approx_equal(result.scaled_sp, 1.0));
    }

    TEST_CASE("sequences with gaps") {
        std::vector<std::string> seqs = {"A-", "A-"};
        auto result = calculate_sp_score(seqs);

        // 只有 1 对；2 列。
        // 列 0：A-A => match=+1
        // 列 1：'-' '-' => gap2 中的 gap-gap，但 gap2S=0，所以这一列贡献 0
        // total_sp=1
        CHECK(approx_equal(result.total_sp, 1.0));
    }

    TEST_CASE("sequences with gap penalties") {
        std::vector<std::string> seqs = {"AA", "A-"};
        auto result = calculate_sp_score(seqs);

        // 只有 1 对；2 列。
        // 列 0：A-A => +1
        // 列 1：A- => gap1（-2）
        // total_sp = 1 + (-2) = -1
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

        // M=4 => pair 数 C(4,2)=6；L=4。
        // 每列 4 个相同碱基：列得分=C(4,2)=6。
        // total_sp = 4 列 * 6 = 24
        // avg_sp = 24/6=4；scaled_sp = 4/4=1。
        CHECK(approx_equal(result.total_sp, 24.0));
        CHECK(approx_equal(result.avg_sp, 4.0));
        CHECK(approx_equal(result.scaled_sp, 1.0));
    }

    TEST_CASE("mixed quality alignment") {
        std::vector<std::string> seqs = {
            "ACGT",
            "ACCT",
            "ACGT"
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

// 中文说明：
// 这里的“性能测试”本质上是带阈值的粗粒度时间断言（类似 smoke benchmark）。
// - 在不同机器/编译器/负载/是否启用 OpenMP 的情况下，耗时波动可能较大。
// - 如果你在 CI 或较慢的环境中遇到偶发失败，可以考虑：
//   1) 调整阈值；或
//   2) 将这些用例单独标记/拆分（doctest 支持按 suite/case 过滤运行）。
// 当前我们保持原阈值不变，仅补充解释。

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

// 中文说明：
// 这一组更像“行为/方向性”测试：不追求精确分值，而是验证结果的趋势是否符合直觉：
// - 高度保守（conserved）应接近满分（scaled_sp 接近 1）
// - 高度可变（variable）应偏负分
// - 插入缺失（indel）会受到 gap1 惩罚而拉低总分
// - SNP 仅局部影响，总体仍应偏正
// - 含 U 的 RNA 与 T 的 DNA 在 normalize 后应等价
// - N（未知碱基）按 gap2 规则处理，会降低相对于完美匹配的得分

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

