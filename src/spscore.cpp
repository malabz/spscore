// sp_omp_kseq.cpp
// Requires: kseq.h (in include path) + zlib
//
// Build:
//   g++ -O3 -march=native -fopenmp sp_omp_kseq.cpp -lz -o sp
//
// Run:
//   ./sp -i alignment.fa[.gz] [--match v] [--mismatch v] [--gap1 v] [--gap2 v]
//
// ============================================================================
// 中文说明（实现思路概览）
// ----------------------------------------------------------------------------
// 本文件实现“SP score（Sum-of-Pairs score）”的计算：
// 1) 输入是一组已经对齐的序列（MSA，多序列比对），要求所有序列长度一致。
// 2) 逐列（column）统计该列中 A/C/G/T/N/'-' 的个数。
// 3) 对于该列，所有序列两两配对（pair）会产生：匹配(match)、错配(mismatch)、
//    缺口相关惩罚(gap1/gap2)。本实现用计数公式一次性算出该列所有 pair 的贡献，
//    避免 O(M^2) 的逐对比较。
// 4) 所有列的分数累加得到 total_sp；再除以 pair 数得到 avg_sp；
//    再除以列数得到 scaled_sp。
//
// 计数约定：cnt[6] = [A,C,G,T,N,'-']。
// - match：AA + CC + GG + TT（同碱基的配对数）
// - mismatch：仅在 {A,C,G,T} 之间统计不同碱基的配对数（例如 A-C、A-G 等）
// - gap1：'-' 与 (A/C/G/T/N) 的配对数（碱基/未知 与 缺口）
// - gap2：缺口-缺口 + N-N + N-(A/C/G/T) 的配对数（“第二类缺口/未知”惩罚）
//   说明：gap2 的定义在不同工具/论文里可能略有差异，本实现按当前单元测试/历史
//   行为来定义。
// ============================================================================

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

// Your local kseq.h
#include "kseq.h"

// kseq needs this macro instantiation for gzFile:
KSEQ_INIT(gzFile, gzread)

// Normalize to {'A','C','G','T','N','-'} with U->T and invalid->N.
//
// 中文说明：
// FASTA 中可能出现小写、U（RNA）、以及各种非标准字符。
// 这里统一归一化到 6 类字符：A/C/G/T/N/-
// - U 被当作 T
// - 任何无法识别的字符一律当作 N（未知）
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

// 计算组合数 C(x,2) = x*(x-1)/2，用于“同类元素两两配对”的数量。
// 例如某列 A 出现 A 次，则 A-A 的配对数是 C(A,2)。
uint64_t choose2_u64(uint64_t x) {
    return (x < 2) ? 0ULL : (x * (x - 1)) / 2ULL;
}

// cnt: [A,C,G,T,N,'-'] in indices 0..5
//
// 中文说明：
// 根据一列中各字符的计数，直接计算该列所有序列两两配对的总得分。
// 这样做的复杂度是 O(1)，比起对 M 条序列做 O(M^2) 的逐对比较快很多。
//
// 计分参数：
// - matchS:   匹配得分（例如 +1）
// - mismatchS:错配得分（例如 -1）
// - gap1S:    gap1 惩罚（例如 -2）
// - gap2S:    gap2 惩罚（例如 0 或其他）
double score_of_counts(const uint32_t cnt[6],
                       double matchS, double mismatchS,
                       double gap1S,  double gap2S) {
    // 为避免计数乘法溢出，转成 64 位
    const uint64_t A = cnt[0], C = cnt[1], G = cnt[2], T = cnt[3], N = cnt[4], D = cnt[5];

    // match: AA + CC + GG + TT
    // 同碱基两两配对：C(A,2)+C(C,2)+...
    const uint64_t match =
        choose2_u64(A) + choose2_u64(C) + choose2_u64(G) + choose2_u64(T);

    // mismatch among {A,C,G,T} with different letters:
    // A*C + A*G + A*T + C*G + C*T + G*T
    // 只统计明确碱基 A/C/G/T 之间的互不相同配对。
    const uint64_t mismatch =
        A*C + A*G + A*T + C*G + C*T + G*T;

    // gap1: '-' with (A/C/G/T/N)
    // 缺口 '-' 与“非缺口字符（包括 N）”的配对数。
    const uint64_t gap1 = (A + C + G + T + N) * D;

    // gap2: gap-gap + N-N + N-(A/C/G/T)
    // gap-gap：C(D,2)
    // N-N：    C(N,2)
    // N-ATGC： (A+C+G+T)*N
    // 注意：这里没有把 N 与 '-' 的配对放在 gap2 中，而是归入 gap1。
    const uint64_t gap2 =
        choose2_u64(D) + choose2_u64(N) + (A + C + G + T) * N;

    // 将各类 pair 数量乘以对应权重，得到该列 SP score 总贡献。
    return match * matchS + mismatch * mismatchS + gap1 * gap1S + gap2 * gap2S;
}

// 读取 fasta（可 gzip 压缩）并把每条序列归一化。
// 中文说明：
// - 本函数读取所有序列到内存中，适合中等规模数据。
// - 返回 vector<string>，每个 string 是一条已经做过 normalize_base 的序列。
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
    // 基本合法性检查：至少两条序列才有 pair。
    if (seqs.size() < 2) {
        throw std::invalid_argument("Error: need at least 2 sequences.");
    }

    const size_t M = seqs.size();
    const size_t L = seqs[0].size();
    if (L == 0) {
        throw std::invalid_argument("Error: empty sequences.");
    }

    // 中文说明：
    // SP score 的计算要求输入是 MSA：所有序列等长（已经对齐）。
    for (size_t i = 1; i < M; ++i) {
        if (seqs[i].size() != L) {
            throw std::invalid_argument(
                "Error: sequences are not the same length (MSA required). "
                "seq[0] length=" + std::to_string(L) +
                ", seq[" + std::to_string(i) + "] length=" + std::to_string(seqs[i].size()));
        }
    }

    // Parallel over columns
    // 中文说明：
    // 各列之间相互独立，所以可以用 OpenMP 做并行归约求和。
    // total 是全局总分，使用 reduction(+:total) 来避免数据竞争。
    double total = 0.0;

    for (long long j = 0; j < (long long)L; ++j) {
        // cnt 存放该列的字符统计：A/C/G/T/N/'-'
        uint32_t cnt[6] = {0,0,0,0,0,0};

        // count column j
        for (size_t i = 0; i < M; ++i) {
            const char c = seqs[i][(size_t)j];
            // 考虑归一化逻辑，确保 U/u 被视为 T，并处理可能的小写或非法字符
            switch (c) {
                case 'A': case 'a': cnt[0]++; break;
                case 'C': case 'c': cnt[1]++; break;
                case 'G': case 'g': cnt[2]++; break;
                case 'T': case 't':
                case 'U': case 'u': cnt[3]++; break;
                case 'N': case 'n': cnt[4]++; break;
                case '-':           cnt[5]++; break;
                default:            cnt[4]++; break; // 非标准字符映射为 N
            }
        }

        // 该列的 pair 总分加入总分
        total += score_of_counts(cnt, matchS, mismatchS, gap1S, gap2S);
    }

    // Pair count = M choose 2
    // 中文说明：
    // avg_sp 表示“平均每一对序列（pair）在全比对上的得分”。
    // pair_num = C(M,2)
    const long double pair_num_ld = (long double)M * (long double)(M - 1) / 2.0L;
    const double pair_num = (double)pair_num_ld;

    SPScoreResult result;
    result.total_sp = total;
    result.avg_sp = total / pair_num;

    // scaled_sp：再除以列数 L，得到“每个 pair、每列”平均得分，便于不同长度对齐间比较。
    result.scaled_sp = result.avg_sp / (double)L;

    return result;
}

// 流式计算 SP score（内存友好版本）
//
// 中文说明：
// 采用单遍扫描策略：
// 1. 流式读取每条序列，累积统计每列的 A/C/G/T/N/'-' 计数
// 2. 读完所有序列后，基于各列的计数一次性计算总得分
//
// 关键优化：
// - 只扫描文件一次
// - 内存占用 O(L)，L 是序列长度，与序列数量 M 无关
// - 每条序列读入后仅更新计数器，然后立即丢弃，不保存在内存中
//
// 时间复杂度：O(M*L)
SPScoreResult calculate_sp_score_streaming(const std::string& path,
                                            double matchS,
                                            double mismatchS,
                                            double gap1S,
                                            double gap2S) {
    gzFile fp = gzopen(path.c_str(), "rb");
    if (!fp) {
        throw std::runtime_error("Cannot open: " + path);
    }

    kseq_t* ks = kseq_init(fp);

    size_t M = 0;  // 已读序列数量
    size_t L = 0;  // 序列长度

    // 列式计数器：column_counts[j][0..5] 记录第 j 列所有序列的 A/C/G/T/N/'-' 计数
    std::vector<std::vector<uint32_t>> column_counts;

    // ========== 第一阶段：流式读取并累积各列计数 ==========
    size_t progress_interval = 1000;  // Update progress every 1000 sequences

    while (kseq_read(ks) >= 0) {
        const size_t seq_len = ks->seq.l;

        // 第一条序列：初始化
        if (M == 0) {
            L = seq_len;
            if (L == 0) {
                kseq_destroy(ks);
                gzclose(fp);
                throw std::invalid_argument("Error: empty sequences.");
            }
            // 初始化列计数器
            column_counts.resize(L, std::vector<uint32_t>(6, 0));
        } else {
            // 验证序列长度一致
            if (seq_len != L) {
                kseq_destroy(ks);
                gzclose(fp);
                throw std::invalid_argument(
                    "Error: sequences are not the same length (MSA required). "
                    "seq[0] length=" + std::to_string(L) +
                    ", seq[" + std::to_string(M) + "] length=" + std::to_string(seq_len));
            }
        }

        // 更新各列的碱基计数
        for (size_t j = 0; j < L; ++j) {
            const char c = normalize_base(ks->seq.s[j]);
            switch (c) {
                case 'A': column_counts[j][0]++; break;
                case 'C': column_counts[j][1]++; break;
                case 'G': column_counts[j][2]++; break;
                case 'T': column_counts[j][3]++; break;
                case 'N': column_counts[j][4]++; break;
                case '-': column_counts[j][5]++; break;
                default:  column_counts[j][4]++; break; // 其他字符视为 N
            }
        }

        ++M;

        // Real-time progress (common style): \r overwrite + flush
        if (M == 1 || (M % progress_interval == 0)) {
            static const char spinner[] = "|/-\\";
            static size_t spin_idx = 0;

            std::cerr
                << '\r'                                   // 回到行首覆盖
                << spinner[spin_idx++ & 3]                // 小转轮
                << " Reading: " << M << " sequences"
                << " (length: " << L << " bp)    "        // 末尾留空，覆盖残留字符
                << std::flush;
        }

    }

    std::cerr
        << "\rDone: " << M << " sequences loaded (length: " << L << " bp)          \n"
        << std::flush;


    kseq_destroy(ks);
    gzclose(fp);

    if (M < 2) {
        throw std::invalid_argument("Error: need at least 2 sequences.");
    }

    // ========== 第二阶段：基于各列计数计算总得分（可并行）==========
    double total = 0.0;

    for (long long j = 0; j < (long long)L; ++j) {
        total += score_of_counts(column_counts[j].data(), matchS, mismatchS, gap1S, gap2S);
    }
    // ========== 计算 avg 和 scaled ==========
    const long double pair_num_ld = (long double)M * (long double)(M - 1) / 2.0L;
    const double pair_num = (double)pair_num_ld;

    SPScoreResult result;
    result.total_sp = total;
    result.avg_sp = total / pair_num;
    result.scaled_sp = result.avg_sp / (double)L;

    return result;
}

// CLI 使用说明
static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " -i <msa.fasta[.gz]> [options]\n"
        << "\n"
        << "Required arguments:\n"
        << "  -i, --input FILE     Input FASTA file (supports .gz compression)\n"
        << "\n"
        << "Optional arguments:\n"
        << "  --match SCORE        Match score (default: 1.0)\n"
        << "  --mismatch SCORE     Mismatch score (default: -1.0)\n"
        << "  --gap1 SCORE         gap1 penalty (default: -2.0)\n"
        << "                       gap1: gap with base/N pairing\n"
        << "  --gap2 SCORE         gap2 penalty (default: 0.0)\n"
        << "                       gap2: gap-gap, N-N, N-base pairing\n"
        << "  -h, --help           Show this help message\n"
        << "\n"
        << "Examples:\n"
        << "  " << prog << " -i alignment.fa\n"
        << "  " << prog << " -i alignment.fa --match 2 --mismatch -3\n"
        << "\n"
        << "Output:\n"
        << "  SP score:    Total score for all columns\n"
        << "  Avg SP:      Average score per pair of sequences\n"
        << "  Scaled SP:   Average score per pair per column (length normalized)\n";
}

// 解析形如 "--match <double>" 的参数。
// idx 会被推进到 value 位置，便于外层 for 循环继续。
static bool parse_double_arg(int& idx, int argc, char** argv, double& out) {
    if (idx + 1 >= argc) return false;
    out = std::strtod(argv[idx + 1], nullptr);
    idx += 1;
    return true;
}

#ifndef NO_MAIN
int main(int argc, char** argv) {
    // 中文说明：
    // 默认参数沿用常见的 SP score 设置：match=+1 mismatch=-1 gap1=-2 gap2=0
    // gap2=0 的含义通常是：gap-gap 或 N-相关不额外计分/罚分（取决于定义）。
    std::string input_path;
    double matchS = 1.0, mismatchS = -1.0, gap1S = -2.0, gap2S = 0.0;
    bool use_streaming = true;  // 是否使用流式模式

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
            std::cerr << "unknown args: " << a << "\n\n";
            usage(argv[0]);
            return 2;
        }
    }

    if (input_path.empty()) {
        std::cerr << "Error: input file (-i) must be specified\n\n";
        usage(argv[0]);
        return 2;
    }


    try {
        SPScoreResult result;

        if (use_streaming) {
            result = calculate_sp_score_streaming(input_path, matchS, mismatchS, gap1S, gap2S);
        } else {
            // 传统模式：一次性加载所有序列
            std::vector<std::string> seqs = load_fasta(input_path);
            result = calculate_sp_score(seqs, matchS, mismatchS, gap1S, gap2S);
        }

        // Output results (using tabs for easier script parsing)
        std::cout.setf(std::ios::fixed);
        std::cout << std::setprecision(6);
        std::cout << "SP score\t" << result.total_sp << "\n";
        std::cout << "Avg SP\t" << result.avg_sp << "\n";
        std::cout << "Scaled SP\t" << result.scaled_sp << "\n";

    } catch (const std::exception& e) {
        // 将错误信息输出到 stderr 并返回非 0
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
#endif // NO_MAIN
