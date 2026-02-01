#ifndef SPSCORE_H
#define SPSCORE_H

#include <cstdint>
#include <vector>
#include <string>

// ============================================================================
// spscore.h（中文注释版）
// ----------------------------------------------------------------------------
// 这个头文件声明了计算 SP score（Sum-of-Pairs score）的核心 API。
//
// SP score 简述：
// - 给定一组已经做过多序列比对（MSA）的序列（要求等长），在每一列上对所有序列做
//   两两配对（pair），按 match/mismatch/gap 规则累加得到总分。
// - 本项目用“计数法”快速计算每一列的 pair 数量：只要知道该列 A/C/G/T/N/'-' 的
//   个数，就能 O(1) 算出该列所有 pair 的分数贡献。
//
// 字符约定：
// - 统一使用 6 类字符：'A','C','G','T','N','-'
// - 'U'（RNA）会被视为 'T'
// - 非法字符会被视为 'N'
// ============================================================================

// Normalize to {'A','C','G','T','N','-'} with U->T and invalid->N.
//
// 中文说明：
// 将输入字符归一化到 {'A','C','G','T','N','-'} 六类。
// - 支持大小写；
// - 'U'/'u' 视为 'T'；
// - 其它未知字符一律转为 'N'。
char normalize_base(char c);

// Calculate binomial coefficient: n choose 2
//
// 中文说明：
// 计算 C(x,2)=x*(x-1)/2，用于求“同类字符两两配对”的数量。
// 例如某列中 A 出现 x 次，则 A-A 的配对数为 C(x,2)。
uint64_t choose2_u64(uint64_t x);

// Calculate SP score for a column given base counts
// cnt: [A,C,G,T,N,'-'] in indices 0..5
//
// 中文说明：
// 计算“单列”的 SP score 总贡献（该列所有序列两两配对的分数总和）。
//
// 输入：
// - cnt[6]：该列各字符计数，按以下顺序：
//   cnt[0]=A, cnt[1]=C, cnt[2]=G, cnt[3]=T, cnt[4]=N, cnt[5]='-'
// - matchS / mismatchS / gap1S / gap2S：对应四类配对的权重（得分/惩罚）
//
// 输出：
// - double：该列所有 pair 的总得分。
//
// 重要说明（语义）：
// - mismatch 只在 A/C/G/T 之间统计不同碱基的配对。
// - '-' 与 (A/C/G/T/N) 的配对计为 gap1。
// - gap2 计入：gap-gap、N-N、以及 N 与 (A/C/G/T) 的配对。
//   这套 gap1/gap2 的分类方式会影响单元测试期望值，请修改前确保同步更新测试。
double score_of_counts(const uint32_t cnt[6],
                       double matchS, double mismatchS,
                       double gap1S,  double gap2S);

// Calculate total SP score for multiple sequence alignment
// Returns: {total_sp, avg_sp, scaled_sp}
//
// 中文说明：
// 用于保存 calculate_sp_score 的输出结果：
// - total_sp：全比对（所有列）总 SP score
// - avg_sp：  total_sp / C(M,2)，即“平均每一对序列”的总分
// - scaled_sp：avg_sp / L，即“每一对序列、每一列”的平均分（长度归一化）
struct SPScoreResult {
    double total_sp;
    double avg_sp;
    double scaled_sp;
};

// 中文说明：
// 计算一组对齐序列（MSA）的 SP score。
//
// 输入：
// - seqs：序列数组，必须满足：
//   1) 至少 2 条序列；
//   2) 每条序列长度完全相同（已经对齐）；
//   3) 建议字符已归一化为 A/C/G/T/N/-（若未归一化，调用者应先处理；或使用 load_fasta）。
// - matchS/mismatchS/gap1S/gap2S：计分参数（提供默认值）。
//
// 输出：
// - SPScoreResult：包含 total/avg/scaled。
//
// 异常：
// - 若 seqs 数量不足、长度不一致、或为空，会抛出 std::invalid_argument。
SPScoreResult calculate_sp_score(const std::vector<std::string>& seqs,
                                  double matchS = 1.0,
                                  double mismatchS = -1.0,
                                  double gap1S = -2.0,
                                  double gap2S = 0.0);

// Load FASTA file (supports gzip)
//
// 中文说明：
// 从 FASTA（支持 .gz）读取序列，并对每个字符执行 normalize_base 归一化。
// 返回的序列可直接用于 calculate_sp_score。
std::vector<std::string> load_fasta(const std::string& path);

// 流式计算 SP score（内存友好版本）
//
// 中文说明：
// 从 FASTA 文件流式读取并计算 SP score，不需要将所有序列加载到内存。
// 采用两遍扫描策略：
// - 第一遍：获取序列数量和长度，验证对齐
// - 第二遍：逐序列读取，逐列累积计数并计算得分
//
// 输入：
// - path：FASTA 文件路径（支持 .gz 压缩）
// - matchS/mismatchS/gap1S/gap2S：计分参数
//
// 输出：
// - SPScoreResult：包含 total/avg/scaled
//
// 优点：
// - 内存占用 O(L)，L 为序列长度，不依赖序列数量
// - 适合处理大量序列（数千到数万条）
//
// 异常：
// - 若文件无法打开、序列数量不足、长度不一致等，会抛出异常
SPScoreResult calculate_sp_score_streaming(const std::string& path,
                                            double matchS = 1.0,
                                            double mismatchS = -1.0,
                                            double gap1S = -2.0,
                                            double gap2S = 0.0);

#endif // SPSCORE_H

