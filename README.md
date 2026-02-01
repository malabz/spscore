# SPScore

A fast, memory-efficient tool to calculate Sum-of-Pairs (SP) scores for Multiple Sequence Alignments (MSA).

## 🚀 Quick Start

### Build
```bash
# Main program
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j

# Test suite
cd test
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
```

### Run Tests
```bash
# Using scripts (Recommended)
./run_tests.sh      
# Manual
./test/build/spscore_test
```

## 💻 Usage

```bash
./spscore -i <alignment.fasta[.gz]> [options]
```

### Options:
- `-i, --input FILE`      Input MSA file (FASTA format, supports .gz)
- `-t, --threads NUM`     Number of threads (default: all cores)
- `--match SCORE`         Match score (default: 1.0)
- `--mismatch SCORE`      Mismatch penalty (default: -1.0)
- `--gap1 SCORE`          Gap-Base penalty (default: -2.0)
- `--gap2 SCORE`          Gap-Gap penalty (default: 0.0)
- `-h, --help`            Display help message

### Example:
```bash
./spscore -i alignment.fasta.gz --threads 8
```

## 🛠 Features
- **O(L) Memory**: Streaming mode processes sequences without loading the entire alignment.
- **Parallel Execution**: Uses OpenMP for multi-threaded column-wise calculations.
- **Robust Parsing**: Built-in FASTA/GZIP support via `kseq.h`.
- **Flexible Scoring**: Customizable penalties for different pair types.


