# SPScore

A fast, memory-efficient tool to calculate Sum-of-Pairs (SP) scores for Multiple Sequence Alignments (MSA).

##  Install (Conda)

```bash
conda install -c malab spscore
```

## Install (cmake)
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


## 🚀 Quick Start


## 💻 Usage

```bash
./spscore -i <alignment.fasta[.gz]> [options]
```

### Options:
- `-i, --input FILE`      Input MSA file (FASTA format, supports .gz)
- `--match SCORE`         Match score (default: 1.0)
- `--mismatch SCORE`      Mismatch penalty (default: -1.0)
- `--gap1 SCORE`          Gap-Base penalty (default: -2.0)
- `--gap2 SCORE`          Gap-Gap penalty (default: 0.0)
- `-h, --help`            Display help message

### Example:
```bash
./spscore -i ./test/data/mt1x.aligned.fasta.gz
```
output:
```
SP score        3715981575.000000
Avg SP  16482.070005
Scaled SP       0.990926
```

### Run Tests
```bash
# Using scripts (Recommended)
./run_tests.sh      
```

## 🛠 Features
- **O(L) Memory**: Streaming mode processes sequences without loading the entire alignment.
- **Robust Parsing**: Built-in FASTA/GZIP support via `kseq.h`.
- **Flexible Scoring**: Customizable penalties for different pair types.
