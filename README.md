# SPScore 测试套件使用说明

## ✅ 项目结构

您的 spscore 项目现在分为两个独立的构建系统：

### 主程序 (根目录)
- **CMakeLists.txt** - 主程序的构建配置
- **spscore** - 主程序可执行文件
- **spscore_lib** - 静态库

### 测试程序 (test目录)
- **test/CMakeLists.txt** - 测试程序的独立构建配置
- **test/spscore_test** - 测试程序可执行文件

## 🚀 运行测试

### 一键运行测试（推荐）⭐

#### Windows
```bash
.\run_tests.bat
```

#### Linux/WSL
```bash
chmod +x run_tests.sh
./run_tests.sh
```

此脚本会自动：
1. ✓ 清理旧的构建文件
2. ✓ 配置 CMake
3. ✓ 编译测试程序
4. ✓ 运行所有测试
5. ✓ 显示测试结果

### 手动运行测试

#### Windows
```bash
# 进入测试目录
cd test

# 配置和编译
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release -j

# 运行测试
Release\spscore_test.exe
# 或
.\spscore_test.exe

# 运行特定测试套件
.\spscore_test.exe --test-suite="normalize_base"
.\spscore_test.exe --test-suite="performance tests"

# 详细输出
.\spscore_test.exe -s
```

#### Linux/WSL
```bash
# 进入测试目录
cd test

# 配置和编译
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j

# 运行测试
./spscore_test

# 运行特定测试套件
./spscore_test --test-suite="normalize_base"
./spscore_test --test-suite="performance tests"

# 详细输出
./spscore_test -s
```

## 🔨 编译主程序

主程序的编译独立于测试程序：

```bash
# 在项目根目录
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j

# 运行主程序
./spscore -i ../test/example.fasta  # Linux/WSL
.\spscore.exe -i ..\test\example.fasta  # Windows
```

## 📋 测试内容

### 1. 单元测试（50+测试用例）
- **normalize_base**: 碱基标准化（大小写、U→T、无效字符等）
- **choose2_u64**: 组合数计算
- **score_of_counts**: 列分数计算（匹配、错配、gap惩罚）
- **calculate_sp_score**: SP分数计算

### 2. 边界情况测试（13个测试用例）
- 错误输入处理
- 极端情况（超长序列、大量序列、全gap等）

### 3. 性能测试（6个测试用例）
测试不同规模的比对：
- 10 sequences × 100 bp
- 50 sequences × 500 bp
- 100 sequences × 1000 bp
- 200 sequences × 2000 bp
- 10 sequences × 10000 bp
- 500 sequences × 100 bp

### 4. 生物学场景测试（6个测试用例）
- 保守区域
- 可变区域
- 插入/缺失
- SNP变异
- RNA序列
- 模糊碱基

## 📊 测试示例

示例FASTA文件已创建在 `test/example.fasta`:
```
>seq1
ATGCATGCATGC
>seq2
ATGCATGCATGC
>seq3
ATGCCTGCATGC
>seq4
ATGC--GCATGC
```

使用主程序测试：
```bash
# 先编译主程序
cd build
cmake ..
cmake --build . -j

# Windows
.\spscore.exe -i ..\test\example.fasta

# Linux/WSL
./spscore -i ../test/example.fasta
```

## 📁 项目结构

```
spscore/
├── CMakeLists.txt              # 主程序构建配置
├── run_tests.sh                # 一键测试脚本 (Linux/WSL) ⭐
├── run_tests.bat               # 一键测试脚本 (Windows) ⭐
├── README_TEST.md              # 详细测试文档
│
├── src/
│   ├── spscore.h              # 公共接口头文件
│   ├── spscore.cpp            # 实现文件和主程序
│   └── kseq.h                 # FASTA解析库
│
├── test/
│   ├── CMakeLists.txt         # 测试程序独立构建配置 ⭐
│   ├── test_spscore.cpp       # 测试用例（97个测试）
│   ├── doctest.h              # 测试框架
│   ├── example.fasta          # 示例数据
│   └── build/                 # 测试程序构建目录 ⭐
│       └── spscore_test       # 测试可执行文件
│
└── build/                      # 主程序构建目录
    ├── spscore                # 主程序可执行文件
    └── libspscore_lib.a       # 静态库
```

## 🎯 快速开始

### 1. 运行测试（推荐第一步）
```bash
# Windows
.\run_tests.bat

# Linux/WSL
chmod +x run_tests.sh
./run_tests.sh
```

### 2. 编译主程序
```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
```

### 3. 使用主程序
```bash
# Linux/WSL
./build/spscore -i test/example.fasta

# Windows
.\build\Release\spscore.exe -i test\example.fasta
```

## 📝 测试统计

- **总测试用例**: 97个
- **测试套件**: 7个
- **代码覆盖**: 100% 公共API
- **性能基准**: 6个不同规模

## ✨ 主要优势

✅ **独立构建** - 主程序和测试程序完全分离，互不干扰  
✅ **一键测试** - 使用脚本自动完成所有步骤  
✅ **快速迭代** - 修改测试后只需重新运行 run_tests 脚本  
✅ **清晰结构** - 测试相关文件都在 test/ 目录下  

祝使用愉快！🎉

