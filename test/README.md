# SPScore 测试程序

这是 spscore 项目的独立测试套件。

## 快速开始

### 方法1：一键运行（推荐）

在**项目根目录**运行：

```bash
# Linux/WSL
./run_tests.sh

# Windows
.\run_tests.bat
```

### 方法2：手动构建和运行

在**当前目录** (test/) 运行：

```bash
# 创建构建目录
mkdir build
cd build

# 配置 CMake
cmake ..

# 编译
cmake --build . -j

# 运行测试
./spscore_test              # Linux/WSL
.\spscore_test.exe          # Windows
```

## 测试选项

```bash
# 运行所有测试
./spscore_test

# 运行特定测试套件
./spscore_test --test-suite="normalize_base"
./spscore_test --test-suite="choose2_u64"
./spscore_test --test-suite="score_of_counts"
./spscore_test --test-suite="calculate_sp_score"
./spscore_test --test-suite="edge cases and error handling"
./spscore_test --test-suite="performance tests"
./spscore_test --test-suite="biological scenarios"

# 详细输出
./spscore_test -s

# 列出所有测试
./spscore_test --list-test-cases

# 运行特定测试用例
./spscore_test --test-case="minimal case*"
```

## 测试内容

- **单元测试** (50+): 测试每个函数的正确性
- **边界测试** (13): 测试极端情况和错误处理
- **性能测试** (6): 测试不同规模的计算性能
- **场景测试** (6): 测试真实生物学应用

总计：**97个测试用例**

## 文件说明

- `CMakeLists.txt` - 测试程序的独立 CMake 配置
- `test_spscore.cpp` - 测试用例源代码
- `doctest.h` - doctest 测试框架头文件
- `example.fasta` - 示例测试数据
- `build/` - 构建输出目录（自动创建）

## 依赖项

- C++17 或更高版本
- CMake 3.28 或更高版本
- zlib
- OpenMP (可选，用于并行加速)

## 故障排除

### 编译错误
```bash
# 清理并重新构建
rm -rf build
mkdir build && cd build
cmake ..
cmake --build . -j
```

### 找不到源文件
确保您在 `test/` 目录下运行 cmake，或使用根目录的 `run_tests.sh` 脚本。

### OpenMP警告
如果系统不支持 OpenMP，测试仍然可以运行，只是速度会慢一些。

## 性能优化

编译 Release 版本以获得最佳性能：

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
```

## 添加新测试

编辑 `test_spscore.cpp`，使用 doctest 宏添加新的测试用例：

```cpp
TEST_CASE("my new test") {
    // 测试代码
    CHECK(1 + 1 == 2);
}
```

重新编译并运行测试即可。

