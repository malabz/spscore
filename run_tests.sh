#!/bin/bash

# SPScore 测试一键运行脚本
# 自动编译和运行测试套件

set -e  # 遇到错误立即退出

echo "========================================="
echo "SPScore 测试套件 - 一键运行"
echo "========================================="
echo ""

# 定义颜色
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 获取脚本所在目录
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# 测试构建目录
TEST_BUILD_DIR="test/build"

# -----------------------------------------------------------------------------
# 测试输出日志
# - 需求："输出测试时候打印的相关信息"。
# - 做法：把测试可执行文件的 stdout/stderr 同时：
#   1) 实时打印到终端；
#   2) 保存到日志文件，便于事后查看/粘贴。
# -----------------------------------------------------------------------------
LOG_DIR="$TEST_BUILD_DIR/logs"
TS="$(date +%Y%m%d_%H%M%S)"
TEST_LOG="$LOG_DIR/spscore_test_${TS}.log"

# 先创建目录（后面步骤 1 会清理 build 目录，因此这里先仅定义变量，创建在步骤 2 之后）

echo "步骤 1/4: 清理旧的构建文件..."
if [ -d "$TEST_BUILD_DIR" ]; then
    rm -rf "$TEST_BUILD_DIR"
    echo -e "${GREEN}✓${NC} 清理完成"
else
    echo -e "${YELLOW}→${NC} 无需清理"
fi
echo ""

echo "步骤 2/4: 配置 CMake..."
mkdir -p "$TEST_BUILD_DIR"
cd "$TEST_BUILD_DIR"
if cmake .. -DCMAKE_BUILD_TYPE=Release; then
    echo -e "${GREEN}✓${NC} CMake 配置成功"
else
    echo -e "${RED}✗${NC} CMake 配置失败"
    exit 1
fi
echo ""

echo "步骤 3/4: 编译测试程序..."
if cmake --build . -j$(nproc 2>/dev/null || echo 4); then
    echo -e "${GREEN}✓${NC} 编译成功"
else
    echo -e "${RED}✗${NC} 编译失败"
    exit 1
fi
echo ""

echo "步骤 4/4: 运行测试..."
mkdir -p "$LOG_DIR"
echo "- 测试输出将同时打印到屏幕，并保存到: $TEST_LOG"
echo "========================================="
# 用 tee 同时输出到终端与日志；PIPESTATUS[0] 取到被 tee 包裹的 ./spscore_test 的退出码
./spscore_test 2>&1 | tee "$TEST_LOG"
TEST_EXIT_CODE=${PIPESTATUS[0]}
echo "========================================="
echo ""

# 返回到项目根目录
cd "$SCRIPT_DIR"

# 显示测试结果
if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${GREEN}✓ 所有测试通过！${NC}"
    echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    echo "测试可执行文件位置: $TEST_BUILD_DIR/spscore_test"
    echo "测试输出日志位置:   $TEST_LOG"
    echo ""
    exit 0
else
    echo -e "${RED}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${RED}✗ 测试失败 (退出代码: $TEST_EXIT_CODE)${NC}"
    echo -e "${RED}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    echo "测试输出日志位置:   $TEST_LOG"
    echo "（可用 'tail -n 200 $TEST_LOG' 快速查看末尾日志）"
    echo ""
    exit $TEST_EXIT_CODE
fi

