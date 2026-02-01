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
echo "========================================="
if ./spscore_test; then
    TEST_EXIT_CODE=0
else
    TEST_EXIT_CODE=$?
fi
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
    echo ""
    exit 0
else
    echo -e "${RED}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${RED}✗ 测试失败 (退出代码: $TEST_EXIT_CODE)${NC}"
    echo -e "${RED}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    exit $TEST_EXIT_CODE
fi

