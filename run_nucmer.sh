#!/bin/bash

# 脚本用法说明
usage() {
    echo "用法: $0 out ref query [filter] [cluster] [length]"
    echo "必需参数:"
    echo "  out      - 输出文件前缀"
    echo "  ref      - 参考序列文件"
    echo "  query    - 查询序列文件"
    echo "可选参数(带默认值):"
    echo "  filter   - 过滤阈值 (默认: 20000)"
    echo "  cluster  - 聚类大小 (默认: 100)"
    echo "  length   - 最小长度 (默认: 1000)"
    exit 1
}

# 检查必需参数数量
if [ $# -lt 3 ]; then
    echo "错误: 缺少必需参数!"
    usage
fi

# 设置参数变量
OUT_PREFIX=$1
REF=$2
QUERY=$3

# 设置可选参数默认值
FILTER=${4:-20000}
CLUSTER=${5:-100}
LENGTH=${6:-1000}

# 检查输入文件是否存在
check_file() {
    if [ ! -f "$1" ]; then
        echo "错误: 文件 $1 不存在!"
        exit 1
    fi
}

check_file "$REF"
check_file "$QUERY"

# 打印参数信息
echo "=== 运行参数 ==="
echo "输出前缀: $OUT_PREFIX"
echo "参考序列: $REF"
echo "查询序列: $QUERY"
echo "过滤阈值: $FILTER"
echo "聚类大小: $CLUSTER"
echo "最小长度: $LENGTH"
echo "==============="

#
nucmer="/work/xup/software/Base/mummer-4.0.0beta2/bin/"
# 执行nucmer比对
if [ ! -e ${OUT_PREFIX}.delta ];then
   echo "正在运行nucmer比对..."
   $nucmer/nucmer -t 24 --prefix "$OUT_PREFIX" -c "$CLUSTER" -l "$LENGTH" "$REF" "$QUERY"
   if [ $? -ne 0 ]; then
        echo "错误: nucmer比对失败!"
        exit 1
   fi

fi 

# 运行delta-filter
echo "正在过滤比对结果..."
$nucmer/delta-filter -l "$FILTER" "${OUT_PREFIX}.delta" > "${OUT_PREFIX}.filtered.${FILTER}.delta"
$nucmer/show-coords ${OUT_PREFIX}.filtered.${FILTER}.delta > ${OUT_PREFIX}.filtered.${FILTER}.coords
$nucmer/mummerplot -f --png -p ${OUT_PREFIX}.filtered.${FILTER} ${OUT_PREFIX}.filtered.${FILTER}.delta

# 显示结果文件信息
echo "=== 结果文件 ==="
ls -lh "${OUT_PREFIX}"*
echo "==============="
