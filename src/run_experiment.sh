#!/bin/bash
# 使い方: ./run_experiment.sh ../experiments/exp05_new

if [ -z "$1" ]; then
    echo "エラー: 実験ディレクトリを指定してください。"
    echo "例: ./run_experiment.sh ../experiments/exp05_new"
    exit 1
fi

EXP_DIR=$1
mkdir -p "$EXP_DIR"

echo "=== コンパイル開始 ==="
g++ -I.. -DNDEBUG -DKV_FASTROUND -O3 hanpuku.cpp -llapack -lblas -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lmpfr -fopenmp -o hanpuku
g++ -I.. -std=c++17 -DNDEBUG -DKV_FASTROUND -O3 resize.cpp -llapack -lblas -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lmpfr -fopenmp -o resize
g++ -I.. -std=c++17 -DNDEBUG -DKV_FASTROUND -O3 verify.cpp -llapack -lblas -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lmpfr -fopenmp -o verify

echo "=== 1. hanpuku 実行 (初期探索) ==="
./hanpuku "$EXP_DIR" > "$EXP_DIR/hanpuku_output.txt"

echo "=== 2. resize 実行 (次数拡張) ==="
./resize "$EXP_DIR" > "$EXP_DIR/resize_output.txt"

echo "=== 3. verify 実行 (精度保証) ==="
./verify "$EXP_DIR" > "$EXP_DIR/verify_output.txt"

echo "=== 4. グラフ生成 ==="
python3 plot_results.py "$EXP_DIR"
python3 plot_bifurcation.py "$EXP_DIR"

echo "=== すべて完了しました！ ==="
echo "結果は $EXP_DIR に保存されています。"
