# 分岐解析と精度保証プロジェクト

このプロジェクトは、精度保証付き数値計算を用いて偏微分方程式の分岐解（非自明解）の探索および数学的な精度保証（真の解の存在証明）を行うためのものです。

## ディレクトリ構成

今後実験を繰り返す際にデータやコードが混ざって混乱しないよう、コード（`src/`）とデータ（`experiments/`）を完全に分離した構成を採用しています。

```text
antitest/
├── src/                      # すべてのソースコードと実行スクリプトの置き場
│   ├── hanpuku.cpp           # 初期探索（ニュートン法）用プログラム
│   ├── resize.cpp            # 近似解の次数拡張用プログラム
│   ├── verify.cpp            # 精度保証（Kantorovichの定理）用プログラム
│   ├── plot_results.py       # 解の形状グラフ生成スクリプト
│   ├── plot_bifurcation.py   # 分岐図生成スクリプト
│   └── run_experiment.sh     # 一括自動実行シェルスクリプト
│
├── experiments/              # 実験ごとのデータを保存するディレクトリ
│   ├── exp01_initial/        # 初期テスト（旧 hanpuku）
│   ├── exp02_fine_mesh/      # 高解像度探索（旧 hanpuku2）
│   ├── exp03_large_eps/      # 大きな ε での精度保証成功（旧 hanpuku3）
│   └── exp04_mpfr10000/      # MPFR 10000 bits での検証（旧 hanpuku4）
│
└── experiment_summary.md     # これまでの実験の目的・設定・結果のまとめ
```

## 新しい実験の実行方法

新しくパラメータ（例: εの範囲や近似次数など）を変更して実験を行いたい場合は、以下の手順で実行してください。

1. `src/` ディレクトリ内のプログラム（`hanpuku.cpp` 等）の必要なパラメータを書き換えます。
2. `src/` ディレクトリに移動し、`run_experiment.sh` を用いて、新しい実験の保存先ディレクトリを指定して実行します。

```bash
cd src/
./run_experiment.sh ../experiments/exp05_new_test
```

スクリプトが自動的に以下を実行します：
- 指定されたディレクトリ（例: `exp05_new_test/`）の作成
- 各C++プログラムのコンパイル
- `hanpuku` の実行（`data/` に初期解を保存）
- `resize` の実行（`80data/` に次数拡張した解を保存）
- `verify` の実行（`results/` に精度保証結果を保存）
- Pythonスクリプトによる各解の形状グラフおよび分岐図の生成

すべての出力データと画像は、指定した実験ディレクトリ（例: `../experiments/exp05_new_test/`）の中に整理されて保存されます。
