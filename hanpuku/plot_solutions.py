import os
import glob
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import io
from collections import defaultdict

DATA_DIR = "data"

# _uh のバイナリファイルをすべて検索
files = glob.glob(os.path.join(DATA_DIR, "*_uh.matrix_kvdd"))

# パラメータ (eps, alpha) ごとに、特徴の違う解（分岐）をグループ化する
grouped_solutions = defaultdict(dict)

for f in files:
    # ファイル名から情報を取り出す (例: eps_0.0350_alpha_0.25_trial_1_H10_3.50103_uh.matrix_kvdd)
    basename = os.path.basename(f)
    name_no_ext = basename.replace(".matrix_kvdd", "")
    
    parts = name_no_ext.split("_")
    try:
        eps = float(parts[1])
        alpha = float(parts[3])
        h10 = float(parts[7])
    except Exception as e:
        continue
        
    # 自明解(H10がほぼ0)は除外する
    if h10 < 1e-3:
        continue
        
    # H10を丸めて、同じ枝(同じ解)とみなせるものをまとめる
    rounded_h10 = round(h10, 2)
    group_key = (eps, alpha)
    
    if rounded_h10 not in grouped_solutions[group_key]:
        grouped_solutions[group_key][rounded_h10] = {
            "path": os.path.join(DATA_DIR, name_no_ext), # vcp::loadは拡張子を自動付与するため拡張子なし
            "h10": h10,
            "eps": eps,
            "alpha": alpha
        }

if not grouped_solutions:
    print("非自明解が見つかりませんでした。")
    exit(0)

# (eps, alpha) のペアごとに別々のグラフを作成する
for (eps, alpha), distinct_sols in grouped_solutions.items():
    plt.figure(figsize=(10, 6))
    
    for rh10, sol in distinct_sols.items():
        print(f"Exporting and plotting data for eps={eps}, alpha={alpha}, H10={sol['h10']}...")
        
        # exporterを実行
        cmd = ["./exporter", sol["path"]]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error running exporter on {sol['path']}")
            continue
            
        # CSV形式の標準出力をPandasで読み込む
        csv_data = result.stdout
        df = pd.read_csv(io.StringIO(csv_data), header=None, names=["x", "y"])
        
        label = f"H10={sol['h10']:.3f}"
        plt.plot(df["x"], df["y"], label=label, linewidth=2)

    plt.xlabel("x", fontsize=14)
    plt.ylabel("u_h(x)", fontsize=14)
    plt.title(f"Non-Trivial Solutions (eps={eps}, alpha={alpha})", fontsize=16)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)
    plt.tight_layout()
    
    # グラフの保存
    filename = f"solutions_graph_eps{eps}_alpha{alpha}.png"
    plt.savefig(filename, dpi=300)
    plt.close()  # 次の画像に前のプロットが残らないように閉じる
    print(f"Graph saved as {filename}")
