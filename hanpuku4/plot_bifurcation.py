import os
import glob
import matplotlib.pyplot as plt

RESULTS_DIR = "results"
data = []

# すべての verification.txt ファイルからデータを抽出
for txt_file in glob.glob(os.path.join(RESULTS_DIR, "*", "verification.txt")):
    eps = None
    h10 = None
    l2 = None
    linf = None
    
    with open(txt_file, 'r') as f:
        for line in f:
            if line.startswith("Epsilon ="):
                # 例: Epsilon = 0.045, Alpha = 0.25
                parts = line.split(",")
                eps_str = parts[0].split("=")[1].strip()
                eps = float(eps_str)
            elif "|| uh ||_H10" in line:
                val = line.split(":")[1].strip()
                h10 = float(val)
            elif "|| uh ||_L2" in line:
                val = line.split(":")[1].strip()
                l2 = float(val)
            elif "|| uh ||_Linf" in line:
                val = line.split(":")[1].strip()
                linf = float(val)
                
    if eps is not None and h10 is not None and l2 is not None and linf is not None:
        data.append({"eps": eps, "H10": h10, "L2": l2, "Linf": linf})

if not data:
    print("データが見つかりませんでした。")
    exit()

# epsilonでソート
data.sort(key=lambda x: x["eps"])

eps_vals = [d["eps"] for d in data]
h10_vals = [d["H10"] for d in data]
l2_vals = [d["L2"] for d in data]
linf_vals = [d["Linf"] for d in data]

# H10ノルムのグラフ
plt.figure(figsize=(8, 6))
plt.scatter(eps_vals, h10_vals, color='red', s=80, edgecolor='black', zorder=5)
plt.xlabel("$\epsilon$", fontsize=16)
plt.ylabel("$H^1_0$ norm", fontsize=16)
plt.title("Bifurcation Diagram: $H^1_0$ norm", fontsize=18)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig("bifurcation_H10.png", dpi=300)
plt.close()

# L2ノルムのグラフ
plt.figure(figsize=(8, 6))
plt.scatter(eps_vals, l2_vals, color='blue', s=80, edgecolor='black', zorder=5)
plt.xlabel("$\epsilon$", fontsize=16)
plt.ylabel("$L^2$ norm", fontsize=16)
plt.title("Bifurcation Diagram: $L^2$ norm", fontsize=18)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig("bifurcation_L2.png", dpi=300)
plt.close()

# Linfノルムのグラフ
plt.figure(figsize=(8, 6))
plt.scatter(eps_vals, linf_vals, color='green', s=80, edgecolor='black', zorder=5)
plt.xlabel("$\epsilon$", fontsize=16)
plt.ylabel("$L^\infty$ norm", fontsize=16)
plt.title("Bifurcation Diagram: $L^\infty$ norm", fontsize=18)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig("bifurcation_Linf.png", dpi=300)
plt.close()

print("分岐図(Bifurcation diagrams)の保存が完了しました。")
