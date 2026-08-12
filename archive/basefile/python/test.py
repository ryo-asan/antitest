import matplotlib.pyplot as plt
import re
import os

# --- 設定: 入力ファイルと出力フォルダを指定してください ---
input_file = os.path.join('0113_2', 'test.txt') # 'data'フォルダ内の'test.txt'
output_folder = 'bifurcation_results_0113_2'     # 保存先
# ------------------------------------------------------

# フォルダの作成
os.makedirs(output_folder, exist_ok=True)

if not os.path.exists(input_file):
    print(f"エラー: {input_file} が見つかりません。")
else:
    with open(input_file, 'r', encoding='utf-8') as f:
        content = f.read()

    # test.txtの特殊な形式に対応する正規表現
    # epsilon:[値]uh_max:[値]
    # || uh ||_H10[値]
    # || uh ||_L2[値]
    # が連続している部分を1セットとして抽出します
    pattern = r'epsilon:([\d.e+-]+)uh_max:([\d.e+-]+)\s*\|\| uh \|\|_H10([\d.e+-]+)\s*\|\| uh \|\|_L2([\d.e+-]+)'
    matches = re.findall(pattern, content)

    eps_vals = []
    linf_vals = []
    h10_vals = []
    l2_vals = []

    for m in matches:
        eps_vals.append(float(m[0]))
        linf_vals.append(float(m[1])) # uh_maxをL-infとして扱う
        h10_vals.append(float(m[2]))
        l2_vals.append(float(m[3]))

    print(f"抽出成功: {len(eps_vals)} 点のデータを取得しました。")

    if len(eps_vals) == 0:
        print("データが抽出されませんでした。正規表現またはファイル内容を確認してください。")
    else:
        # プロット用関数（点のみ）
        def save_plot(x, y, ylabel, title, filename, color):
            plt.figure(figsize=(8, 6))
            plt.plot(x, y, 'o', markersize=3, color=color, alpha=0.5)
            plt.xlabel(r'$\epsilon$')
            plt.ylabel(ylabel)
            plt.title(title)
            plt.grid(True, linestyle=':', alpha=0.6)
            save_path = os.path.join(output_folder, filename)
            plt.savefig(save_path, dpi=300)
            plt.close()
            print(f"保存完了: {save_path}")

        # 3つのグラフを出力
        save_plot(eps_vals, linf_vals, r'$L^\infty$ norm', 
                  r'Bifurcation Diagram: $\epsilon$ vs $L^\infty$', 'bifurcation_linf.png', 'blue')
        
        save_plot(eps_vals, l2_vals, r'$L^2$ norm', 
                  r'Bifurcation Diagram: $\epsilon$ vs $L^2$', 'bifurcation_l2.png', 'orange')
        
        save_plot(eps_vals, h10_vals, r'$H^1_0$ norm', 
                  r'Bifurcation Diagram: $\epsilon$ vs $H^1_0$', 'bifurcation_h10.png', 'green')