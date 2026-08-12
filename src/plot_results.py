import os
import pandas as pd
import matplotlib.pyplot as plt
import sys

base_dir = sys.argv[1] if len(sys.argv) > 1 else "."
RESULTS_DIR = os.path.join(base_dir, "results")

if not os.path.exists(RESULTS_DIR):
    print("Results directory not found.")
    exit(0)

# Iterate through all subdirectories in results
for item in os.listdir(RESULTS_DIR):
    item_path = os.path.join(RESULTS_DIR, item)
    if os.path.isdir(item_path):
        csv_file = os.path.join(item_path, "graph_data.csv")
        if os.path.exists(csv_file):
            df = pd.read_csv(csv_file, header=None, names=["x", "y"])
            
            plt.figure(figsize=(10, 6))
            plt.plot(df["x"], df["y"], linewidth=2, color="blue")
            plt.xlabel("x", fontsize=14)
            plt.ylabel("u_h(x)", fontsize=14)
            
            # Extract info from directory name for title
            # Format: eps_0.0300_alpha_0.25_trial_27_H10_1.20345
            title = "Approximate Solution"
            parts = item.split("_")
            if len(parts) >= 8:
                eps = parts[1]
                alpha = parts[3]
                h10 = parts[7]
                title = f"eps={eps}, alpha={alpha}, H10={h10}"
                
            plt.title(title, fontsize=16)
            plt.grid(True)
            plt.tight_layout()
            
            plt.savefig(os.path.join(item_path, "graph.png"), dpi=300)
            plt.close()
            print(f"Plotted graph for {item}")
