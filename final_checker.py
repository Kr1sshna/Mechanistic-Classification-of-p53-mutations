import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

def main():
    num = 0
    my_data = pd.read_csv("scan_results.csv")
    foldx_data = pd.read_csv("foldxresults.csv")
    
    if os.path.exists("scan_results.csv"):
        print("scan_results.csv exists")
        num+=1
    if os.path.exists("foldxresults.csv"):
        print("foldxresults.csv exists\n")
        num+=1

    if num != 2: #checks if scan_results.csv and foldxresults.csv exists
        print("File doesnt exist, error")
        return

    merged = pd.merge(my_data, foldx_data, on="Mutation")
    
    print("--- MERGED DATA ---")
    print(merged)
    merged.to_csv("merged_data.csv", index=False)
    

    plt.figure(figsize=(12, 12))
    
    # Scatter Plot
    ax=sns.scatterplot(
        data=merged, 
        x="FoldX_DDG", 
        y="Delta_Delta_E",
        hue = "Type",
        palette="Set2",
        s=80
    )
    
    # Labels
    texts = []
    for i in range(merged.shape[0]):
        texts.append(
            plt.text(
                merged.FoldX_DDG[i]+0.2,
                merged.Delta_Delta_E[i], 
                merged.Mutation[i], 
                fontsize=12
            )
        )

    adjust_text( #uses adjust text library to avoid overlapping texts in the plot
        texts,
        ax=ax,
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5)
    )

    plt.title("Correlation Graph", fontsize=16)   #title,x-axis label and y-axis label of the plot
    plt.xlabel("FoldX Predicted $\Delta\Delta G$ (kcal/mol)", fontsize=14)
    plt.ylabel("Calculated Electrostatic Destabilization $\Delta\Delta E$ (kcal/mol)", fontsize=14)
    
    # Draw Lines
    plt.axhline(0, color='gray', linestyle='--')
    plt.axvline(0, color='gray', linestyle='--')
    plt.grid(True, alpha=0.3)
    
    plt.savefig("final_correlation.png") #saves the plot
    plt.show()

if __name__ == "__main__":
    main()
