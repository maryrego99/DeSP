import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

TIMESTAMP = datetime.now().strftime("_%d-%m_%H.%M.%S")

def plot_recovery_vs_coverage(csv_file):
    df = pd.read_csv(csv_file) 

    df['mean_coverage'] = pd.to_numeric(df['mean_coverage'], errors='coerce')
    df['oligo_recovery(seq)'] = pd.to_numeric(df['oligo_recovery(seq)'], errors='coerce')
    df['oligo_recovery(decode)'] = pd.to_numeric(df['oligo_recovery(decode)'], errors='coerce')
    df['dropout(seq)'] = pd.to_numeric(df['dropout(seq)'], errors='coerce')
    df['dropout(decode)'] = pd.to_numeric(df['dropout(decode)'], errors='coerce')

    alpha = df['α'].iloc[0]
    subs_rate = df['subs_rate'].iloc[0]

    plt.figure(figsize=(12, 6))

    if any(df['decode_success'] == 'Yes'):
        first_success = df[df['decode_success'] == 'Yes'].iloc[0]
        success_coverage = first_success['mean_coverage']
        success_recovery = first_success['oligo_recovery(decode)']
        success = df[df['decode_success'] == 'Yes']
        fail = df[df['decode_success'] == 'No']
        # plt.plot(success['mean_coverage'], success['oligo_recovery(seq)'], 'o-', label='Decoded: Yes', color='green', zorder=3)
        plt.plot(success['mean_coverage'], success['oligo_recovery(decode)'], 'o-', label='Decoded: Yes', color='green', zorder=3)
    else:
        success_coverage = "N/A. Not decoded successfully"
        success_recovery = "N/A"
        fail = df[df['decode_success'] == 'No']

    
    plt.plot(df['mean_coverage'], df['oligo_recovery(seq)'], 'o-', label='Recovery (Synthesis)', color='blue')
    plt.plot(df['mean_coverage'], df['oligo_recovery(decode)'], 'o-', label='Recovery (Decode)', color='purple')
    plt.xlabel('Mean Coverage')
    plt.ylabel('% Oligo Recovery')
    plt.title('% Oligo Recovery: Synthesis vs Decode (RS)')

    # plt.xticks(np.arange(0, 10.5 + 0.5, 0.5))
    # plt.xlim(0, 10.5)
    plt.xticks(np.arange(0, 11, 1))
    plt.xlim(0, 11)
    plt.yticks(np.arange(0, 111, 10))
    plt.ylim(40, 110)

    if any(df['decode_success'] == 'Yes'):
        first_success = df[df['decode_success'] == 'Yes'].iloc[0]
        success_coverage = first_success['mean_coverage']
        success_recovery = first_success['oligo_recovery(decode)']
    else:
        success_coverage = "N/A"
        success_recovery = "N/A"

    param_text = (
        f'α = {alpha}\n'
        f'Substitution Error Rate = {subs_rate}\n'
        f'First successful decode:\n'
        # f'  Coverage ≈ {success_coverage:.2f}\n'
        # f'  Recovery ≈ {success_recovery:.2f}%'
        f'  Coverage ≈ {success_coverage}\n'
        f'  Recovery ≈ {success_recovery}%'
    )
    plt.gca().text(
        0.70, 0.25, param_text,
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray')
    )

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    # plt.show()
    plt.savefig("coverage-analysis/seq-depth/visualizations/recovery_vs_coverage/rs_2_a="+str(alpha)+".png")

def plot_oligo_copy_distributions(syn_file, pcr_file, seq_file, save_path=None):
    
    
    syn_df = pd.read_csv(syn_file)
    pcr_df = pd.read_csv(pcr_file)
    seq_df = pd.read_csv(seq_file)

    #merge files by oligo index
    merged_df = syn_df.merge(pcr_df, on='oligo_index', suffixes=('_syn', '_pcr'))
    merged_df = merged_df.merge(seq_df, on='oligo_index')
    merged_df.rename(columns={'total_copies': 'total_copies_seq'}, inplace=True)
    dropout_indices = merged_df[merged_df['total_copies_seq'] == 0]['oligo_index']

    print(merged_df.head(10))
    print(merged_df.tail(10))

    #histogram
    plt.figure(figsize=(14, 6))
    bins = 50

    plt.hist(merged_df['total_copies_syn'], bins=bins, alpha=0.4, label='Synthesis', color='skyblue')
    # plt.hist(merged_df['total_copies_pcr'], bins=bins, alpha=0.4, label='After PCR', color='orange')
    plt.hist(merged_df['total_copies_seq'], bins=bins, alpha=0.4, label='After Sequencing', color='green')

    plt.xlabel('Number of Copies per Oligo')
    plt.ylabel('Number of Oligos')
    plt.title('Oligo Copy Distribution Across Pipeline Stages')
    plt.legend()
    plt.grid(True)

    if save_path:
        plt.savefig(save_path+"oligo-hist"+TIMESTAMP+".png", bbox_inches='tight')
        print(f"Saved plot to {save_path}")
    else:
        plt.show()


    # PCR
    plt.figure(figsize=(16, 6))
    plt.plot(merged_df['oligo_index'], merged_df['total_copies_pcr'], label='After PCR', marker='s', linewidth=1)
    plt.scatter(dropout_indices, [0]*len(dropout_indices), color='red', s=10, label='Dropouts', zorder=10)

    plt.xlabel('Oligo Index')
    plt.ylabel('Number of Copies')
    plt.title('Copy Count per Oligo After PCR')
    plt.legend()
    plt.grid(alpha=0.5)

    if save_path:
        plt.savefig(save_path+"oligo_dist_PCR"+TIMESTAMP+".png")
        print(f"Saved plot to {save_path+"oligo_dist_PCR"+TIMESTAMP+".png"}")
    else:
        plt.show()

    # Syn and Seq
    plt.figure(figsize=(18, 6))
    plt.plot(merged_df['oligo_index'], merged_df['total_copies_syn'], 
             label='Synthesis', color='blue', marker='o', markersize=4, linewidth=1, alpha=0.8)
    
    plt.plot(merged_df['oligo_index'], merged_df['total_copies_seq'], 
             label='After Sequencing', color='darkorange', marker='^', markersize=4, linewidth=1, alpha=0.8)
    
    
    plt.scatter(dropout_indices, [0]*len(dropout_indices), color='red', s=10, label='Dropouts', zorder=10)

    plt.title("Per-Oligo Copy Count: Synthesis vs Sequencing", fontsize=14)
    plt.xlabel("Oligo Index", fontsize=12)
    plt.ylabel("Number of Copies", fontsize=12)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    # plt.grid(True)

    if save_path:
        plt.savefig(save_path+"oligo_dist_across_stages"+TIMESTAMP+".png")
        print(f"Saved plot to {save_path+"oligo_dist_across_stage"+TIMESTAMP+".png"}")
    else:
        plt.show()

    


if __name__ == "__main__":
    # plot_oligo_copy_distributions(
    #     syn_file="coverage-analysis/seq-depth/files/syn_dna_counts.csv",
    #     pcr_file="coverage-analysis/seq-depth/files/pcr_copy_counts.csv",
    #     seq_file="coverage-analysis/seq-depth/files/seq_copy_counts.csv",
    #     save_path="coverage-analysis/seq-depth/visualizations/"
    # )

    plot_recovery_vs_coverage("coverage-analysis/seq-depth/files/coverage_metrics_rs_a=0.1.csv")