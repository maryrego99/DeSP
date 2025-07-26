import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from matplotlib.lines import Line2D

def clean_data(csv_file):
    df = pd.read_csv(csv_file) 
    df['mean_coverage'] = pd.to_numeric(df['mean_coverage'], errors='coerce')
    df['dropout(seq)'] = pd.to_numeric(df['dropout(seq)'], errors='coerce')
    df['dropout(decode)'] = pd.to_numeric(df['dropout(decode)'], errors='coerce')
    df['dropout(seq)'] = pd.to_numeric(df['dropout(seq)'], errors='coerce')
    df['dropout(decode)'] = pd.to_numeric(df['dropout(decode)'], errors='coerce')

    df = df.sort_values(by='mean_coverage')

    group_size = 3
    mean_rows = []

    for i in range(0, len(df) - group_size+1, group_size):
        group = df.iloc[i:i+group_size]
        success_count = (group['decode_success'] == 'Yes').sum()
        decode_success = 'Yes' if success_count >= 2 else 'No'
        mean_rows.append({
            'mean_coverage': group['mean_coverage'].mean(),
            'dropout_seq': group['dropout(seq)'].mean(),
            'dropout_decode': group['dropout(decode)'].mean(),
            'decode_success': decode_success
        })
    mean_df = pd.DataFrame(mean_rows)
    
    return df, mean_df

def get_ecc_label(filename):
    if "crc_grand" in filename.lower():
        return "CRC-Grand"
    elif "crc" in filename.lower():
        return "CRC"
    elif "optimized_grand" in filename.lower():
        return "Optimized-GRAND"
    elif "rs" in filename.lower():
        return "RS"
    return "Unknown"

def plot_dropout(mean_df, alpha, subs_rate, output_path, ecc_label):
    plt.figure(figsize=(10, 6))

    plt.plot(mean_df['mean_coverage'], mean_df['dropout_seq'], label='Dropout (seq)', marker='o', color='blue')
    plt.plot(mean_df['mean_coverage'], mean_df['dropout_decode'], label='Dropout (decode)', marker='x', color='purple')

    for i in range(1, len(mean_df)):
        x_vals = mean_df['mean_coverage'].iloc[i-1:i+1]
        y_vals = mean_df['dropout_decode'].iloc[i-1:i+1]
        color = 'green' if mean_df['decode_success'].iloc[i] == 'Yes' else 'purple'
        plt.plot(x_vals, y_vals, marker='x', color=color)


    param_text = f'α = {alpha}\nSubstitution Error Rate = {subs_rate}'
    plt.gca().text(0.70, 0.95, param_text,
                   transform=plt.gca().transAxes,
                   fontsize=10,
                   verticalalignment='top',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray'))


    plt.xlabel('Mean Coverage')
    plt.ylabel('% Oligo Dropout')
    plt.title(f'% Oligo Dropout: Synthesis vs Decode ({ecc_label})')

    plt.xticks(np.arange(0, 11, 1))
    plt.xlim(0, 11)
    plt.yticks(np.arange(0, 111, 10))
    plt.ylim(0, 110)

    success_flag = any(mean_df['decode_success'] == 'Yes')
    handles,labels = plt.gca().get_legend_handles_labels()

    if success_flag:
        success_line = Line2D([0], [0], color='green', linestyle='-', marker='x', label='Decode Successful')
        handles.append(success_line)
        labels.append('Decode Successful')

    plt.legend(handles=handles, labels=labels)

    plt.grid(True)
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path)
    plt.close()

def process_all_files(input_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    csv_files = glob.glob(os.path.join(input_folder, "*.csv"))

    for csv_file in csv_files:
        df, mean_df = clean_data(csv_file)

        alpha = df['α'].iloc[0]
        subs_rate = df['subs_rate'].iloc[0]

        base_name = os.path.splitext(os.path.basename(csv_file))[0] 
        print(f"Basename: {base_name}")
        ecc_label = get_ecc_label(base_name)

        if "coverage_metrics_" in base_name:
            output_label = base_name.split("coverage_metrics_")[1]
        else:
            output_label = base_name

        output_file = output_label + ".png"

        output_path = os.path.join(output_folder, output_file)
        plot_dropout(mean_df, alpha, subs_rate, output_path, ecc_label)

if __name__ == "__main__":
    input_folder = "coverage-analysis/seq-depth/files/oligo_recovery"
    output_folder = "coverage-analysis/seq-depth/visualizations/dropout_vs_coverage"
    process_all_files(input_folder, output_folder)