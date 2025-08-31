import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from matplotlib.lines import Line2D

def clean_data(csv_file):
    df = pd.read_csv(csv_file) 
    df['mean_coverage'] = pd.to_numeric(df['mean_coverage'],errors='coerce')
    df['dropout(seq)'] = pd.to_numeric(df['dropout(seq)'],errors='coerce')
    df['dropout(decode)'] = pd.to_numeric(df['dropout(decode)'],errors='coerce')
    df['dropout(seq)'] = pd.to_numeric(df['dropout(seq)'],errors='coerce')
    df['dropout(decode)'] = pd.to_numeric(df['dropout(decode)'],errors='coerce')

    # df = df.sort_values(by='mean_coverage')

    group_size = 5
    mean_rows = []

    for i in range(0, len(df) - group_size+1, group_size):
        group = df.iloc[i:i+group_size]
        success_count = (group['decode_success'] == 'Yes').sum()
        decode_success = 'Yes' if success_count >= 3 else 'No'
        mean_rows.append({
            'mean_coverage': group['mean_coverage'].mean(),
            'dropout_seq': group['dropout(seq)'].mean(),
            'dropout_decode': group['dropout(decode)'].mean(),
            'decode_success': decode_success
        })
    mean_rows_df = pd.DataFrame(mean_rows)
    # print(mean_rows_df)
    
    return df, mean_rows_df

def get_ecc_label(filename):
    if "rs" in filename.lower():
        return "RS"
    elif "crc-only" in filename.lower():
        return "CRC-Only"
    elif "bitwise" and "bruteforce" in filename.lower():
        return "Bit-wise Bruteforce"
    elif "bitwise" and "topk" in filename.lower():
        return "Bit-wise heuristic with Top K"
    elif "bitwise" and "nok" in filename.lower():
        return "Bit-wise heuristic without K limit"
    elif "bitwise" and "ends" in filename.lower():
        return "Bit-wise heuristic edge indices"
    elif "basewise" and "bruteforce" in filename.lower():
        return "Base-wise Bruteforce"
    elif "basewise" and "topk" in filename.lower():
        return "Base-wise heuristic with Top K"
    elif "basewise" and "nok" in filename.lower():
        return "Base-wise heuristic without K limit"
    elif "basewise" and "ends" in filename.lower():
        return "Base-wise heuristic edge indices"
    return "Unknown"

def plot_dropout(mean_rows_df, alpha, subs_rate, output_path, ecc_label):
    print(ecc_label)
    plt.figure(figsize=(10, 6))

    plt.plot(mean_rows_df['mean_coverage'], mean_rows_df['dropout_seq'], label='Dropout (seq)', marker='o', color='blue')
    plt.plot(mean_rows_df['mean_coverage'], mean_rows_df['dropout_decode'], label='Dropout (decode)', marker='x', color='purple')

    for i in range(1, len(mean_rows_df)):
        x_values = mean_rows_df['mean_coverage'].iloc[i-1:i+1]
        # print(f"\nX: {x_values}")
        y_values = mean_rows_df['dropout_decode'].iloc[i-1:i+1]
        # print(f"Y: {y_values}")
        color = 'green' if mean_rows_df['decode_success'].iloc[i-1] == 'Yes' else 'purple'
        # if mean_rows_df['decode_success'].iloc[i-1] == 'Yes':
        #     print(f"{mean_rows_df['decode_success'].iloc[i-1]} vlaue is yes")
        plt.plot(x_values, y_values, marker='x', color=color)


    # param_text = f'α = {alpha}\nSubstitution Error Rate = {subs_rate}'
    # plt.gca().text(0.70, 0.95, param_text,
    #                transform=plt.gca().transAxes,
    #                fontsize=14,
    #                verticalalignment='top',
    #                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray'))


    plt.xlabel('Mean Coverage', fontsize=16)
    plt.ylabel('% Oligo Dropout', fontsize=16)
    # plt.title(f'% Oligo Dropout: Synthesis vs Decode ({ecc_label})')
    # plt.title(f'{ecc_label}')

    plt.xticks(np.arange(0, 11, 1), fontsize=12)
    plt.xlim(0, 11)
    plt.yticks(np.arange(0, 111, 10), fontsize=12)
    plt.ylim(0, 110)

    success_flag = any(mean_rows_df['decode_success'] == 'Yes')
    handles,labels = plt.gca().get_legend_handles_labels()

    if success_flag:
        success_line = Line2D([0], [0], color='green', linestyle='-', marker='x', label='Decode Successful')
        handles.append(success_line)
        labels.append('Decode Successful')

    plt.legend(handles=handles, labels=labels, fontsize=12)

    plt.grid(True)
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path)
    plt.close()

def parse_all_files(input_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    csv_files = glob.glob(os.path.join(input_folder, "*.csv"))

    for csv_file in csv_files:
        df, mean_rows_df = clean_data(csv_file)

        alpha = df['α'].iloc[0]
        subs_rate = df['subs_rate'].iloc[0]

        base_name = os.path.splitext(os.path.basename(csv_file))[0] 
        print(f"Basename: {base_name}")
        ecc_label = get_ecc_label(base_name)

        if "coverage_metrics_" in base_name:
            output_label=base_name.split("coverage_metrics_")[1]
        else:
            output_label=base_name

        output_file=output_label+".pdf"

        output_path=os.path.join(output_folder, output_file)
        plot_dropout(mean_rows_df, alpha, subs_rate, output_path, ecc_label)

if __name__ == "__main__":
    input_folder = "IO/Output/final_csvs"
    output_folder = "IO/Output/coverage_images"
    parse_all_files(input_folder, output_folder)