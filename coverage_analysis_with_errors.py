
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from math import exp
from Encode.DNAFountain import DNAFountain, Glass
from Encode.Helper_Functions import preprocess
from Model.config import DEFAULT_PASSER
from Model.Model import Synthesizer, Decayer, PCRer, Sampler, Sequencer

def run_multiple_coverage_experiments(file_path, alpha=0.5, rs=4, n_runs=30, max_reads=600):
    # Use the same parameters as in Streamlit
    arg = DEFAULT_PASSER
    arg.syn_number = 30
    arg.syn_sub_prob = 0.05 / 3
    arg.syn_yield = 0.98
    arg.pcrc = 12
    arg.pcrp = 0.8
    arg.sam_ratio = 0.001
    arg.seq_depth = 5
    arg.seq_TM = arg.seq_TM  # already defaulted to NGS or NNP

    data, pad = preprocess(file_path, 20)
    N = len(data)

    all_curves = []
    read_counts = []
    dropout_counts = []
    failed_runs = 0

    for run in range(n_runs):
        print(f"Run {run+1}/{n_runs}...")
        f = DNAFountain(data, alpha, rs=rs)
        f.encode()
        dna_file = f"{file_path}_simu_run{run}.dna"
        f.save(dna_file)

        # Error simulation
        with open(dna_file) as file:
            dnas = file.readlines()
        in_dnas = [dna.strip() for dna in dnas]

        dnas_syn = Synthesizer(arg)(in_dnas)
        dnas_dec = Decayer(arg)(dnas_syn)
        dnas_pcr = PCRer(arg)(dnas_dec)
        dnas_sam = Sampler(arg)(dnas_pcr)
        dnas_seq = Sequencer(arg)(dnas_sam)

        noisy_dna_file = f"{file_path}_simu_noisy_run{run}.dna"
        with open(noisy_dna_file, "w") as f_out:
            f_out.writelines("\n".join(dnas_seq))

        # Decode
        g = Glass(noisy_dna_file, chunk_num=N, rs=rs)
        ret, solve_num, lineRead, chunksDone, errors, coverage_vs_reads, chunk_seen = g.decode()

        while len(coverage_vs_reads) < max_reads:
            coverage_vs_reads.append(coverage_vs_reads[-1])
        all_curves.append(coverage_vs_reads)

        if N in coverage_vs_reads:
            read_counts.append(coverage_vs_reads.index(N))
        else:
            failed_runs += 1

        dropout_counts.append(chunk_seen.count(0))

    avg_coverage = np.mean(all_curves, axis=0)
    theoretical_reads = np.arange(1, len(avg_coverage) + 1)
    expected_coverage = [N * (1 - exp(-r / N)) for r in theoretical_reads]

    plt.figure(figsize=(10, 5))
    plt.plot(avg_coverage, label='Average coverage', color='blue')
    plt.plot(theoretical_reads, expected_coverage, '--', color='green', label='Theoretical n·ln(n)')
    plt.xlabel("Read count")
    plt.ylabel("Unique oligos covered")
    plt.title(f"Average Coverage Curve (α={alpha}, n={n_runs})")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    os.makedirs("coverage-analysis/coverage_vs_reads", exist_ok=True)
    plt.savefig(f"coverage-analysis/coverage_vs_reads/{file_path[6:]}_{n_runs}.png", dpi=300)

    plt.figure()
    plt.hist(dropout_counts, bins=range(0, max(dropout_counts)+2), color='red', alpha=0.6)
    plt.title("Dropout Distribution Across Runs")
    plt.xlabel("Number of Chunks Never Seen")
    plt.ylabel("Frequency")
    plt.tight_layout()
    os.makedirs("coverage-analysis/dropout-dist", exist_ok=True)
    plt.savefig(f"coverage-analysis/dropout-dist/{file_path[6:]}_{n_runs}_histogram.png", dpi=300)

    avg_reads = np.mean(read_counts) if read_counts else float('nan')
    success_rate = (len(read_counts) / n_runs) * 100
    avg_dropouts = np.mean(dropout_counts)

    print(f"Avg reads to full coverage: {avg_reads:.2f}")
    print(f"Success rate: {success_rate:.1f}% ({n_runs - failed_runs} / {n_runs})")
    print(f"Avg dropouts: {avg_dropouts:.2f}")

if __name__ == "__main__":
    start_time = time.time()
    run_multiple_coverage_experiments("files/lena.jpg", alpha=0.25, rs=4, n_runs=20)
    print("Run time: %s seconds" % (time.time() - start_time))
