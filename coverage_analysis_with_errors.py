
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from math import exp
from Encode.DNAFountain import DNAFountain, Glass
from Encode.Helper_Functions import preprocess
from Model.config import DEFAULT_PASSER
from Model.Model import Synthesizer, Decayer, PCRer, Sampler, Sequencer
from Analysis.Analysis import dna_chunk, plot_oligo_number_distribution, plot_error_distribution, save_simu_result

def run_multiple_coverage_experiments(file_path, alpha=0.5, rs=4, n_runs=30, max_reads=600):
    from copy import deepcopy
    from Model.config import DEFAULT_PASSER

    # Make a safe copy so the global DEFAULT_PASSER isn't affected
    arg = deepcopy(DEFAULT_PASSER)

    # Hardcoded simulation parameters
    arg.syn_number = 30             # Number of strands synthesized per oligo
    arg.syn_sub_prob = 0.05 / 3     # Substitution error rate (divided across 3 types)
    arg.syn_yield = 0.98            # Synthesis yield

    arg.pcrc = 12                   # PCR cycles
    arg.pcrp = 0.8                  # PCR efficiency

    arg.sam_ratio = 0.005           # Sampling ratio
    arg.seq_depth = 5               # Sequencing depth

    # Choose sequencing platform manually (NGS = Illumina, NNP = Nanopore)
    from Model.config import TM_NGS
    arg.seq_TM = TM_NGS             # or: TM_NNP

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
        dnas_pcr = PCRer(N=12, p=0.8)(dnas_dec)
        print("Sampling probability =", arg.sam_ratio, type(arg.sam_ratio))
        dnas_sam = Sampler(p=0.005)(dnas_pcr)
        dnas_seq = Sequencer(arg)(dnas_sam)

        noisy_dna_file = f"{file_path}_simu_noisy_run{run}.dna"
        save_simu_result(dnas_seq,noisy_dna_file)

        # Decode
        g = Glass(noisy_dna_file, chunk_num=N, rs=rs)
        ret, solve_num, lineRead, chunksDone, errors, coverage_vs_reads, chunk_seen = g.decode()
        if ret == 0: 
            print('Decoding succeeded!')
        else:
            print('Decoding Failed-.-')

        if coverage_vs_reads:
            while len(coverage_vs_reads) < max_reads:
                coverage_vs_reads.append(coverage_vs_reads[-1])
        else:
            coverage_vs_reads = [0] * max_reads  # or whatever fallback makes sense
        all_curves.append(coverage_vs_reads)

        if N in coverage_vs_reads:
            read_counts.append(coverage_vs_reads.index(N))
        else:
            failed_runs += 1

        dropout_counts.append(chunk_seen.count(0))

    avg_coverage = np.mean(all_curves, axis=0)
    theoretical_reads = np.arange(1, len(avg_coverage) + 1)
    expected_coverage = [N * (1 - exp(-r / N)) for r in theoretical_reads]

    # plt.figure(figsize=(10, 5))
    # plt.plot(avg_coverage, label='Average coverage', color='blue')
    # plt.plot(theoretical_reads, expected_coverage, '--', color='green', label='Theoretical n·ln(n)')
    # plt.xlabel("Read count")
    # plt.ylabel("Unique oligos covered")
    # plt.title(f"Average Coverage Curve (α={alpha}, n={n_runs})")
    # plt.legend()
    # plt.grid(True)
    # plt.tight_layout()
    # os.makedirs("coverage-analysis/coverage_vs_reads", exist_ok=True)
    # plt.savefig(f"coverage-analysis/coverage_vs_reads/{file_path[6:]}_{n_runs}.png", dpi=300)

    plt.figure(figsize=(10, 5))
    plt.plot(avg_coverage, label='Average coverage', color='blue')

    # Theoretical curve
    plt.plot(theoretical_reads, expected_coverage, '--', color='green', label='Theoretical n·ln(n)')

    # Add horizontal line for N
    plt.axhline(y=N, color='gray', linestyle='--', label='Full Coverage Target (N)')

    # Mark full coverage if achieved
    if N in avg_coverage:
        idx = list(avg_coverage).index(N)
        plt.scatter(idx, N, color='red', label=f'Full Coverage at {idx} reads')

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
    run_multiple_coverage_experiments("files/lena.jpg", alpha=0.7, rs=4, n_runs=16)
    print("Run time: %s seconds" % (time.time() - start_time))
