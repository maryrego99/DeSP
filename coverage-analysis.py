# per oligo analysis
import numpy as np
import matplotlib.pyplot as plt
import datetime
import time
from math import exp
from Encode.DNAFountain import DNAFountain, Glass
from Encode.Helper_Functions import preprocess

def run_multiple_coverage_experiments(file_path, alpha=0.5, rs=4, n_runs=100, max_reads=600):
    # Step 1: Preprocess input file into chunks
    data, pad = preprocess(file_path, 20)
    N = len(data)
    
    all_curves = []
    read_counts = []
    dropout_counts = []
    failed_runs = 0

    for run in range(n_runs):
        # Step 2: Encode with DNA Fountain
        f = DNAFountain(data, alpha, rs=rs)
        f.encode()
        dna_file = f"{file_path}_simu_run{run}.dna"

        # # simulate droput
        # # Save only a fraction of DNA strands to simulate poor sampling
        # all_dna = [d[0] for d in f.dna_dl]  # Get list of DNA strings
        # np.random.seed(run)  # Ensure reproducibility across runs
        # sampled_dna = np.random.choice(all_dna, size=int(0.7 * len(all_dna)), replace=False)  # 60% sampling
        # dna_file = f"{file_path}_simu_run{run}.dna"
        # with open(dna_file, "w") as out:
        #     out.writelines(d + "\n" for d in sampled_dna)

        f.save(dna_file)

        # Step 3: Decode and collect coverage data
        g = Glass(dna_file, chunk_num=N, rs=rs)
        ret, solve_num, lineRead, chunksDone, errors, coverage_vs_reads, chunk_seen = g.decode()
        if ret == 0: 
            print("Decoding succeeded!")
        else:
            print("'Decoding Failed-.-'")

        # Pad coverage curve if needed
        while len(coverage_vs_reads) < max_reads:
            coverage_vs_reads.append(coverage_vs_reads[-1])

        all_curves.append(coverage_vs_reads)

        # Read count to full coverage
        if N in coverage_vs_reads:
            read_counts.append(coverage_vs_reads.index(N))
        else:
            failed_runs += 1

        # # Dropout: how many chunks never seen
        # seen_chunks = set()
        # for i, val in enumerate(coverage_vs_reads):
        #     if i == 0 or val > coverage_vs_reads[i - 1]:
        #         seen_chunks.update(range(coverage_vs_reads[i - 1] if i > 0 else 0, val))
        # dropout_counts.append(N - len(seen_chunks))
        # print(f"dropout counts: {chunk_seen.count(0)}")
        dropout_counts.append(chunk_seen.count(0))

    # Average curve
    avg_coverage = np.mean(all_curves, axis=0)

    # Theoretical curve
    theoretical_reads = np.arange(1, len(avg_coverage) + 1)
    expected_coverage = [N * (1 - exp(-r / N)) for r in theoretical_reads]

    # Plot coverage
    plt.figure(figsize=(10, 5))
    plt.plot(avg_coverage, label='Average coverage', color='blue')
    plt.plot(theoretical_reads, expected_coverage, '--', color='green', label='Theoretical n·ln(n)')
    plt.xlabel("Read count")
    plt.ylabel("Unique oligos covered")
    plt.title(f"Average Coverage Curve (α={alpha}, n={n_runs})")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    # plt.show()
    # timestamp = datetime.datetime.now().strftime("%d-%m-%Y_%H:%M:%S")
    # plt.savefig(f"coverage_vs_reads/{file_path[6:]}_{timestamp}.png", dpi=300)
    plt.savefig(f"coverage-analysis/coverage_vs_reads/{file_path[6:]}_{n_runs}.png", dpi=300)

    # Plot dropout histogram
    plt.figure()
    plt.hist(dropout_counts, bins=range(0, max(dropout_counts)+2), color='red', alpha=0.6)
    plt.title("Dropout Distribution Across Runs")
    plt.xlabel("Number of Chunks Never Seen")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.show()
    plt.savefig(f"coverage-analysis/dropout-dist/{file_path[6:]}_{n_runs}_histogram.png", dpi=300)

    # Print summary
    if read_counts:
        avg_reads = np.mean(read_counts)
        success_rate = (len(read_counts) / n_runs) * 100
    else:
        avg_reads = float('nan')
        success_rate = 0.0

    avg_dropouts = np.mean(dropout_counts)
    print(f"Avg reads to full coverage: {avg_reads:.2f}")
    print(f"Success rate: {success_rate:.1f}% ({n_runs - failed_runs} / {n_runs})")
    print(f"Avg dropouts: {avg_dropouts:.2f}")

# Example call:
# run_coverage_analysis("files/lena.jpg", alpha=0.5, rs=4, n_runs=100)

start_time = time.time()
run_multiple_coverage_experiments("files/lena.jpg", alpha=0.5, rs=4, n_runs=10)
print("Run time: %s seconds" % (time.time() - start_time))