import numpy as np
import matplotlib.pyplot as plt
import datetime
import time
from Encode.DNAFountain import DNAFountain, Glass
from Encode.Helper_Functions import preprocess

def run_multiple_coverage_experiments(file_path, alpha=0.5, rs=4, n_runs=30):
    # Step 1: Preprocess input file into chunks
    data, pad = preprocess(file_path, 20)
    N = len(data)
    all_curves = []

    for run in range(n_runs):
        # Step 2: Encode with DNA Fountain
        f = DNAFountain(data, alpha, rs=rs)
        f.encode()
        dna_file = f"{file_path}_simu_run{run}.dna"
        f.save(dna_file)

        # Step 3: Decode and collect coverage data
        g = Glass(dna_file, chunk_num=N, rs=rs)
        ret, solve_num, lineRead, chunksDone, errors, coverage_vs_reads = g.decode()
        
        # Pad with last value to align lengths
        while len(coverage_vs_reads) < 600:
            coverage_vs_reads.append(coverage_vs_reads[-1])
        all_curves.append(coverage_vs_reads)

    # Step 4: Compute average coverage across runs
    avg_coverage = np.mean(all_curves, axis=0)

    # Step 5: Plot
    plt.figure(figsize=(8, 4))
    for curve in all_curves[:5]:  # plot a few sample runs
        plt.plot(curve, color='gray', alpha=0.3)
    plt.plot(avg_coverage, color='blue', label='Average coverage')
    plt.xlabel("Read count")
    plt.ylabel("Unique oligos covered")
    plt.title(f"Average Coverage Curve (α={alpha}, n={n_runs})")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    # plt.show()
    # timestamp = datetime.datetime.now().strftime("%d-%m-%Y_%H:%M:%S")
    # plt.savefig(f"coverage_vs_reads/{file_path[6:]}_{timestamp}.png", dpi=300)
    plt.savefig(f"coverage_vs_reads/{file_path[6:]}_{n_runs}.png", dpi=300)

start_time = time.time()
run_multiple_coverage_experiments("files/lena.jpg", alpha=0.5, rs=4, n_runs=500)
print("Run time: %s seconds" % (time.time() - start_time))