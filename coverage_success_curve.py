from Analysis.Fountain_analyzer import FT_Analyzer
from Model.Model import DNA_Channel_Model
from Model.config import DEFAULT_PASSER
# import Model.config as config
import matplotlib.pyplot as plt
import numpy as np

# Parameters
file_name = "lena.jpg"       # input file path (should be in /files folder)
rs_length = 4                # RS parity symbols
chunk_size = 20              # bytes per chunk
num_trials = 5              # number of simulations per alpha
alphas = np.arange(0.1, 1.05, 0.05)  # alpha range from 0.0 to 1.0 (=> 1.0x to 2.0x coverage)

# error sim
arg = DEFAULT_PASSER
arg.syn_sub_prob = 0.1 / 3


# Model
Model = DNA_Channel_Model(None, DEFAULT_PASSER)

# Results storage
coverage_levels = []
success_rates = []
avg_reads_per_success = []

for alpha in alphas:
    successes = 0
    read_counts = []

    print(f"\n▶ Running simulations for alpha = {alpha:.2f} (Coverage = {1 + alpha:.2f}x)")
    analyzer = FT_Analyzer(file_name, Model, alpha=alpha, rs_length=rs_length, chunk_size=chunk_size)

    for i in range(num_trials):
        print(f"  Trial {i+1}/{num_trials}...", end=' ')
        analyzer.encode()  # Reset encoder to ensure fresh encoding each time
        analyzer.simu()
        reads_used = analyzer.decode(save=False)

        if reads_used != -1:
            successes += 1
            read_counts.append(reads_used)
            print(f"Success ({reads_used} reads)")
        else:
            print("Failed")

    coverage_levels.append(round(1 + alpha, 2))
    success_rate = successes / num_trials
    success_rates.append(success_rate)
    avg_reads = np.mean(read_counts) if read_counts else None
    avg_reads_per_success.append(avg_reads)

# Plot: Success rate vs. Coverage
plt.figure(figsize=(8,5))
plt.plot(coverage_levels, success_rates, marker='o', color='green')
plt.title("Decoding Success Rate vs. Coverage")
plt.xlabel("Coverage (x)")
plt.ylabel("Success Rate")
plt.grid(True)
plt.ylim(0, 1.05)
plt.tight_layout()
plt.savefig(f"coverage-analysis/coverage_success_curve/{file_name}_er-{arg.syn_sub_prob*3}.png", dpi=300)

# Plot: Average Reads (only for successful decodes)
plt.figure(figsize=(8,5))
plt.plot(coverage_levels, avg_reads_per_success, marker='o', color='blue')
plt.title("Avg. Reads Needed for Successful Decoding vs. Coverage")
plt.xlabel("Coverage (x)")
plt.ylabel("Average Reads Used")
plt.grid(True)
plt.tight_layout()
plt.show()