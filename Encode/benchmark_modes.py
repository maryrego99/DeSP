import os
import time
import logging
from DNAFountain import DNAFountain, Glass

# Set up logging (optional)
logging.basicConfig(level=logging.INFO)

# ----------------- INPUT CONFIG ----------------- #
INPUT_FILE = "dummy.txt"
CHUNK_SIZE = 32  # Adjust to match your config

with open(INPUT_FILE, "rb") as f:
    raw = f.read()

# Pad if needed
pad = (CHUNK_SIZE - len(raw) % CHUNK_SIZE) % CHUNK_SIZE
raw += b'\0' * pad

# Split into chunks
chunks = [raw[i:i+CHUNK_SIZE] for i in range(0, len(raw), CHUNK_SIZE)]

# ----------------- BENCHMARK FUNCTION ----------------- #
def benchmark_mode(rs_flag, description):
    print(f"\n==== Benchmark: {description} (rs = {rs_flag}) ====")
    
    # ------------------- Encoding ------------------- #
    encoder = DNAFountain(
        file_in=chunks,
        alpha=0.5,
        rs=rs_flag
    )

    start = time.time()
    encoder.encode()
    encoder.save("test_out.dna")
    encode_time = time.time() - start
    print(f"[✓] Encoding complete in {encode_time:.2f} sec")

    # ------------------- Decoding ------------------- #
    decoder = Glass(
        in_file_name="test_out.dna",
        chunk_num=len(chunks),
        rs=rs_flag,
        chunk_size=CHUNK_SIZE
    )

    ret, solve_num, line_count, chunks_done, errors = decoder.decode()
    print(f"[✓] Decoding result: {chunks_done}/{len(chunks)} chunks recovered")
    print(f"[✓] Lines read: {line_count} | Errors: {errors}")

    return {
        "mode": description,
        "chunks_recovered": chunks_done,
        "lines_read": line_count,
        "errors": errors,
        "encode_time": encode_time,
        "success": chunks_done == len(chunks)
    }

# ----------------- MODES TO TEST ----------------- #
results = []

results.append(benchmark_mode(0, "Fountain only"))
results.append(benchmark_mode(-1, "LDPC only"))
results.append(benchmark_mode(10, "RS only (10 symbols)"))  # adjust as needed
results.append(benchmark_mode(-2, "Hybrid (LDPC + RS)"))

# ----------------- PRINT SUMMARY ----------------- #
print("\n===== Summary =====")
for r in results:
    print(f"{r['mode']:20} | Chunks: {r['chunks_recovered']}/{len(chunks)} | Errors: {r['errors']} | Time: {r['encode_time']:.2f}s | {'✅' if r['success'] else '❌'}")