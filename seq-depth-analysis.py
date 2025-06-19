import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from Analysis.Fountain_analyzer import FT_Analyzer
from Model.Model import Synthesizer, Decayer, PCRer, Sampler, Sequencer
from Model.config import DEFAULT_PASSER
from Encode.Helper_Functions import preprocess
from Encode.DNAFountain import DNAFountain, Glass
from Analysis.Analysis import save_simu_result


def log_coverage_metrics(input_file, seq_counts_file, total_oligos, params, decoded_success, out_csv="coverage_metrics.csv"):
    with open(seq_counts_file, "r") as f:
        next(f)  #exclude the header
        counts = [int(line.strip().split(",")[1]) for line in f]

    total_reads = sum(counts)
    nonzero_oligos = sum(1 for c in counts if c > 0)
    dropout_oligos = sum(1 for c in counts if c == 0)

    mean_coverage = total_reads/total_oligos if total_oligos > 0 else 0
    percent_seen = 100 * nonzero_oligos/total_oligos if total_oligos > 0 else 0
    dropout_rate = 100 * dropout_oligos/total_oligos if total_oligos > 0 else 0

    row = {
        "input_file": input_file,
        "file": os.path.basename(seq_counts_file),
        "alpha": params.get("alpha"),
        "pcrc": params.get("pcrc"),
        "pcrp": params.get("pcrp"),
        "sampling_ratio": params.get("sam_ratio"),
        "subs_rate": params.get("subs_rate"),
        "rs": params.get("rs"),
        "total_oligos": total_oligos,
        "total_reads": total_reads,
        "mean_coverage": round(mean_coverage, 2),
        "%_oligos_seen": round(percent_seen, 2),
        "dropout_rate": round(dropout_rate, 2),
        "decoded_successfully": "Yes" if decoded_success else "No"
    }

    file_exists = os.path.isfile(out_csv)
    with open(out_csv, "a", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=row.keys())
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

    print(f"Coverage metrics logged to {out_csv}")

def analyze_oligo_coverage(file_path, alpha, subs_rate):
    rs = 0
    chunk_size = 20 # increase -> lesser oligos

    #input file into chunks
    data, pad = preprocess(file_path, chunk_size)
    print(file_path, ' loaded and split into ', len(data), ' data chunks.')
    N = len(data)  #no. of chunks

    #simulation args
    arg = DEFAULT_PASSER
    arg.syn_number = 30
    arg.syn_sub_prob = subs_rate / 3
    arg.syn_yield = 0.99
    arg.pcrc = 14
    arg.pcrp = 0.9 #shouldnt change
    arg.sam_ratio = 0.9
    arg.seq_depth = 10 # vary ex:1-5, different recovery rates
    from Model.config import TM_NGS, TM_NNP
    arg.seq_TM = TM_NGS # TM_NNP

    #encoding
    f = DNAFountain(data, alpha, rs=rs)
    good, tries = f.encode() #with rs=0, the encoded.dna file is shorter than files/lena.dna as that had rs encoded
    dna_file = f"{file_path}_encoded.dna"
    f.save(dna_file)
    print('Data encoded into ' ,good, ' DNA strands after ', tries, ' tries.')
    print('Saved to ', dna_file)
    # designed_oligos = [entry[0] for entry in f.dna_dl]

    #Error simulation
    with open(dna_file) as file:
        dnas = file.readlines()
    # in_dnas = [dna.strip() for dna in dnas]
    in_dnas = [dna.split('\n')[0] for dna in dnas]
    print(dna_file, ' loaded: ', len(in_dnas), ' strands of length ', len(in_dnas[0]))

    dnas_syn = Synthesizer(arg)(in_dnas)
    with open("coverage-analysis/seq-depth/files/syn_dna_counts.csv", "w") as f:
        f.write("oligo_index,total_copies\n")
        for i, barcode in enumerate(dnas_syn):
            f.write(f"{i},{barcode['num']}\n")
    dnas_dec = Decayer(arg)(dnas_syn)
    dnas_pcr = PCRer(N=12, p=0.8)(dnas_dec)

    with open("coverage-analysis/seq-depth/files/pcr_copy_counts.csv", "w") as f:
        f.write("oligo_index,total_copies\n")
        for i, barcode in enumerate(dnas_pcr):
            f.write(f"{i},{barcode['num']}\n")

    #to check oligo copies
    # pcr_out_file = f"{file_path}_after_pcr.dna"
    # save_simu_result(dnas_pcr, pcr_out_file)

    dnas_sam = Sampler(p=arg.sam_ratio)(dnas_pcr)
    dnas_seq = Sequencer(arg)(dnas_sam)
    # print(f"dnas_seq: {dnas_seq}")
    with open("coverage-analysis/seq-depth/files/seq_copy_counts.csv", "w") as f:
        f.write("oligo_index,total_copies\n")
        for i, barcode in enumerate(dnas_seq):
            f.write(f"{i},{barcode['num']}\n")

    noisy_dna_file = f"{file_path}_errors.dna"
    save_simu_result(dnas_seq, noisy_dna_file)
    print('Simulation results saved to ', noisy_dna_file)


    #Decoding
    print('Trying to decode from sequencing readouts.')
    g = Glass(noisy_dna_file, chunk_num=N, rs=rs)
    ret, _, total_reads, chunks_done, errors, coverage_vs_reads, chunk_seen = g.decode()
    decoded_file = f"{file_path}_decoded.jpg"
    if ret ==0:
        print("Decoding successful")
        g.save(decoded_file)
    else:
        print("Decoding failed.")
    
    # crc_pass = g.crc_pass
    # crc_fail = g.crc_fail

    decoded_success = (ret == 0)

    params = {
        "alpha": alpha,
        "pcrc": arg.pcrc,
        "pcrp": arg.pcrp,
        "sam_ratio": arg.sam_ratio,
        "subs_rate": subs_rate,
        "rs": rs,
    }

    log_coverage_metrics(
        input_file = file_path.split("files/")[-1],
        seq_counts_file="coverage-analysis/seq-depth/files/seq_copy_counts.csv",
        total_oligos=len(in_dnas),
        params=params,
        decoded_success = decoded_success,
        out_csv="coverage-analysis/seq-depth/files/coverage_metrics.csv"
    )


if __name__ == "__main__":
    analyze_oligo_coverage("coverage-analysis/seq-depth/files/lena.jpg", alpha=0.7, subs_rate=0.003)
