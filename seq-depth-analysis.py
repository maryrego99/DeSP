import os
import csv
import time
import numpy as np
import matplotlib.pyplot as plt
from Analysis.Fountain_analyzer import FT_Analyzer
from Model.Model import Synthesizer, Decayer, PCRer, Sampler, Sequencer
from Model.config import DEFAULT_PASSER, TM_NGS
from Encode.Helper_Functions import preprocess
from DNAFountain.Fountain import DNAFountain
from DNAFountain.Glass import Glass
from ECC.NoECCDecoder import NoECCDecoder
from ECC.RSDecoder import ReedSolomonDecoder
from ECC.CRCDecoder import CRCDecoder
from ECC.CRCGrandDecoder import CRCGrandDecoder
from ECC.ecc_encoders import crc32_encoder, make_rs_encoder, no_encoder
from Analysis.Analysis import save_simu_result


def log_metrics_to_csv(row, out_csv):
    file_exists = os.path.isfile(out_csv)
    with open(out_csv, "a", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=row.keys())
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)
    print(f"Logged metrics to {out_csv}")


def log_data_recovery_metrics(file_path, total_chunks, recovered_chunks, out_csv):
    recovery_rate = 100 * recovered_chunks / total_chunks if total_chunks > 0 else 0
    row = {
        "input_file": file_path,
        "total_chunks": total_chunks,
        "recovered_chunks": recovered_chunks,
        "recovery_rate(%)": round(recovery_rate, 2),
        "dropout_rate(%)": round(100 - recovery_rate, 2)
    }
    log_metrics_to_csv(row, out_csv)


def plot_data_chunk_recovery(chunk_seen, save_path):
    if not chunk_seen or len(chunk_seen) == 0:
        print("Warning: chunk_seen is empty. Skipping plot.")
        return
    plt.figure(figsize=(10, 2))
    plt.imshow([chunk_seen], cmap='Greens', aspect='auto')
    plt.xlabel("Chunk Index")
    plt.yticks([])
    plt.title("Recovered Chunks (green = recovered)")
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
    print(f"Chunk recovery plot saved to {save_path}")


def log_coverage_metrics(input_file, seq_counts_file, total_oligos, params, decoded_success, decode_time, out_csv):
    with open(seq_counts_file, "r") as f:
        next(f)  #exclude the header
        counts = [int(line.strip().split(",")[1]) for line in f]

    total_reads = sum(counts)
    nonzero_oligos = sum(1 for c in counts if c > 0)
    dropout_oligos = sum(1 for c in counts if c == 0)

    mean_coverage = total_reads/total_oligos if total_oligos > 0 else 0
    percent_seen = 100 * nonzero_oligos/total_oligos if total_oligos > 0 else 0
    dropout_rate = 100 * dropout_oligos/total_oligos if total_oligos > 0 else 0

    # decode_recovery_rate = 100 * params.get("crc_pass")/total_oligos if total_oligos > 0 else 0
    # decode_dropout_rate = 100 * params.get("crc_fail")/total_oligos if total_oligos > 0 else 0

    decode_recovery_rate = 100 * params.get("oligos_used_rs_decode", 0) / total_oligos if total_oligos > 0 else 0
    decode_dropout_rate = 100 - decode_recovery_rate
    # print(f"Post sequencing - Non zero oligos: {nonzero_oligos}, Dropout oligos: {dropout_oligos}, Total Oligos: {total_oligos}")

    row = {
        "input_file": input_file,
        "α": params.get("alpha"),
        "sam_ratio": params.get("sam_ratio"),
        "subs_rate": params.get("subs_rate"),
        "total_oligos": total_oligos,
        "total_reads": total_reads,
        "mean_coverage": round(mean_coverage, 2),
        "oligo_recovery(seq)": round(percent_seen, 2),
        "dropout(seq)": round(dropout_rate, 2),
        "oligo_recovery(decode)": round(decode_recovery_rate, 2),
        "dropout(decode)": round(100 - decode_recovery_rate, 2),
        "decode_success": "Yes" if decoded_success else "No",
        "decode_time": round(decode_time, 2)
    }
    log_metrics_to_csv(row, out_csv)


def log_grand_coverage_metrics(input_file, seq_counts_file, total_oligos, params, decoded_success,
                               grand_pass, grand_fail, chunks_recovered, total_chunks, decode_time, out_csv):
    
    with open(seq_counts_file, "r") as f:
        next(f)  #exclude the header
        counts = [int(line.strip().split(",")[1]) for line in f]

    total_reads = sum(counts)
    nonzero_oligos = sum(1 for c in counts if c > 0)
    dropout_oligos = sum(1 for c in counts if c == 0)

    mean_coverage = total_reads/total_oligos if total_oligos > 0 else 0
    percent_seen = 100 * nonzero_oligos/total_oligos if total_oligos > 0 else 0
    dropout_rate = 100 * dropout_oligos/total_oligos if total_oligos > 0 else 0

    decode_recovery_rate = 100 * params.get("oligos_used_rs_decode", 0) / total_oligos if total_oligos > 0 else 0
    decode_dropout_rate = 100 - decode_recovery_rate
    
    grand_success = 100 * grand_pass / (grand_pass + grand_fail) if (grand_pass + grand_fail) else 0
    data_recovered = 100 * chunks_recovered / total_chunks if total_chunks else 0

    row = {
        "input_file": input_file,
        "α": params.get("alpha"),
        "sam_ratio": params.get("sam_ratio"),
        "subs_rate": params.get("subs_rate"),
        "rs": params.get("rs"),
        "total_oligos": total_oligos,
        "total_reads": total_reads,
        "mean_coverage": round(mean_coverage, 2),
        "oligo_recovery(seq)": round(percent_seen, 2),
        "dropout(seq)": round(dropout_rate, 2),
        "oligo_recovery(decode)": round(decode_recovery_rate, 2),
        "dropout(decode)": round(decode_dropout_rate, 2),
        "GRAND_pass_oligos": grand_pass,
        "GRAND_fail_oligos": grand_fail,
        "GRAND_success_rate": round(grand_success, 2),
        "chunks_recovered": chunks_recovered,
        "total_chunks": total_chunks,
        "data_recovery": round(data_recovered, 2),
        "decode_success": "Yes" if decoded_success else "No",
        "decode_time": round(decode_time, 2)
    }
    log_metrics_to_csv(row, out_csv)


def analyze_oligo_coverage(file_path, alpha, ecc_type="CRC_GRAND", subs_rate=0.003, seq_depth=10):
    ecc_label = ecc_type.lower()
    decoder_map = {
        "none": NoECCDecoder(),
        "rs": ReedSolomonDecoder(rs_len=4),
        "crc-only": CRCDecoder(),
        "crc_grand": CRCGrandDecoder(max_flips=2)
    }
    encoder_map = {
        "none": no_encoder,
        "rs": make_rs_encoder(rs_len=4),
        "crc-only": crc32_encoder,
        "crc_grand": crc32_encoder
    }

    decoder = decoder_map[ecc_label]
    encoder = encoder_map[ecc_label]

    chunk_size = 20
    data, pad = preprocess(file_path, chunk_size)
    print(file_path, ' loaded and split into ', len(data), ' data chunks.')
    N = len(data)  #no. of chunks

    #simulation args
    arg = DEFAULT_PASSER
    arg.syn_number = 30
    arg.syn_sub_prob = subs_rate / 3
    # arg.syn_ins_prob = 0.1s
    # arg.syn_del_prob = 0.003
    arg.syn_yield = 0.99
    arg.seq_depth = seq_depth
    arg.seq_TM = TM_NGS

    #encoding
    f = DNAFountain(data, alpha, ecc_encoder=encoder)
    good, tries = f.encode()
    dna_file = f"{file_path}_encoded_{ecc_label}_a{alpha}.dna"
    f.save(dna_file)
    print('Data encoded into ' ,good, ' DNA strands after ', tries, ' tries.')
    print('Saved to ', dna_file)
    # designed_oligos = [entry[0] for entry in f.dna_dl]

    #Error simulation
    with open(dna_file) as file:
        dnas = file.readlines()
    # in_dnas = [dna.strip() for dna in dnas]
    in_dnas = [dna.split('\n')[0] for dna in dnas]
    oligos_generated = len(in_dnas)
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

    seq_file = f"coverage-analysis/seq-depth/files/seq_copy_counts_{ecc_label}_a{alpha}.csv"
    with open(seq_file, "w") as f:
        f.write("oligo_index,total_copies\n")
        for i, barcode in enumerate(dnas_seq):
            f.write(f"{i},{barcode['num']}\n")

    noisy_dna_file = f"{file_path}_errors_{ecc_label}_a{alpha}.dna"
    save_simu_result(dnas_seq, noisy_dna_file)
    print('Simulation results saved to ', noisy_dna_file)


    #Decoding
    print('Trying to decode from sequencing readouts.')
    start_decode_time = time.time()
    g = Glass(noisy_dna_file, chunk_num=N, ecc_decoder=decoder)
    ret, _, _, chunks_done, _, _, chunk_seen, chunks = g.decode()
    decode_time = time.time() - start_decode_time
    print(f"Decoding time: {decode_time}")

    if ret == 0:
        g.save(f"{file_path}_decoded_{ecc_label}_a{alpha}.jpg")
    elif chunks_done > 0:
        g.save_partial_if_valid("partial_output", chunks, pad=0)
    else:
        print("No chunks recovered; nothing to save.")
        print("Decoding failed.")

    used_oligos = len(g.seen_seeds)

    recovered_oligo_ratio = g.recovered_droplets/oligos_generated * 100
    print(f"% Recovered Droplets(Oligos): {recovered_oligo_ratio:.2f}")
    decoded_success = (ret == 0)

    params = {
        "alpha": alpha,
        "sam_ratio": arg.sam_ratio,
        "subs_rate": subs_rate,
        "oligos_used_rs_decode": used_oligos
    }

    short_file = os.path.basename(file_path)

    if ecc_label == "crc_grand":
        log_grand_coverage_metrics(
            input_file=short_file,
            seq_counts_file=seq_file,
            total_oligos=oligos_generated,
            params=params,
            decoded_success=decoded_success,
            grand_pass=g.grand_pass,
            grand_fail=g.grand_fail,
            chunks_recovered=sum(chunk_seen),
            total_chunks=len(chunk_seen),
            decode_time=decode_time,
            out_csv=f"coverage-analysis/seq-depth/files/final_results/{ecc_label}_a{alpha}.csv"
        )
    else:
        log_coverage_metrics(
            input_file=short_file,
            seq_counts_file=seq_file,
            total_oligos=oligos_generated,
            params=params,
            decoded_success=decoded_success,
            decode_time=decode_time,
            out_csv=f"coverage-analysis/seq-depth/files/final_results/{ecc_label}_a{alpha}.csv"
        )

    log_data_recovery_metrics(
        file_path=short_file,
        total_chunks=N,
        recovered_chunks=chunks_done,
        out_csv=f"coverage-analysis/seq-depth/files/data_recovery_{ecc_label}_a{alpha}.csv"
    )

    # print(f"Chunks seen: {chunk_seen}")
    # print(set(chunk_seen))
    print("Total chunks:", len(chunk_seen))
    print("Recovered chunks:", sum(chunk_seen))
    print(f"% Data Recovered: {sum(chunk_seen)/len(chunk_seen)*100}")

    plot_data_chunk_recovery(chunk_seen, f"coverage-analysis/seq-depth/files/data_recovery/chunk_recovery_{ecc_label}_a{alpha}.png")


if __name__ == "__main__":
    # for ecc in ["rs", "crc", "crc_grand"]:
    for i in np.arange(0.5, 10.5, 0.5):
        for _ in range(3):
            analyze_oligo_coverage(
                "coverage-analysis/seq-depth/files/lena.jpg",
                alpha=0.5,
                ecc_type="crc",
                subs_rate=0.0035,
                seq_depth=i
            )

    # analyze_oligo_coverage("coverage-analysis/seq-depth/files/lena.jpg", alpha=0.5, ecc_type="crc_grand", subs_rate=0.0035, seq_depth=0.5)