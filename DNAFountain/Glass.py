import zlib
import binascii
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import combinations
from Encode.Helper_Functions import *
from DNAFountain.Droplet import Droplet
from Encode.RPNG import *
from DNAFountain.CRCHandler import *
from ECC.ECCDecoder import ECCDecoder
from ECC.NoECCDecoder import NoECCDecoder
from ECC.RSDecoder import ReedSolomonDecoder
from ECC.CRCDecoder import CRCDecoder
from ECC.CRCGrandDecoder import CRCGrandDecoder
from ECC.ecc_encoders import crc32_encoder, make_rs_encoder, no_encoder
from Model.config import TM_NGS, TM_NNP
from collections import Counter
# from DNAFountain.CRCHandler import heuristic_grand_crc_repair

#----------------------------------------------------Glass-------------------------------------------------#        
class Glass:
    def __init__(self, in_file_name, chunk_num, header_size = 4, 
                 rs = 0, c_dist = 0.1, delta = 0.5, 
                flag_correct = True, gc = 0.05, max_homopolymer = 3, 
                max_hamming = 100, chunk_size = 32, exDNA = False, np = False, truth = None, ecc_decoder=None):
        
        self.entries = []
        self.droplets = set()
        self.num_chunks = chunk_num
        self.chunks = [None] * self.num_chunks
        self.header_size = header_size
        self.chunk_size = chunk_size
        self.exDNA = exDNA
        self.np = np
        self.chunk_to_droplets = defaultdict(set)
        self.done_segments = set()
        self.truth = truth
        self.in_file_name = in_file_name
        self.max_hamming = max_hamming

        self.PRNG = PRNG(K = self.num_chunks, delta = delta, c = c_dist, np = np)

        self.ecc_decoder = ecc_decoder or NoECCDecoder()

        # self.rs = rs
        # self.RSCodec = None
        self.correct = flag_correct
        self.seen_seeds = set()
        
        # if self.rs > 0:
        #     self.RSCodec = RSCodec(rs)
    
    def add_dna(self, dna_string):
        # transfer data to int
        data = dna_to_int_array(dna_string)
         
        if isinstance(self.ecc_decoder, CRCGrandDecoder):
            #GRAND uses CRC decoder itself for check
            crc_decoder = CRCDecoder()
            flag, data_corrected = crc_decoder.decode(data, original_dna=dna_string)    #first CRC check
        else:
            flag, data_corrected = self.ecc_decoder.decode(data)

        if flag == -1:
            return -1, None
        
        # split seed and payload
        seed_array = data_corrected[:self.header_size]
        seed = sum([   int(x)*256**i        for i, x in enumerate(seed_array[::-1])   ])
        payload = data_corrected[self.header_size:]
        self.add_seed(seed)
        
        # decode
        self.PRNG.set_seed(seed)
        blockseed, d, ix_samples = self.PRNG.get_src_blocks_wrap() #reconstruct the linear combination
        # print(blockseed,d,ix_samples,payload)

        #create droplet based on ecc used
        if isinstance(self.ecc_decoder, CRCGrandDecoder) or isinstance(self.ecc_decoder, CRCDecoder):
            ecc_encoder=crc32_encoder
        if isinstance(self.ecc_decoder, ReedSolomonDecoder):
            ecc_encoder=make_rs_encoder(rs_len=4)
        if isinstance(self.ecc_decoder, NoECCDecoder):
            ecc_encoder=no_encoder
        d = Droplet(payload, seed, ix_samples, ecc_encoder=ecc_encoder) #reconstruct droplet
        self.addDroplet(d) 
        return seed, data

    def addDroplet(self, droplet):
        self.droplets.add(droplet)
        for chunk_num in droplet.num_chunks:
            self.chunk_to_droplets[chunk_num].add(droplet) #we document for each chunk all connected droplets        
        self.updateEntry(droplet) #one round of message passing

        
    def updateEntry(self, droplet):
        #removing solved segments from droplets
        for chunk_num in (droplet.num_chunks & self.done_segments):
            droplet.data = xor(droplet.data,self.chunks[chunk_num]) #subtract (ie. xor) the value of the solved segment from the droplet.
            droplet.num_chunks.remove(chunk_num)
            self.chunk_to_droplets[chunk_num].discard(droplet) #remove the edge between droplet and input segment.

        #solving segments when the droplet have exactly 1 segment
        if len(droplet.num_chunks) == 1: #the droplet has only one input segment
            lone_chunk = droplet.num_chunks.pop()
            self.chunks[lone_chunk] = droplet.data #assign the droplet value to the input segment (=entry[0][0])
            self.done_segments.add(lone_chunk) #add the lone_chunk to a data structure of done segments.
            self.droplets.discard(droplet) #remove the edge between the droplet and input segment
            self.chunk_to_droplets[lone_chunk].discard(droplet) #remove the edge between the input segment and the droplet
            #update other linked droplets
            for other_droplet in self.chunk_to_droplets[lone_chunk].copy():
                self.updateEntry(other_droplet)

    def add_seed(self, seed):
        self.seen_seeds.add(seed)

    def len_seen_seed(self):
        return len(self.seen_seeds)

    def isDone(self):
        if self.num_chunks - len(self.done_segments) > 0:
            return None 
        return True

    def chunksDone(self):
        return len(self.done_segments)

    def String(self):
        res = ''
        for x in self.chunks:
            res += ''.join(map(chr, x))
        return res
    
    def StringNoPadding(self):
        return self.String().rstrip('\0')

    def removePadding(self,pad):
        if pad != -1:
            self.chunks[-1] = self.chunks[-1][:-pad]

        crp = []
        for b in self.chunks[-1]:
            if 0 == b:
                break 
            crp.append(b)
        self.chunks[-1] = crp
        return crp
    
    def check_image_header(self, binary_data):
        if binary_data.startswith(b'\xFF\xD8'):
            return 'jpg'
        elif binary_data.startswith(b'\x89PNG\r\n\x1a\n'):
            return 'png'
        return None

    def save_partial_if_valid(self, file_name_base, chunks, pad=-1):
        if pad != -1 and pad > 0:
            chunks = chunks[:-pad]

        binary_data = b''.join(bytes(c) for c in chunks if c is not None)
        img_type = self.check_image_header(binary_data)

        if img_type:
            file_name = f"{file_name_base}.{img_type}"
            with open(file_name, 'wb') as f:
                f.write(binary_data)
            print(f"Partial image with valid header saved to: {file_name}")
            with open("IO/Input/lena.jpg", "rb") as f:
                original_header = f.read(1024)

            patched = original_header + binary_data[1024:]

            with open("IO/Output/patched_partial.jpg", "wb") as f:
                f.write(patched)
        else:
            print("Partial output does not contain a valid image header — not saved.")
            print("Decoding failed.")
    
    def save(self,file_name, pad = -1):
        self.removePadding(pad)
        with open(file_name,'wb') as f:
            for c in self.chunks:
                f.write(bytes(c))
#             logging.info('saved')
            print('Decoding Successful. Image saved')
            f.close()
        
    def binString(self):
        bs = b''
        for c in self.chunks:
            bs += bytes(c)
        return bs
    
    def bchunks(self):
        chunks = []
        for c in self.chunks:
            chunks.append(bytes(c))
        return chunks
    
    def print_chunks(self):
        print(self.chunks)
        
    def display_chunks(self):
        i = 0
        not_none = []
        for x in self.chunks:
            print(i,''.join(map(chr, x)))
            i+=1
            if x!= None:
                not_none.append(i)
        return not_none


    def log_error_profile(self, counter, label, output_path):
        if counter is None or not counter:
            return

        print(f"\nTop {label} positions flipped in successful GRAND repairs:")
        print(counter.most_common(10))

        positions = sorted(counter.keys())
        frequencies = [counter[pos] for pos in positions]

        plt.figure(figsize=(12, 5))
        plt.bar(positions, frequencies)
        plt.xlabel(f"{label.capitalize()} Position", fontsize=16)
        plt.ylabel("Frequency in Successful GRAND Repairs", fontsize=16)
        # plt.title(f"Position-Based Error Profile ({label.capitalize()} Level)")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()
        plt.savefig(output_path)
        plt.close()
        print("error-profile saved")
    
    
    def finalize_decoding(self, line, solve_num, errors, coverage_vs_reads, chunk_seen,
                       crc_pass, crc_fail, grand_pass, grand_fail,
                       repaired_strands, attempted_grand,error_bit_position_counter=None, error_base_position_counter=None):
        
        if error_bit_position_counter is not None:
            self.log_error_profile(
                counter=error_bit_position_counter,
                label="bit",
                output_path="coverage-analysis/visualizations/error-profile/bit-error-profile.pdf"
            )

        if error_base_position_counter is not None:
            self.log_error_profile(
            counter=error_base_position_counter,
            label="base",
            output_path="coverage-analysis/visualizations/error-profile/base-error-profile.pdf"
        )

        print(f"Originally CRC Pass: {crc_pass}, CRC Fail: {crc_fail}, Total Reads from synthesis: {line}")
        # usable_ratio = crc_pass / (crc_pass + crc_fail)
        # print(f"Usable droplet ratio: {usable_ratio:.2%}")
        if isinstance(self.ecc_decoder, CRCGrandDecoder):
            print(f"Attempted GRAND on strands: {len(attempted_grand)}")
            print(f"Repaired strands using GRAND: {len(repaired_strands)}")
            print(f"GRAND Pass: {grand_pass}, GRAND Fail: {grand_fail}")
            print(f"GRAND Success Rate: {grand_pass/(grand_pass+grand_fail)*100:.2f}")
        # print(repaired_strands)
        print(f"Valid Droplets: {len(self.valid_droplets)}")
        self.recovered_droplets = len(self.valid_droplets)

        status_code = 0 if self.isDone() else -1
        return status_code, solve_num, line, self.chunksDone(), errors, coverage_vs_reads, chunk_seen, self.chunks

    def decode(self):
        f = open(self.in_file_name, 'r')
        line = 0
        errors = 0
        solve_num = []
        crc_pass = 0
        crc_fail = 0
        grand_pass = 0
        grand_fail = 0
        self.valid_droplets = []
        repaired_strands = []
        attempted_grand = []

        chunk_seen = [0] * self.num_chunks
        coverage_vs_reads = []

        error_bit_position_counter = Counter()
        error_base_position_counter = Counter()

        while True:
            try:
                dna = f.readline().rstrip('\n')
            except:
                return self.finalize_decoding(
                    line, solve_num, errors, coverage_vs_reads, chunk_seen,
                    crc_pass, crc_fail, grand_pass, grand_fail,
                    repaired_strands, attempted_grand, error_bit_position_counter, error_base_position_counter
                )

            if len(dna) == 0:
                return self.finalize_decoding(
                    line, solve_num, errors, coverage_vs_reads, chunk_seen,
                    crc_pass, crc_fail, grand_pass, grand_fail,
                    repaired_strands, attempted_grand, error_bit_position_counter, error_base_position_counter
                )

            line += 1
            seed, data = self.add_dna(dna) # first error check - RS or CRC - based on this seed is set

            if isinstance(self.ecc_decoder, ReedSolomonDecoder):
                if seed == -1:
                    errors += 1

            if isinstance(self.ecc_decoder, CRCDecoder):
                if seed == -1:
                    errors += 1
                    crc_fail += 1  # CRC failed
                else:
                    crc_pass += 1  # CRC passed
                    self.valid_droplets.append((seed, data))

            if isinstance(self.ecc_decoder, CRCGrandDecoder):
                if seed == -1:
                    # first CRC check failed
                    crc_fail += 1
                    attempted_grand.append((seed, data))

                    repaired_dna = self.ecc_decoder.repair(dna, error_bit_position_counter)
                    
                    if repaired_dna:
                        # repaired_strands.append(repaired_dna)
                        seed, data = self.add_dna(repaired_dna) # second crc check after repair
                        if seed == -1:
                            # if second failed too (after repair)
                            errors += 1
                            grand_fail += 1
                        else:
                            grand_pass += 1
                            self.valid_droplets.append((seed, data))
                            repaired_strands.append(repaired_dna)
                    else:
                        errors += 1
                        grand_fail += 1
                else:
                    crc_pass += 1
                    self.valid_droplets.append((seed, data))

            self.crc_pass = crc_pass
            self.crc_fail = crc_fail

            self.grand_pass = grand_pass
            self.grand_fail = grand_fail
            

            if line % 200 == 0:
                pass

            # if line == 1:
            #     chunk_seen = [0] * self.num_chunks
            #     coverage_vs_reads = []

            if seed != -1:
                self.PRNG.set_seed(seed)
                # blockseed, d, ix_samples = self.PRNG.get_src_blocks_wrap()
                # for chunk_id in ix_samples:
                #     chunk_seen[chunk_id] = 1
                chunk_seen = [1 if i in self.done_segments else 0 for i in range(self.num_chunks)]
                coverage_vs_reads.append(sum(chunk_seen))
            solve_num.append(self.chunksDone())

            if self.isDone():
                return self.finalize_decoding(
                    line, solve_num, errors, coverage_vs_reads, chunk_seen,
                    crc_pass, crc_fail, grand_pass, grand_fail,
                    repaired_strands, attempted_grand, error_bit_position_counter, error_base_position_counter
                )


