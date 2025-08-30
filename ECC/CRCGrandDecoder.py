
import zlib
from .ECCDecoder import ECCDecoder
from Encode.Helper_Functions import dna_to_int_array
from DNAFountain.CRCHandler import *  

class CRCGrandDecoder(ECCDecoder):
    def __init__(self, max_flips=2, grand_variant='bitwise_brute', top_k=20, TM_matrix=None):
        self.max_flips = max_flips
        self.grand_variant = grand_variant
        self.top_k = top_k
        self.TM_matrix = TM_matrix
    

    def repair(self, dna_string, error_bit_position_counter):
        variant = self.grand_variant

        if variant.startswith("bitwise"):
            subtype = variant.replace("bitwise_","")
            if subtype == "heuristic_edges":
                return bitwise_heuristic_edges(dna_string, max_flips=self.max_flips,position_error_counter=error_bit_position_counter)
            elif subtype == "heuristic_topk":
                return bitwise_heuristic_topk(dna_string, max_flips=self.max_flips, top_k=self.top_k, position_error_counter=error_bit_position_counter)
            elif subtype == "heuristic_nok":
                return bitwise_heuristic_nok(dna_string, max_flips=self.max_flips, position_error_counter=error_bit_position_counter)
            else:
                return bitwise_brute(dna_string, max_flips=self.max_flips)
                

        
        elif variant.startswith("basewise"):
            subtype = variant.replace("basewise_","")
            if subtype == "heuristic_edges":
                return basewise_heuristic_edges(dna_string, max_flips=self.max_flips,TM_matrix=self.TM_matrix,position_error_counter=error_bit_position_counter)
            elif subtype == "heuristic_topk":
                return basewise_heuristic_topk(dna_string, max_flips=self.max_flips, TM_matrix=self.TM_matrix, top_k=self.top_k, position_error_counter=error_bit_position_counter)
            elif subtype == "heuristic_nok":
                return basewise_heuristic_nok(dna_string, max_flips=self.max_flips, TM_matrix=self.TM_matrix, position_error_counter=error_bit_position_counter)
            else:
                return basewise_brute(dna_string, max_flips=self.max_flips)