from reedsolo import RSCodec
from .ECCDecoder import ECCDecoder
from Encode.Helper_Functions import rs_decode

class ReedSolomonDecoder(ECCDecoder):
    def __init__(self, rs_len):
        self.rs_len = rs_len
        self.rs_codec = RSCodec(rs_len)
    
    def decode(self, data, original_dna=None):
        flag, data_corrected = rs_decode(data, self.rs_codec)
        if flag == -1:
            return -1, None
        
        return 0, data_corrected