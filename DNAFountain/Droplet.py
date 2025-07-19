import zlib
from Encode.Helper_Functions import byte_to_dna

#----------------------------------------------------Droplet-------------------------------------------------#  
class Droplet:
    def __init__(self, data, seed, num_chunks = None, rs = 0, rs_obj = None, degree = None, ecc_encoder=None):
        #num_chunks is a list of the orignal packets numbers used to xor 
        #rs is the number of Reed Solomon symbols to add to the message

        self.data = data
        self.seed = seed
        self.num_chunks = set(num_chunks)
        # self.rs = rs
        # self.rs_obj = rs_obj
        self.degree = degree

        self.DNA = None

        self.ecc_encoder = ecc_encoder
    
    def toDNA(self, flag = None):
        #this function wraps the seed, data payload, and Reed Solomon (or CRC).
        if self.DNA is not None:
            return self.DNA
        self.DNA = byte_to_dna(self._package())
        return self.DNA
    
    def chunkStr(self):
        num = 0
        s = ''
        for i in self.num_chunks:
            if(6 == num):
                s += '...'
                break
            s+= str(i) + ' '
        num += 1
        return s
        
    def _package(self):
        #this function converts the seed to a list of 4bytes HARD CODED!!!
        #adds the seed to the data (list of integers)
        #computes a reed solomon on the seed+data.
        #returns everything.

        seed_ord = self.seed.to_bytes(4, byteorder = 'big')
        message = seed_ord + bytes(self.data)

        if self.ecc_encoder is not None:
            message = self.ecc_encoder(message)

        return message