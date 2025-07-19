import logging
from reedsolo import RSCodec
from Encode.Helper_Functions import *
from Encode.RPNG import *
from DNAFountain.Droplet import Droplet
    
#----------------------------------------------------Fountain-------------------------------------------------#      
class DNAFountain:

    def __init__(self, 
                file_in,
                alpha, 
                stop = None,
                rs = 0, 
                c_dist = 0.1, 
                delta = 0.5, 
                scanner = None,
                ecc_encoder=None,
                ecc_name="crc32",
                ecc_bytes=4
                ):

        #alpha is the redundency level
        #stop is whether we have a limit on the number of oligos
        #chunk_size and file_size are in bytes
        #rs is the number of bytes for reed-solomon error correcting code over gf(2^8).
        #c_dist is a parameter of the degree distribution
        #delta is a parameter of the degree distribution
        #np: should we use numpy random number generator? Faster, but incompatible for previous versions
        #max_homopolymer: the largest homopolymer allowed
        #gc: the allowable range of gc +- 50%

        #data:
        self.file_in = file_in
        self.chunk_size = len(file_in[0])
        self.num_chunks = len(file_in)
        
        #reduancy:
        self.alpha = alpha
        self.stop = stop
        self.final = self.calc_stop()

        #random mnumber generator
        self.lfsr = lfsr(lfsr32s(), lfsr32p()) #starting an lfsr with a certain state and a polynomial for 32bits.
        self.lfsr_l = len('{0:b}'.format( lfsr32p())) - 1 #calculate the length of lsfr in bits 
        self.seed = self.lfsr.__next__()

        self.PRNG = PRNG(K = self.num_chunks, delta = delta, c = c_dist, np = False) #creating the solition distribution object
        self.PRNG.set_seed(self.seed)

        #error correcting code:
        # self.rs = rs #the number of symbols (bytes) to add
        # self.rs_obj = RSCodec(self.rs) #initalizing an reed solomon object

        #biological screens:
        self.scanner = scanner
        if self.scanner == None:
            self.scanner = Scanner()

        self.ecc_encoder = ecc_encoder
        self.ecc_name = ecc_name
        self.ecc_bytes = ecc_bytes # to append

        self.tries = 0 #number of times we tried to create a droplet
        self.good = 0 #droplets that were screened successfully.
        
        self.oligo_l = self.calc_oligo_length()
        # store the generated droplets
        self.dna_df = None
        self.dna_dl = []
    
    
    def calc_oligo_length(self):
        #return the number of nucleotides in an oligo:
        bits = self.chunk_size * 8 + self.lfsr_l + (self.ecc_bytes * 8)  # CRC-32 = 4 bytes = 32 bits
        return bits / 4


    def calc_stop(self):
        if self.stop is not None:
            return self.stop
        stop = int(self.num_chunks*(1+self.alpha))+1
        return stop

    def droplet(self):
        #creating a droplet.
        data = None

        d, num_chunks = self.rand_chunk_nums() #creating a random list of segments.

        # print(num_chunks)
        for num in num_chunks: #iterating over each segment
            if data is None: #first round. data payload is empty.
                data = self.chunk(num) #just copy the segment to the payload.
            else: #more rounds. Xor the new segments with the payload.
                data = xor(data,self.chunk(num))  #map(operator.xor, data, self.chunk(num))
        
        # print(data)
        self.tries +=  1 

        #we have a droplet:
        return Droplet(data = data, 
                       seed = self.seed,
                       num_chunks = num_chunks,
                       degree = d,
                       ecc_encoder=self.ecc_encoder)

    def chunk(self, num):
        #return the num-th segment from the file
        return self.file_in[num]

    #-------------------generate random chunk numebers----------------# 
    def updateSeed(self):
        #This function creates a fresh seed for the droplet and primes the solition inverse cdf sampler
        self.seed = self.lfsr.__next__() #deploy one round of lfsr, and read the register.
        self.PRNG.set_seed(self.seed) #update the seed with the register

    def rand_chunk_nums(self):
        #This funcation returns a subset of segments based on the solition distribution.
        #It updates the lfsr to generates a new seed.
        self.updateSeed() #get a fresh seed and prime the solition inverse cdf sampler.
        blockseed, d, ix_samples = self.PRNG.get_src_blocks_wrap()
        return d, ix_samples #return a list of segments.

    #----------------screen generated droplets----------------------#
    def screen(self, droplet):
        if self.scanner.Pass(droplet.toDNA()):
            self.good += 1
            dna = droplet.toDNA()
            degree = droplet.degree
            chunk_str = droplet.chunkStr()
            seed = droplet.seed
            self.dna_dl.append([dna,seed,degree,chunk_str])
            return 1
        return 0

    def save(self,file_name = 'out.dna'):
        with open(file_name, 'w') as f:
            # f.write('Fountain code\n')
            # f.write('CN: ' + str(self.num_chunks) +'\n')
            # f.write('CL: ' + str(self.chunk_size) + '\n')
            # f.write('RS: ' + str(self.rs) + '\n')
            f.writelines('\n'.join([d[0] for d in self.dna_dl]))
            f.close()

    def encode(self):
        self.dl = []
        self.tries = 0
        self.good = 0
        while self.good < self.final:
            self.screen(self.droplet())
            if self.tries%2000 == 0:
                logging.info("generate %d chunks after %d tries",self.good, self.tries)
                # print("generate %d chunks after %d tries" % (self.good, self.tries))
                
        logging.info("Finish generating %d chunks after %d tries", self.good,self.tries)
        # print("Finish generating %d chunks after %d tries"% (self.good,self.tries))
        return self.good, self.tries    