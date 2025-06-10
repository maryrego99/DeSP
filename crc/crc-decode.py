import zlib
import random

class FakeDroplet:
    def __init__(self, payload: bytes, seed: int):
        self.seed = seed
        self.payload = payload
        self.crc = None
        self.dna = None

    def add_crc32(self):
        self.crc = zlib.crc32(self.payload)
        return self.payload + self.crc.to_bytes(4, 'big')

    def to_dna(self):
        data_with_crc = self.add_crc32()
        self.dna = self.byte_to_fake_dna(data_with_crc)
        return self.dna

    def corrupt(self, error_rate=0.01):
        data = bytearray(self.payload + self.crc.to_bytes(4, 'big'))
        for i in range(len(data)):
            if random.random() < error_rate:
                data[i] ^= 0xFF  # Flip bits
        return bytes(data)

    def verify_crc(self, received):
        data_part = received[:-4]
        crc_part = int.from_bytes(received[-4:], 'big')
        calc_crc = zlib.crc32(data_part)
        return calc_crc == crc_part

    def byte_to_fake_dna(self, b):
        dna_map = ['A', 'C', 'G', 'T']
        return ''.join(dna_map[byte % 4] for byte in b)

### Simulation Parameters
n_droplets = 100
error_rate = 0.1

# Create N fake droplets
original_data_list = [f"chunk{i}".encode() for i in range(n_droplets)]
droplets = [FakeDroplet(data, seed=i) for i, data in enumerate(original_data_list)]

valid_reads = []
failed_reads = []

for droplet in droplets:
    dna = droplet.to_dna()
    received = droplet.corrupt(error_rate=error_rate)
    if droplet.verify_crc(received):
        valid_reads.append(droplet.payload)
    else:
        failed_reads.append(droplet.payload)

### Decoding summary
print(f"Total droplets: {n_droplets}")
print(f"Valid reads (passed CRC): {len(valid_reads)}")
print(f"Failed reads (CRC failed): {len(failed_reads)}")

# Simulated decoding: unique payloads recovered
recovered = set(valid_reads)
print(f"Unique payloads recovered after filtering: {len(recovered)}")


max_reads = 300
n_droplets = 100
error_rate = 0.1
seen_payloads = set()
reads_attempted = 0

while len(seen_payloads) < n_droplets and reads_attempted < max_reads:
    idx = random.randint(0, n_droplets - 1)
    d = FakeDroplet(f"chunk{idx}".encode(), seed=idx)
    d.to_dna()
    received = d.corrupt(error_rate)
    reads_attempted += 1

    if d.verify_crc(received):
        seen_payloads.add(d.payload)

print(f"Reached full coverage of {n_droplets} droplets using {reads_attempted} valid reads")