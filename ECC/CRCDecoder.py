import zlib
from .ECCDecoder import ECCDecoder

class CRCDecoder(ECCDecoder):
    def decode(self, data, original_dna=None):
        if len(data) < 4:
            return -1, None

        received_crc = int.from_bytes(data[-4:], byteorder='big')
        data_wo_crc = data[:-4]
        calc_crc = zlib.crc32(bytes(data_wo_crc))

        if calc_crc != received_crc:
            return -1, None

        return 0, data_wo_crc