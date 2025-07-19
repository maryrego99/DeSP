from .ECCDecoder import ECCDecoder

class NoECCDecoder(ECCDecoder):
    def decode(self, data):
        return 0, data      # No correction