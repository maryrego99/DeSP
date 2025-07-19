import zlib
from reedsolo import RSCodec

def crc32_encoder(data:bytes) -> bytes:
    crc_val = zlib.crc32(data)
    return data + crc_val.to_bytes(4, byteorder='big')

def no_encoder(data:bytes) -> bytes:
    return data

def make_rs_encoder(rs_len:int):
    rs = RSCodec(rs_len)
    
    def rs_encoder(data:bytes) -> bytes:
        return rs.encode(data)
    
    return rs_encoder