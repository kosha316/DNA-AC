import psutil
from Crypto import Random
from Crypto.Cipher import AES
from Crypto.Random import get_random_bytes
from Crypto.PublicKey import RSA
from Crypto.Cipher import PKCS1_v1_5
from Crypto.Cipher import ARC4
from ecdsa import SigningKey, VerifyingKey, NIST256p
from PIL import Image
import tenseal as ts
import numpy as np
import hashlib
import hmac
from phe import paillier
import random
import os
import trans_type as trans
import seq_methods
import coder
import process as ps
import time
import securitytest as st


def add_to_size(byte_s):
    while len(byte_s) % 16 != 0:
        byte_s += ' '.encode()
    # print(len(byte_s))
    return byte_s


def remove_blank(byte_s):
    return byte_s.strip()


def encry_aes(plain_path, key_bytes, cipher_path, iv, seg_length):
    with open(plain_path, "rb") as f:
        plain_bytes = f.read()
        plain_bytes = coder.ran1bit(plain_bytes)
        plain_bytes = add_to_size(plain_bytes)
        aes_object = AES.new(key_bytes, AES.MODE_CBC, iv)
        cipher = aes_object.encrypt(plain_bytes)
        f.close()

    cipher_list = trans.bytes_to_fasta(cipher, seg_length)
    cipher_list = seq_methods.add_index(cipher_list)
    print(len(cipher_list))
    cipher_list = ps.homo_encode(cipher_list, 5)
    cipher_list = coder.seq_pad(cipher_list, seg_length + index_length, 'AGCT')
    seq_methods.write_seq(cipher_path, cipher_list)

    return 0


def decry_aes(decry_path, iv, key_bytes, cipher_path):
    cipher_list = seq_methods.read_seq(cipher_path)
    cipher_list = coder.remove_pad(cipher_list, 'AGCT')
    cipher_list = ps.homo_decode(cipher_list, 5)
    list_len = len(cipher_list)
    cipher_list = seq_methods.remove_index(cipher_list, list_len)
    cipher_bytes = trans.fasta_to_bytes(cipher_list)

    decry = AES.new(key_bytes, AES.MODE_CBC, iv)
    decry_bytes = decry.decrypt(cipher_bytes)
    decry_bytes = remove_blank(decry_bytes)

    with open(decry_path, "wb") as f:
        f.write(decry_bytes)
        f.close()


def arc4_encry(plain_path, key_byte, cipher_path, seg_length):
    with open(plain_path, "rb") as f:
        plain_bytes = f.read()
        # 雪崩测试
        plain_bytes = coder.ran1bit(plain_bytes)
        arc4_object = ARC4.new(key_byte)
        cipher = arc4_object.encrypt(plain_bytes)
        f.close()

    cipher_list = trans.bytes_to_fasta(cipher, seg_length)
    cipher_list = seq_methods.add_index(cipher_list)
    print(len(cipher_list))
    cipher_list = ps.homo_encode(cipher_list, 5)
    cipher_list = coder.seq_pad(cipher_list, seg_length + index_length, 'AGCT')
    seq_methods.write_seq(cipher_path, cipher_list)

    return 0


def arc4_decry(decry_path, key_byte, cipher_path):
    cipher_list = seq_methods.read_seq(cipher_path)
    cipher_list = coder.remove_pad(cipher_list, 'AGCT')
    cipher_list = ps.homo_decode(cipher_list, 5)
    list_len = len(cipher_list)
    cipher_list = seq_methods.remove_index(cipher_list, list_len)
    cipher_bytes = trans.fasta_to_bytes(cipher_list)

    arc4_object = ARC4.new(key_byte)
    decry_bytes = arc4_object.decrypt(cipher_bytes)

    with open(decry_path, "wb") as f:
        f.write(decry_bytes)
        f.close()
    return 0

def encry_rsa(plain_path, key_byte, cipher_path):
    with open(plain_path, "rb") as f:
        plain_bytes = f.read()
        plain_bytes = add_to_size(plain_bytes)
        f.close()
    rsa_key = RSA.importKey(key_byte)
    cipher = PKCS1_v1_5.new(rsa_key)
    cipher_bytes = bytes()
    for i in range(0, len(plain_bytes), 200):
        cipher_bytes += (cipher.encrypt(plain_bytes[i:i + 200]))
    with open(cipher_path, "wb") as f:
        f.write(cipher_bytes)
        f.close()


def decry_rsa(decry_path, key_byte, cipher_path):
    with open(cipher_path, "rb") as f:
        cipher_bytes = f.read()
        f.close()
    rsa_key = RSA.importKey(key_byte)
    plain = PKCS1_v1_5.new(rsa_key)
    decry_bytes = bytes()
    for i in range(0, len(cipher_bytes), 256):
        decry_bytes += (plain.decrypt(cipher_bytes[i:i + 256], "error"))
    with open(decry_path, "wb") as f:
        f.write(decry_bytes)
        f.close()


def add_hmac(plain_path, key_byte, cipher_path):
    with open(plain_path, "rb") as f:
        plain_bytes = f.read()
        hmac_object = hmac.new(key_byte, plain_bytes, digestmod=hashlib.sha256)
        hmac_value = hmac_object.digest()
        f.close()
    with open(cipher_path, "wb") as f:
        f.write(hmac_value)
        f.write(plain_bytes)
        f.close()
    return 0


def verify_hmac(plain_path, key_byte):
    with open(plain_path, "rb") as f:
        ori_hmac_value = f.read(32)
        plain_bytes = f.read()
        f.close()
    hmac_object = hmac.new(key_byte, plain_bytes, digestmod=hashlib.sha256)
    plain_hmac_value = hmac_object.digest()
    # print(ori_hmac_value)
    # print(plain_hmac_value)
    return hmac.compare_digest(ori_hmac_value, plain_hmac_value)


def ecc_encry(plain_path, key_byte, cipher_path):
    with open(plain_path, "rb") as f:
        plain_bytes = f.read()
        hash_obj = hashlib.sha256(plain_bytes).digest()
        ecc_sign = key_byte.sign(hash_obj)
        f.close()
    with open(cipher_path, "wb") as f:
        f.write(ecc_sign)
        f.write(plain_bytes)
        f.close()
    return 0


def ecc_verify(plain_path, key_byte):
    with open(plain_path, "rb") as f:
        ecc_sign = f.read(64)
        plain_bytes = f.read()
        f.close()
    hash_obj = hashlib.sha256(plain_bytes).digest()
    ver_key = key_byte.get_verifying_key()
    valid = ver_key.verify(ecc_sign, hash_obj)
    return valid


def he_encry(plain_path, context, cipher_path):
    # Read image into Python
    img = Image.open(plain_path)
    # Convert image to matrix
    pixels = list(img.getdata())
    matrix = np.array(pixels)
    cipher1 = ts.ckks_tensor(context, matrix)
    np.save(cipher_path, cipher1)
    return 0


def paillier_encry(plain_path, key_pair, cipher_path):
    public_key = key_pair[1]
    with open(plain_path, "rb") as f:
        plain_bytes = f.read()
        plain_bytes = bytearray(plain_bytes)
        encrypted_number_list = [public_key.encrypt(x) for x in plain_bytes]
        f.close()
    with open(cipher_path, "wb") as f:
        f.write(bytes(encrypted_number_list))
        f.close()
    return 0


def get_context():
    context = ts.context(ts.SCHEME_TYPE.CKKS, 8192, coeff_mod_bit_sizes=[22, 21, 21, 21, 21, 21, 21, 21, 21, 21])
    context.global_scale = pow(2, 21)
    context.generate_galois_keys()
    return context


def get_key_pair():
    public_key, private_key = paillier.generate_paillier_keypair()
    return [private_key, public_key]


def rsa_generate_key():
    rsa = RSA.generate(2048, Random.new().read)
    pri_key = rsa.exportKey()
    pub_key = rsa.public_key().exportKey()
    return [pri_key, pub_key]


def ecc_generate_key():
    pri_key = SigningKey.generate(curve=NIST256p)
    return pri_key


if __name__ == "__main__":
    current_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "secure")

    keys_path = os.path.join(current_path, "aes", "keys.txt")
    context_path = os.path.join(current_path, "aes", "context.txt")
    plain_path = os.path.join(current_path, "input", "x_ray.png")
    recon_path = os.path.join(current_path, "aes", "output-results.txt")
    cipher_fasta_path = os.path.join(current_path, "aes", "x_ray.txt")
    decry_path = os.path.join(current_path, "aes", "x_ray_decry.png")
    decry_key_path = os.path.join(current_path, "aes", "decry_key.txt")
    err_cipher_path = os.path.join(current_path, "aes", "x_ray_ks.txt")
    code_mode = AES.MODE_CBC
    #

    #
    segment_length = 128
    index_length = 8
    #
    # encry_aes(plain_path, key, cipher_path)
    # decry_aes(decry_path, key, cipher_path)
    #
    # rsa_key = rsa_generate_key()
    # encry_rsa(plain_path, keys[1], cipher_path)
    # decry_rsa(decry_fasta_path, keys[0], cipher_fasta_path, 64)
    #
    # add_hmac(plain_path, key, cipher_path)
    # print(verify_hmac(cipher_path, key))
    #
    # key_bytes = seq_methods.arc4_keys_write(keys_path)
    # key_bytes = seq_methods.arc4_keys_read(keys_path)
    # iv, key_bytes = seq_methods.aes_keys_write(keys_path)
    iv, key_bytes = seq_methods.aes_keys_read(keys_path)


    # process = psutil.Process(os.getpid())
    # start_time = time.time()
    # start_memory = process.memory_info().rss
    encry_aes(plain_path, key_bytes, err_cipher_path, iv, segment_length)
    # decry_aes(decry_path, iv, key_bytes, recon_path)

    # arc4_encry(plain_path, key_bytes, err_cipher_path, segment_length)
    # arc4_decry(decry_path, key_bytes, cipher_fasta_path)
    # end_memory = process.memory_info().rss
    # end_time = time.time()
    # print(f"内存使用：{round((end_memory - start_memory) / 1024 / 1024, 4)}MB")
    # print(f"运行时间：{round(end_time - start_time, 4)}S")


    # ecc_encry(plain_path, ecc_key, cipher_path)
    # print(ecc_verify(cipher_path, ecc_key))
    # paillier_encry(plain_path, key_pair, cipher_path)
    # st.caculate_difference_file(cipher_fasta_path, err_cipher_path)
    st.caculate_nbcr(cipher_fasta_path, err_cipher_path)
    st.caculate_baci(cipher_fasta_path, err_cipher_path)
    # print("file content is same?    ", trans.compare_files(plain_path, decry_path))
