from Crypto.Cipher import AES
from Crypto.Random import get_random_bytes
import random
import os
import secure.trans_type as trans
# from Crypto.Util.Padding import pad
import operator
from Crypto import Random
from Crypto.PublicKey import RSA
from Crypto.Cipher import PKCS1_v1_5
from Crypto.Cipher import ARC4
from ecdsa import SigningKey, VerifyingKey, NIST256p
import hmac
import hashlib
import securitytest as st
import json
import process as ps
from Homopolymer import encoder, decoder
import itertools
import math
import editdistance
import time
import psutil
import biotest
import coder

dic1 = ['A', 'T', 'C', 'G']
exor_table = [
    ['A', 'T', 'C', 'G'],
    ['T', 'A', 'G', 'C'],
    ['C', 'G', 'A', 'T'],
    ['G', 'C', 'T', 'A']
]


def read_seq(file_path):
    with open(file_path, "r", encoding='utf-8') as f:
        fasta_segments = f.read()
        fasta_segments = fasta_segments.strip('\n')
        seq_list = fasta_segments.split('\n')
        f.close()
        return seq_list


def write_seq(file_path, seq_list):
    with open(file_path, 'w', encoding='utf-8') as f:
        for item in seq_list:
            f.write(item)
            f.write('\n')
    f.close()


def encry_aes(plain_path, IV, key_byte, cipher_path, seg_length, index_length, mode_type):
    plain_list = coder.encode(plain_path, seg_length, index_length, pad_str='AGCT')
    # print(plain_list)
    plain_bytes = trans.fasta_to_bytes(plain_list)
    # plain_bytes = ran1bit(plain_bytes)
    size_pad_set = size_pad(plain_bytes)
    plain_bytes = size_pad_set['data']
    size_pad_num = size_pad_set['pad_num']
    cipher = AES.new(key_byte, mode_type, IV)
    cipher_bytes = cipher.encrypt(plain_bytes)
    cipher_bytes_size = len(cipher_bytes)
    seq_pad_set = seq_pad(cipher_bytes, seg_length + index_length)
    cipher_bytes = seq_pad_set['data']
    # seq_pad_num = seq_pad_set['pad_num']

    cipher_list = trans.bytes_to_fasta(cipher_bytes, seg_length + index_length)

    # print(len(cipher_list))
    cipher_list = add_index(cipher_list)
    # print(len(cipher_list))
    cipher_list = ps.homo_encode(cipher_list, 5)
    cipher_list = coder.seq_pad(cipher_list, seg_length + index_length, 'AGCT')
    write_seq(cipher_path, cipher_list)
    # ps.find_homo(cipher_list, 4)
    return size_pad_num, cipher_bytes_size


def decry_aes(plain_path, IV, key_byte, cipher_path, seg_length, index_length, mode_type, size_pad_num,
              cipher_bytes_size):
    cipher_list = read_seq(cipher_path)
    cipher_list = coder.remove_pad(cipher_list, 'AGCT')
    cipher_list = ps.homo_decode(cipher_list, 5)
    list_len = len(cipher_list)
    cipher_list = remove_index(cipher_list, list_len)

    cipher_bytes = trans.fasta_to_bytes(cipher_list)
    cipher_bytes = satisfy_bound(cipher_bytes, cipher_bytes_size)
    plain = AES.new(key_byte, mode_type, IV)
    plain_bytes = plain.decrypt(cipher_bytes)
    plain_bytes = remove_pad(plain_bytes, size_pad_num)
    plain_list = trans.bytes_to_fasta(plain_bytes, seg_length + index_length)
    # print(plain_list)

    # print(len(plain_list))
    coder.decode(plain_list, index_length, 'AGCT', plain_path)


def arc4_encry(plain_path, key_byte, cipher_path, seg_length, index_length):
    plain_list = coder.encode(plain_path, seg_length, index_length, pad_str='AGCT')
    plain_bytes = trans.fasta_to_bytes(plain_list)

    size_pad_set = size_pad(plain_bytes)
    plain_bytes = size_pad_set['data']
    size_pad_num = size_pad_set['pad_num']
    arc4_object = ARC4.new(key_byte)
    cipher_bytes = arc4_object.encrypt(plain_bytes)
    cipher_bytes_size = len(cipher_bytes)
    seq_pad_set = seq_pad(cipher_bytes, seg_length + index_length)
    cipher_bytes = seq_pad_set['data']

    cipher_list = trans.bytes_to_fasta(cipher_bytes, seg_length + index_length)
    cipher_list = add_index(cipher_list)
    print(len(cipher_list))
    cipher_list = ps.homo_encode(cipher_list, 5)
    cipher_list = coder.seq_pad(cipher_list, seg_length + index_length, 'AGCT')
    write_seq(cipher_path, cipher_list)
    # ps.find_homo(cipher_list, 4)
    return size_pad_num, cipher_bytes_size


def arc4_decry(plain_path, key_byte, cipher_path, seg_length, size_pad_num, cipher_bytes_size, index_length):
    cipher_list = read_seq(cipher_path)
    cipher_list = coder.remove_pad(cipher_list, 'AGCT')
    cipher_list = ps.homo_decode(cipher_list, 5)
    list_len = len(cipher_list)
    cipher_list = remove_index(cipher_list, list_len)
    cipher_bytes = trans.fasta_to_bytes(cipher_list)
    cipher_bytes = satisfy_bound(cipher_bytes, cipher_bytes_size)
    arc4_object = ARC4.new(key_byte)
    plain_bytes = arc4_object.decrypt(cipher_bytes)
    plain_bytes = remove_pad(plain_bytes, size_pad_num)
    plain_list = trans.bytes_to_fasta(plain_bytes, seg_length + index_length)
    coder.decode(plain_list, index_length, 'AGCT', plain_path)


def exa_encry(key_path, plain_path, cipher_path, seq_length, decry_key_path, index_length):
    plain_list = coder.encode(plain_path, seq_length, index_length, pad_str='AGCT')
    key_len = len(plain_list) * (seq_length + index_length)
    with open(key_path, 'rb') as kf:
        fullkey = kf.read(int(key_len / 8))

    fullkey = bytes_to_binary_string(fullkey)
    key_list = string_to_sublists(fullkey, seq_length + index_length)

    decry_key, cipher_list = exor(key_list, plain_list, 5)
    cipher_list = exa_index(cipher_list, 5)
    print(len(cipher_list))
    write_seq(cipher_path, cipher_list)
    write_seq(decry_key_path, decry_key)


def exa_decry(decry_key_path, plain_path, cipher_path, index_length):
    cipher_list = read_seq(cipher_path)
    list_len = len(cipher_list)
    cipher_list = exa_remove_index(cipher_list, list_len, 5)
    decry_key_list = read_seq(decry_key_path)
    plain_list = rexor(decry_key_list, cipher_list)
    print(len(plain_list))
    coder.decode(plain_list, index_length, 'AGCT', plain_path)


def encry_rsa(plain_path, key_byte, cipher_path, seg_length):
    plain_bytes = trans.fasta_to_bytes(plain_path)
    # print(len(plain_bytes))
    size_pad_set = size_pad(plain_bytes)
    plain_bytes = size_pad_set['data']
    size_pad_num = size_pad_set['pad_num']
    rsa_key = RSA.importKey(key_byte)
    cipher = PKCS1_v1_5.new(rsa_key)
    cipher_bytes = bytearray()
    for i in range(0, len(plain_bytes), 200):
        cipher_bytes += (cipher.encrypt(plain_bytes[i:i + 200]))
    cipher_bytes_size = len(cipher_bytes)
    # print(cipher_bytes_size)
    seq_pad_set = seq_pad(cipher_bytes, seg_length)
    cipher_bytes = seq_pad_set['data']
    # print(len(cipher_bytes))

    trans.bytes_to_fasta(cipher_bytes, cipher_path, seg_length)
    return size_pad_num, cipher_bytes_size


def decry_rsa(plain_path, key_byte, cipher_path, seg_length, size_pad_num, cipher_bytes_size):
    cipher_bytes = trans.fasta_to_bytes(cipher_path)
    cipher_bytes = satisfy_bound(cipher_bytes, cipher_bytes_size)
    print(len(cipher_bytes))
    rsa_key = RSA.importKey(key_byte)
    plain = PKCS1_v1_5.new(rsa_key)
    plain_bytes = bytearray()
    for i in range(0, len(cipher_bytes), 256):
        print(i)
        plain_bytes += (plain.decrypt(cipher_bytes[i:i + 256], "error"))
    plain_bytes = remove_pad(plain_bytes, size_pad_num)
    trans.bytes_to_fasta(plain_bytes, plain_path, seg_length)


def add_hmac(plain_path, key_byte, cipher_path, seg_length):
    plain_bytes = trans.fasta_to_bytes(plain_path)
    hmac_object = hmac.new(key_byte, plain_bytes, digestmod=hashlib.sha256)
    hmac_value = hmac_object.digest()
    cipher_bytes = hmac_value + plain_bytes
    cipher_bytes = seq_pad(cipher_bytes, seg_length)
    trans.bytes_to_fasta(cipher_bytes, cipher_path, seg_length)
    return 0


def verify_hmac(cipher_path, key_byte):
    cipher_bytes = trans.fasta_to_bytes(cipher_path)
    cipher_bytes = remove_blank(cipher_bytes)
    plain_bytes = cipher_bytes[32:]
    ori_hmac_value = cipher_bytes[:32]
    hmac_object = hmac.new(key_byte, plain_bytes, digestmod=hashlib.sha256)
    plain_hmac_value = hmac_object.digest()
    return hmac.compare_digest(ori_hmac_value, plain_hmac_value)


def ecc_encry(plain_path, key_byte, cipher_path, seg_length):
    plain_bytes = trans.fasta_to_bytes(plain_path)
    hash_obj = hashlib.sha256(plain_bytes).digest()
    ecc_sign = key_byte.sign(hash_obj)
    cipher_bytes = ecc_sign + plain_bytes
    cipher_bytes = seq_pad(cipher_bytes, seg_length)
    trans.bytes_to_fasta(cipher_bytes, cipher_path, seg_length)
    return 0


def ecc_verify(cipher_path, key_byte):
    cipher_bytes = trans.fasta_to_bytes(cipher_path)
    cipher_bytes = remove_blank(cipher_bytes)
    plain_bytes = cipher_bytes[64:]
    ecc_sign = cipher_bytes[:64]
    hash_obj = hashlib.sha256(plain_bytes).digest()
    ver_key = key_byte.get_verifying_key()
    valid = ver_key.verify(ecc_sign, hash_obj)
    return valid


def string_to_sublists(input_string, sublist_size):
    # 计算可以生成多少个指定长度的子字符串
    num_full_sublists = len(input_string) // sublist_size
    # 创建一个列表，包含完整的子字符串
    sublists = [input_string[i:i + sublist_size] for i in range(0, num_full_sublists * sublist_size, sublist_size)]
    # 处理剩余的字符，如果它们不足以构成一个完整的子字符串
    remaining_chars = input_string[num_full_sublists * sublist_size:]
    if remaining_chars:
        sublists.append(remaining_chars)
    return sublists


def bytes_to_binary_string(data):
    binary_string = ''.join(format(byte, '08b') for byte in data)
    base_string = []
    for i in range(0, len(binary_string)):
        if binary_string[i] == '0':
            base_string.append('A')
        if binary_string[i] == '1':
            base_string.append('G')
    return base_string


def base_xor2(base1, base2):
    num1 = dic1.index(base1)
    num2 = dic1.index(base2)
    result_base = exor_table[num1][num2]
    return result_base


def seq_xor(state_seq, info_seq):
    result_seq = []
    for i in range(0, len(info_seq)):
        result_base = base_xor2(info_seq[i], state_seq[i])
        result_seq.append(result_base)
    result_seq = ''.join(result_seq)
    return result_seq


def caculate_seq_gc(seq):
    gc_count = 0.00
    for base in seq:
        if base == 'G' or base == 'C':
            gc_count += 1
    gc_percent = gc_count / len(seq)
    return gc_percent


def seq_gc_upper(key_seq, info_seq):
    key_seq_list = list(key_seq)
    result_seq = list(info_seq)
    for i in range(0, len(info_seq)):
        if result_seq[i] == 'G' or result_seq[i] == 'C':
            if key_seq_list[i] == 'A':
                key_seq_list[i] = 'C'
            # if key_seq_list[i] == 'T':
            #     key_seq_list[i] = 'G'
            # if key_seq_list[i] == 'C':
            #     key_seq_list[i] = 'A'
            if key_seq_list[i] == 'G':
                key_seq_list[i] = 'T'
            result_base = base_xor2(key_seq_list[i], info_seq[i])
            result_seq[i] = result_base
            seq_gc = caculate_seq_gc(''.join(result_seq))
            if seq_gc < 0.55:
                break
    result_seq = ''.join(result_seq)
    decry_keyseq = ''.join(key_seq_list)
    return decry_keyseq, result_seq


def seq_gc_control(key_seq, info_seq):
    seq_gc = caculate_seq_gc(info_seq)
    decry_keyseq, result_seq = key_seq, info_seq
    if seq_gc > 0.55:
        decry_keyseq, result_seq = seq_gc_upper(key_seq, info_seq)
    if seq_gc < 0.45:
        decry_keyseq, result_seq = seq_gc_lower(key_seq, info_seq)
    return decry_keyseq, result_seq


def seq_gc_lower(key_seq, info_seq):
    key_seq_list = list(key_seq)
    result_seq = list(info_seq)
    for i in range(0, len(info_seq)):
        if result_seq[i] == 'A' or result_seq[i] == 'T':
            if key_seq_list[i] == 'A':
                key_seq_list[i] = 'C'
            # if key_seq_list[i] == 'T':
            #     key_seq_list[i] = 'G'
            # if key_seq_list[i] == 'C':
            #     key_seq_list[i] = 'A'
            if key_seq_list[i] == 'G':
                key_seq_list[i] = 'T'
            result_base = base_xor2(key_seq_list[i], info_seq[i])
            result_seq[i] = result_base
            seq_gc = caculate_seq_gc(''.join(result_seq))
            if seq_gc > 0.45:
                break
    result_seq = ''.join(result_seq)
    decry_keyseq = ''.join(key_seq_list)
    return decry_keyseq, result_seq


def seq_exor(key_seq, info_seq, homolength):
    key_seq_list = list(key_seq)
    result_seq = list(seq_xor(key_seq, info_seq))
    for i in range(homolength, len(result_seq)):
        check_str = result_seq[i - homolength:i + 1]
        if len(set(check_str)) == 1:
            # print(result_seq)
            if key_seq_list[i] == 'A':
                key_seq_list[i] = 'C'
            if key_seq_list[i] == 'T':
                key_seq_list[i] = 'G'
            result_base = base_xor2(key_seq_list[i], info_seq[i])
            result_seq[i] = result_base
    result_seq = ''.join(result_seq)
    decry_keyseq = ''.join(key_seq_list)
    return decry_keyseq, result_seq


def exor(key_list, info_list, homolength):
    result_list = []
    decry_key_list = []
    for i in range(0, len(info_list)):
        decry_key_seq, result_seq = seq_exor(key_list[i], info_list[i], homolength)
        # decry_key_seq, result_seq = seq_gc_control(decry_key_seq, result_seq)
        result_list.append(result_seq)
        decry_key_list.append(decry_key_seq)
    return decry_key_list, result_list


def rexor(key_list, info_list):
    result_list = []
    for i in range(0, len(info_list)):
        if len(info_list[i]) < len(key_list[i]):
            info_list[i] = info_list[i] + 'A' * (len(key_list[i]) - len(info_list[i]))
        if len(info_list[i]) > len(key_list[i]):
            info_list[i] = info_list[i][0:len(key_list[i])]
        result_seq = seq_xor(key_list[i], info_list[i])
        result_list.append(result_seq)
    return result_list


def size_pad(byte_s):
    pad_num = 0
    while len(byte_s) % 16 != 0:
        byte_s += ' '.encode()
        pad_num += 1
    # print(len(byte_s))
    return {'data': byte_s, 'pad_num': pad_num}


# def remove_seq_pad(byte_s, num):

def seq_pad(byte_s, seg_length):
    pad_num = 0
    while len(byte_s) % (seg_length / 4) != 0:
        byte_s += '!'.encode()
        pad_num += 1
        # print(len(byte_s))
    return {'data': byte_s, 'pad_num': pad_num}


def remove_pad(byte_s, pad_num):
    new_bytes = byte_s
    if pad_num != 0:
        new_bytes = byte_s[:-pad_num]
    return new_bytes


def satisfy_bound(byte_s, cipher_bytes_size):
    while len(byte_s) < cipher_bytes_size:
        byte_s += '!'.encode()
    return byte_s[:cipher_bytes_size]


def ecc_generate_key():
    pri_key = SigningKey.generate(curve=NIST256p)
    return pri_key


ecc_key = ecc_generate_key()


def aes_keys_write(keys_path):
    iv = get_random_bytes(16)
    key_bytes = get_random_bytes(32)
    # print(iv, key_bytes)
    with open(keys_path, 'w') as f:
        iv_seq = trans.bytes2seq(iv)
        key_seq = trans.bytes2seq(key_bytes)
        f.write(iv_seq)
        f.write('\n')  # 添加换行符
        f.write(key_seq)

    return iv, key_bytes


def aes_keys_read(keys_path):
    with open(keys_path, 'r') as f:
        iv_seq = f.readline().strip()  # 读取第一行并去除换行符
        key_seq = f.readline().strip()  # 读取第二行并去除换行符
        iv = trans.fasta_to_bytes(iv_seq)
        key_bytes = trans.fasta_to_bytes(key_seq)
    # print(iv, key_bytes)
    return iv, key_bytes


def arc4_keys_write(keys_path):
    key_bytes = get_random_bytes(32)
    with open(keys_path, 'w') as f:
        key_seq = trans.bytes2seq(key_bytes)
        f.write(key_seq)
    return key_bytes


def arc4_keys_read(keys_path):
    with open(keys_path, 'r') as f:
        key_seq = f.readline().strip()
        key_bytes = trans.fasta_to_bytes(key_seq)
    return key_bytes


def rsa_keys_write(keys_path):
    rsa = RSA.generate(2048, Random.new().read)
    pri_key = rsa.exportKey()
    pub_key = rsa.public_key().exportKey()
    with open(keys_path, 'wb') as f:
        f.write(pri_key)
    return pub_key


def rsa_keys_read(keys_path):
    with open(keys_path, 'rb') as f:
        pri_key = f.read()
    return pri_key


def write_context(size_pad_num, cipher_bytes_size, context_path):
    with open(context_path, 'w') as f:
        f.write(str(size_pad_num))
        f.write('\n')  # 添加换行符
        f.write(str(cipher_bytes_size))


def read_context(context_path):
    with open(context_path, 'r') as f:
        size_pad_num = int(f.readline().strip())  # 读取第一行并去除换行符
        cipher_bytes_size = int(f.readline().strip())  # 读取第二行并去除换行符
    return size_pad_num, cipher_bytes_size


def log_base_4(x):
    return math.ceil(math.log(x) / math.log(4))


def add_index(seq_list):
    outlist = []
    list_len = len(seq_list)
    index_len = log_base_4(list_len)
    # print(index_len)
    dic_index = []
    for item in itertools.product('GATC', repeat=index_len):
        dic_index.append(''.join(item))
    for i in range(0, list_len):
        newseq = dic_index[i] + seq_list[i]
        outlist.append(newseq)
    return outlist


def remove_index(seq_list, list_len):
    outlist = []
    index_len = log_base_4(list_len)
    dic_index = []
    for item in itertools.product('GATC', repeat=index_len):
        dic_index.append(''.join(item))
    for i in range(0, list_len):
        index_str = dic_index[i]
        matched_seq = next((s for s in seq_list if s.startswith(index_str)), None)
        if matched_seq is None:
            matched_seq = find_most_similar(index_str, seq_list)
        outlist.append(matched_seq[index_len:])

    return outlist


def not_contains_n_consecutive_chars(string, homolength):
    # 遍历字符串，检查每个可能的子串
    for i in range(homolength, len(string)):
        check_str = string[i - homolength:i]
        # 如果找到n个连续的char，则返回False
        if len(set(check_str)) == 1:
            return False
    return True


def exa_index(seq_list, homolength):
    outlist = []
    list_len = len(seq_list)
    index_len = log_base_4(list_len)
    # print(index_len)
    dic_index = []
    for item in itertools.product('GATC', repeat=index_len):
        if not_contains_n_consecutive_chars(item, homolength):
            dic_index.append(''.join(item))
    for i in range(0, list_len):
        newseq = dic_index[i] + seq_list[i]
        outlist.append(newseq)
    return outlist


def exa_remove_index(seq_list, list_len, homolength):
    outlist = []
    index_len = log_base_4(list_len)
    dic_index = []
    for item in itertools.product('GATC', repeat=index_len):
        if not_contains_n_consecutive_chars(item, homolength):
            dic_index.append(''.join(item))
    for i in range(0, list_len):
        index_str = dic_index[i]
        matched_seq = next((s for s in seq_list if s.startswith(index_str)), None)
        if matched_seq is None:
            matched_seq = find_most_similar(index_str, seq_list)
        outlist.append(matched_seq[index_len:])

    return outlist


def compute_seq_runlength(seq):
    max_run_length = 1
    count = 1
    for i in range(len(seq) - 1):
        if seq[i] == seq[i + 1]:
            count += 1
        else:
            if count > max_run_length:
                max_run_length = count
            count = 1
    return max_run_length


def find_most_similar(target, string_list):
    # 初始化最小距离和最相似的字符串
    min_distance = float('inf')
    most_similar = None

    for string in string_list:
        # 计算Levenshtein距离
        distance = editdistance.distance(target, string)
        # 如果当前距离小于最小距离，则更新最小距离和最相似的字符串
        if distance < min_distance:
            min_distance = distance
            most_similar = string

    return most_similar


def gc_seq_compute(seq_list):
    seq_gc_satis = 0.00

    for seq in seq_list:
        gc_count = 0.00
        for base in seq:
            if base == 'G' or base == 'C':
                gc_count += 1
        gc_percent = gc_count / len(seq)
        if 0.6 > gc_percent > 0.4:
            seq_gc_satis += 1

    seq_gc = round(seq_gc_satis / len(seq_list), 4)
    return seq_gc


def run_length_compute(seq_list):
    max_run_length = 0
    for seq in seq_list:
        seq_run_length = compute_seq_runlength(seq)
        if seq_run_length > max_run_length:
            # print(seq)
            max_run_length = seq_run_length

    return max_run_length





if __name__ == "__main__":
    current_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "secure")

    keys_path = os.path.join(current_path, "ds-aes", "keys.txt")
    context_path = os.path.join(current_path, "ds-aes", "context.txt")
    plain_path = os.path.join(current_path, "input", "x_ray.png")
    recon_path = os.path.join(current_path, "ds-aes", "output-results.txt")
    cipher_fasta_path = os.path.join(current_path, "ds-aes", "x_ray.txt")
    decry_path = os.path.join(current_path, "ds-aes", "x_ray_decry.png")
    decry_key_path = os.path.join(current_path, "ds-aes", "decry_key.txt")
    err_cipher_path = os.path.join(current_path, "ds-aes", "x_ray_ks.txt")
    code_mode = AES.MODE_CBC
    #

    #
    segment_length = 128
    index_length = 8

    # aes encry
    # iv, key_bytes = aes_keys_write(keys_path)
    # print(iv, key_bytes)
    # size_pad_num, cipher_bytes_size = encry_aes(plain_fasta_path, iv, key_bytes, cipher_fasta_path, segment_length,
    #                                             index_length, code_mode)
    # write_context(size_pad_num, cipher_bytes_size, context_path)
    iv, key_bytes = aes_keys_read(keys_path)
    # key_bytes = ran1bit(key_bytes)
    # print(iv, key_bytes)
    # process = psutil.Process(os.getpid())
    # start_time = time.time()
    # start_memory = process.memory_info().rss
    size_pad_num, cipher_bytes_size = encry_aes(plain_path, iv, key_bytes, err_cipher_path, segment_length,
                                                index_length, code_mode)
    # write_context(size_pad_num, cipher_bytes_size, context_path)
    # size_pad_num, cipher_bytes_size = read_context(context_path)
    # decry_aes(decry_path, iv, key_bytes, recon_path, segment_length, index_length, code_mode, size_pad_num,
    #                     cipher_bytes_size)
    # end_memory = process.memory_info().rss
    # end_time = time.time()
    # print(f"内存使用：{round((end_memory - start_memory) / 1024 / 1024, 4)}MB")
    # print(f"运行时间：{round(end_time - start_time, 4)}S")
    # write_context(size_pad_num, cipher_bytes_size, context_path)
    #
    # # aes decry
    # size_pad_num, cipher_bytes_size = read_context(context_path)
    # decry_aes(decry_fasta_path, iv, key_bytes, recon_path, segment_length, code_mode, size_pad_num,
    #           cipher_bytes_size, 1456)

    # rsa encry
    # pub_key = rsa_keys_write(keys_path)
    # size_pad_num, cipher_bytes_size = encry_rsa(plain_fasta_path, pub_key, cipher_fasta_path, segment_length)
    # write_context(size_pad_num, cipher_bytes_size, context_path)
    # rsa decry
    # pri_key = rsa_keys_read(keys_path)
    # size_pad_num, cipher_bytes_size = read_context(context_path)
    # decry_rsa(decry_fasta_path, pri_key, recon_path, segment_length, size_pad_num, cipher_bytes_size)

    # add_hmac(plain_fasta_path, key_bytes, cipher_fasta_path, segment_length)
    # print(verify_hmac(cipher_fasta_path, key_bytes))
    #
    # ecc_encry(plain_fasta_path, ecc_key, cipher_fasta_path, segment_length)
    # print(ecc_verify(cipher_fasta_path, ecc_key))

    # arc4 encry
    # key_bytes = arc4_keys_write(keys_path)
    # size_pad_num, cipher_bytes_size = arc4_encry(plain_path, key_bytes, cipher_fasta_path, segment_length, index_length)
    # write_context(size_pad_num, cipher_bytes_size, context_path)
    #
    # key_bytes = arc4_keys_read(keys_path)

    # process = psutil.Process(os.getpid())
    # start_time = time.time()
    # start_memory = process.memory_info().rss
    # size_pad_num, cipher_bytes_size = arc4_encry(plain_path, key_bytes, err_cipher_path, segment_length, index_length)
    # write_context(size_pad_num, cipher_bytes_size, context_path)
    # size_pad_num, cipher_bytes_size = read_context(context_path)
    # arc4_decry(decry_path, key_bytes, recon_path, segment_length, size_pad_num, cipher_bytes_size, index_length)
    # end_memory = process.memory_info().rss
    # end_time = time.time()
    # print(f"内存使用：{round((end_memory - start_memory) / 1024 / 1024, 4)}MB")
    # print(f"运行时间：{round(end_time - start_time, 4)}S")

    # #
    #
    # exa encry
    # process = psutil.Process(os.getpid())
    # start_time = time.time()
    # start_memory = process.memory_info().rss
    # exa_encry(keys_path, plain_path, err_cipher_path, segment_length, decry_key_path, index_length)

    # cipher_list = read_seq(cipher_fasta_path)
    # print(run_length_compute(cipher_list))
    # print(gc_seq_compute(cipher_list))
    # exa_decry(decry_key_path, decry_path, recon_path, index_length)
    # end_memory = process.memory_info().rss
    # end_time = time.time()
    # print(f"内存使用：{round((end_memory - start_memory) / 1024 / 1024, 4)}MB")
    # print(f"运行时间：{round(end_time - start_time, 4)}S")
    # st.count_difference_file(recon_path, cipher_fasta_path)
    # st.caculate_difference_file(cipher_fasta_path, err_cipher_path)
    st.caculate_nbcr(cipher_fasta_path, err_cipher_path)
    st.caculate_baci(cipher_fasta_path, err_cipher_path)
