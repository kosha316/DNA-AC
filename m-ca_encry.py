import random

import mealy
import cellular_automata
from operator import *
import os
import eccorection
import securitytest
import hashlib
import trans_type as trans
import time
import psutil
import coder
import seq_methods

dic1 = ['A', 'T', 'C', 'G']


class MealyCellularAutomata(object):
    machine = mealy.Matter()
    ca = cellular_automata.CellularAutomata()

    def work(self, plain_seqlist, key_seq, key_path, homo_length, index_length):
        cipher_seqlist = []
        rnums_list = []

        for sequence in plain_seqlist:
            # Mealy模块
            mid_sequence = self.machine.work2(sequence, index_length)
            # CA模块
            end_sequence, rnums_str = self.ca.seq_encry(mid_sequence, key_seq, homo_length, 0)
            rnums_list.append(rnums_str)
            cipher_seqlist.append(end_sequence)
        rnums_seq = ''.join(rnums_list)
        write_rnums(key_path, rnums_seq)
        return cipher_seqlist

    def rework(self, cipher_seqlist, key_path, index_length, num_length):
        plain_seqlist = []
        key_seq = read_hash_key(key_path)
        rnums_list = get_rnums(key_path, num_length)
        count = 0
        for sequence in cipher_seqlist:
            # CA模块
            mid_sequence = self.ca.seq_decry(sequence, key_seq, 0, rnums_list[count])
            # Mealy模块
            end_sequence = self.machine.rework2(mid_sequence, index_length)
            plain_seqlist.append(end_sequence)
            count += 1
        return plain_seqlist

    # def work2(self, plain_seqlist, key_seq, homo_length, index_length):
    #     cipher_list = []
    #     self.machine.initiate_state = atcg2int(key_seq)
    #     # print(key_seq)
    #     # Mealy模块
    #     column_list = get_column(plain_seqlist)
    #     for i in range(index_length, len(column_list)):
    #         new_column = self.machine.work(column_list[i], 0, i)
    #         column_list[i] = new_column
    #     seq_list = get_column(column_list)
    #     new_key = int2atcg(self.machine.end_state)
    #
    #     # CA模块
    #
    #     for mid_sequence in seq_list:
    #         end_sequence = self.ca.seq_encry(mid_sequence, new_key, homo_length, index_length)
    #         cipher_list.append(end_sequence)
    #
    #     return cipher_list, new_key
    #
    # def rework2(self, cipher_seqlist, key_seq, index_length):
    #     plain_seqlist = []
    #     seqnum = 0
    #     # print(len(self.ca.list_cyclenum))
    #     # CA模块
    #     for sequence in cipher_seqlist:
    #         if len(sequence) < len(key_seq):
    #             sequence = sequence + 'A' * (len(key_seq) - len(sequence))
    #         if len(sequence) > len(key_seq):
    #             sequence = sequence[0:len(key_seq)]
    #         mid_sequence = self.ca.seq_decry(sequence, key_seq, index_length)
    #         plain_seqlist.append(mid_sequence)
    #         seqnum += 1
    #
    #     # Mealy模块
    #     column_list = get_column(plain_seqlist)
    #     for i in range(index_length, len(column_list)):
    #         new_column = self.machine.rework(column_list[i], i - index_length, 0)
    #         column_list[i] = new_column
    #     plain_seqlist = get_column(column_list)
    #
    #     return plain_seqlist


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


def get_column(seq_list):
    column_list = []
    tempseq = []
    for i in range(len(seq_list[0])):
        for j in range(len(seq_list)):
            tempseq.append(seq_list[j][i])
        tempseq = ''.join(tempseq)
        column_list.append(tempseq)
        tempseq = []
    return column_list


def int2atcg(numlist):
    baselist = []
    for num in numlist:
        base = dic1[num]
        baselist.append(base)
    baselist = ''.join(baselist)
    return baselist


def atcg2int(seq):
    numlist = []
    for base in seq:
        num = dic1.index(base)
        numlist.append(num)
    return numlist


def write_rnums(rnums_path, rnums_seq):
    with open(rnums_path, 'a') as f:
        f.write(rnums_seq)


def get_rnums(file_path, length):
    with open(file_path, 'r', encoding='utf-8') as file:  # 根据需要设置正确的编码
        lines = file.readlines()
        # 检查是否有至少两行数据
        if len(lines) >= 2:
            s = lines[1].strip()
    rnums_list = [s[i:i + length] for i in range(0, len(s), length)]
    # print(rnums_list)
    return rnums_list


def base_rotate(info_base, state_base):
    state_num = dic1.index(state_base)
    info_index = (dic1.index(info_base) + state_num) % 4
    new_base = dic1[info_index]
    return new_base


def random1base(seq):
    index = random.randint(0, 128)
    num = random.randint(0, 4)
    rbase = base_rotate(seq[index], dic1[num])
    # print(seq[index], rbase)
    return seq[:index] + rbase + seq[index + 1:]


def get_hash_key(file_path, key_path):
    # 创建一个SHA-512哈希对象
    sha512 = hashlib.sha512()
    # 需要哈希的数据
    with open(file_path, 'rb') as f:
        data = f.read()
        data = coder.ran1bit(data)
    # 使用update方法更新哈希对象
    sha512.update(data)
    # 使用hexdigest方法获取十六进制表示的哈希值
    hash_value = sha512.hexdigest()
    # 将哈希值转换为字符串
    keyseq = trans.bytes2seq(bytes.fromhex(hash_value))
    # print(hash_value.encode())
    with open(key_path, 'w') as f:
        f.write(keyseq)
        f.write('\n')
    return keyseq


def read_hash_key(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:  # 根据需要设置正确的编码
        lines = file.readlines()
    # print(lines[0].strip())
    return lines[0].strip()


def MCA_encry(MCA, plain_path, cipher_path, key_seq, key_path, homo_length, index_length, seq_length):
    plain_list = coder.encode(plain_path, seq_length, index_length, pad_str='AGCT')
    cipher_list = MCA.work(plain_list, key_seq, key_path, homo_length, index_length)
    write_seq(cipher_path, cipher_list)


def MCA_encry_compa(MCA, plain_path, cipher_path, key_seq, key_path, homo_length, index_length):
    plain_list = read_seq(plain_path)
    cipher_list = MCA.work(plain_list, key_seq, key_path, homo_length, index_length)
    write_seq(cipher_path, cipher_list)


def MCA_decry(MCA, decry_path, cipher_path, key_path, index_length):
    cipher_list = read_seq(cipher_path)
    plain_list = MCA.rework(cipher_list, key_path, index_length, num_length=3)
    coder.decode(plain_list, index_length, 'AGCT', decry_path)


def MCA_decry_compa(MCA, decry_path, cipher_path, key_path, index_length):
    cipher_list = read_seq(cipher_path)
    plain_list = MCA.rework(cipher_list, key_path, index_length, num_length=3)
    write_seq(decry_path, plain_list)


if __name__ == "__main__":
    current_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "secure")
    plain_path = os.path.join(current_path, "input", "x_ray.png")
    fasta_path = os.path.join(current_path, "input",  "x_ray.txt")
    cipher_path = os.path.join(current_path, "output",  "x_ray.txt")
    decry_path = os.path.join(current_path, "output",  "x_ray_decry.png")
    recon_path = os.path.join(current_path, "output",  "output-results.txt")
    key_path = os.path.join(current_path, "output",  "keys.txt")
    # rnums_path = os.path.join(current_path, "output", "other_code", "rnums.txt")
    # err_plain_path = os.path.join(current_path, "input", "err1base_file", '1.txt')
    err_cipher_path = os.path.join(current_path, "output", "x_ray_ks.txt")

    MCA = MealyCellularAutomata()
    # ECC = eccorection.ecc_code()

    # print(len(key_seq))
    # plain_list = read_seq(plain_path)
    # cipher_list = read_seq(cipher_path)
    key_seq = get_hash_key(plain_path, key_path)
    # key_seq = random1base(key_seq)
    #
    process = psutil.Process(os.getpid())
    start_time = time.time()
    start_memory = process.memory_info().rss
    # cipher_list = MCA.work(plain_list, key_seq, 4, 8)
    # decry_list = MCA.rework(cipher_list, key_seq, 8)
    MCA_encry(MCA, plain_path, err_cipher_path, key_seq, key_path, 4, 8, 128)
    # MCA_encry_compa(MCA, fasta_path, cipher_path, key_seq, key_path, 4, 8)
    # MCA_decry_compa(MCA, decry_path, recon_path, key_path, 8)
    # MCA_decry(MCA, decry_path, recon_path, key_path, 8)
    end_memory = process.memory_info().rss
    end_time = time.time()
    print(f"内存使用：{round((end_memory - start_memory) / 1024 / 1024, 4)}MB")
    print(f"运行时间：{round(end_time - start_time, 4)}S")
    # write_seq(cipher_path, cipher_list)
    # write_seq(decry_path, decry_list)

    #
    # eccencode_list = ECC.encode(cipher_list, 'ReedSolomon')
    # print(len(eccencode_list[0]))
    # write_seq(tihuan_cipher_path, cipher_list)
    # recon_list = read_seq(recon_path)
    # eccdecode_list = ECC.decode(recon_list, 'ReedSolomon')
    # print(eq(cipher_list, eccdecode_list))
    # print(len(eccdecode_list[0]))
    # cipher_list = read_seq(cipher_path)
    # recon_list = read_seq(recon_path)
    # decry_list = MCA.rework(recon_list, key_seq, 0)
    # write_seq(decry_path, decry_list)
    # securitytest.count_similarity_file(cipher_path, decry_path)
    # securitytest.count_similarity_list(plain_list, decry_list)
    # securitytest.caculate_difference_file(cipher_path, err_cipher_path)
    securitytest.caculate_nbcr(cipher_path, err_cipher_path)
    securitytest.caculate_baci(cipher_path, err_cipher_path)
