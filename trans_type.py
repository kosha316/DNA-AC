import sys
import operator
import numpy
import os

f2b_rule = {'A': 0b00, 'T': 0b01, 'C': 0b10, 'G': 0b11}
b2f_rule = {0b00: 'A', 0b01: 'T', 0b10: 'C', 0b11: 'G'}


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


def fasta_to_bytes(fasta_list):
    fasta_segments = ''.join(fasta_list)
    dec_segments = []
    new_bytes = []
    for base in fasta_segments:
        if base in {'A', 'T', 'C', 'G'}:
            dec_segments.append(f2b_rule[base])
    # print(dec_segments)
    # print("---------------分隔--------------------")

    byte_s = bytearray(dec_segments)
    # print(len(byte_s))
    # while len(byte_s) % 4 != 0:
    #     byte_s.append(' ')

    for i in range(len(byte_s)):
        if i % 4 == 3 and i != 0:
            new_bytes.append(byte_s[i - 3] * 64 + byte_s[i - 2] * 16 + byte_s[i - 1] * 4 + byte_s[i])
    # print(count)
    # print(new_bytes)
    # print("---------------分隔--------------------")
    # print(len(new_bytes))
    # print("---------------分隔--------------------")
    # print(byte_segments.__sizeof__())
    return bytearray(new_bytes)


def bytes_to_fasta(byte_s, segment_length):
    new_bytes = []
    # print(byte_s)
    for i in range(len(byte_s)):
        new_bytes.append(int((byte_s[i] & 0b11000000) / 64))
        new_bytes.append(int((byte_s[i] & 0b00110000) / 16))
        new_bytes.append(int((byte_s[i] & 0b00001100) / 4))
        new_bytes.append(int(byte_s[i] & 0b00000011))
    # print(new_bytes)
    fasta_segments = []
    new_dec_segments = list(new_bytes)
    # print(operator.eq(dec_segments,new_dec_segments))
    count = 0
    for base in new_dec_segments:
        # if base in {0, 1, 2, 3}:
        # if base not in {0, 1, 2, 3}:
        #     print("base not in 0,1,2,3")
        #     sys.exit()
        fasta_segments.append(b2f_rule[base])
        count += 1
        if count % segment_length == 0:
            fasta_segments.append('\n')
    str_tmp = ''.join(fasta_segments)
    str_tmp = str_tmp.strip('\n')
    # print(len(str_tmp))
    # print(dec_segments.__sizeof__())
    # print(byte_segments.__sizeof__())
    fasta_list = str_tmp.split('\n')

    return fasta_list


def bytes2seq(byte_s):
    new_bytes = []
    # print(byte_s)
    for i in range(len(byte_s)):
        new_bytes.append(int((byte_s[i] & 0b11000000) / 64))
        new_bytes.append(int((byte_s[i] & 0b00110000) / 16))
        new_bytes.append(int((byte_s[i] & 0b00001100) / 4))
        new_bytes.append(int(byte_s[i] & 0b00000011))
    # print(new_bytes)
    fasta_segments = []
    new_dec_segments = list(new_bytes)
    for base in new_dec_segments:
        # if base in {0, 1, 2, 3}:
        # if base not in {0, 1, 2, 3}:
        #     print("base not in 0,1,2,3")
        #     sys.exit()
        fasta_segments.append(b2f_rule[base])
    return ''.join(fasta_segments)






def compare_files(file1, file2):
    with open(file1, 'rb') as f1, open(file2, 'rb') as f2:
        while True:
            byte1 = f1.read(1)
            byte2 = f2.read(1)
            if byte1 != byte2:
                return False
            if not byte1:
                break
    return True


if __name__ == "__main__":
    current_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    plain_fasta_path = os.path.join(current_path, "secure", "arc4", "monalisa.txt")
    bin_path = os.path.join(current_path, "secure", "arc4", "monalisa.bin")

    byte_segments = fasta_to_bytes(read_seq(plain_fasta_path))
    with open(bin_path, 'wb') as f:
        f.write(byte_segments)
    f.close()
    # print(byte_segments)
    # bytes_to_fasta(byte_segments, cipher_fasta_path, 128)
    # print(compare_files(plain_fasta_path, cipher_fasta_path))
