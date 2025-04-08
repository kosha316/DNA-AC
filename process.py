import os
from Bio.Seq import Seq
import itertools
from operator import *
import securitytest as st

dic3 = ['A', 'T', 'C', 'G']


class Constants:
    dic1 = list()
    dic2 = list()

    def __init__(self):
        for item in itertools.product('GATC', repeat=4):
            self.dic1.append(''.join(item))
        for item in itertools.product('TCGA', repeat=3):
            self.dic2.append(''.join(item))


# 序列关系计数实现
def seq_substr_count(seq_list, str_len):
    substr_count = {}
    homo_length = 5
    repeat_count = 0
    reverse_count = 0
    homo_count = 0
    list1 = list()
    list2 = list()
    for item in seq_list:
        sequence = Seq(item)
        reverse_seq = sequence.reverse_complement()
        for i in range(0, len(sequence) - str_len):
            reverse_str = reverse_seq[i: i + str_len]
            if reverse_str not in list2:
                list2.append(reverse_str)
                tmpcount = sequence.count(reverse_str)
                if tmpcount > 0:
                    reverse_count += tmpcount
                else:
                    i = i + str_len

        for i in range(0, len(sequence) - str_len):
            str_tmp = sequence[i: i + str_len]
            if str_tmp not in list1:
                list1.append(str_tmp)
                tmpcount = sequence.count(str_tmp) - 1
                if tmpcount > 0:
                    repeat_count += tmpcount
                else:
                    i = i + str_len

        for j in range(0, len(sequence) - homo_length):
            str_tmp = sequence[j: j + homo_length + 1]
            if str_tmp == 'A' * (homo_length + 1) or str_tmp == 'T' * (homo_length + 1) or str_tmp == 'C' * (
                    homo_length + 1) or str_tmp == 'G' * (homo_length + 1):
                homo_count += 1
                j += homo_length + 1

    substr_count['repeat'] = repeat_count
    substr_count['reverse_comp'] = reverse_count
    substr_count['homo'] = homo_count
    return substr_count


def seq_substr_process(sequence: Seq, str_len):
    list1 = list()
    list2 = list()
    reverse_seq = sequence.reverse_complement()

    # # 处理互补子序列
    # for i in range(0, len(sequence) - str_len):
    #     reverse_str = reverse_seq[i: i + str_len]
    #     if reverse_str not in list2:
    #         list2.append(reverse_str)
    #         tmpcount = sequence.count(reverse_str)
    #         if tmpcount > 0:
    #             reverse_count += tmpcount
    #         else:
    #             i = i + str_len
    # # 处理重复子序列
    # for i in range(0, len(sequence) - str_len):
    #     str_tmp = sequence[i: i + str_len]
    #     if str_tmp not in list1:
    #         list1.append(str_tmp)
    #         tmpcount = sequence.count(str_tmp) - 1
    #         if tmpcount > 0:
    #             repeat_count += tmpcount
    #         else:
    #             i = i + str_len

    return 0


# 添加索引
def add_index(seq_list, index_length):
    index_dic = list()
    for item in itertools.product('CTGA', repeat=index_length):
        index_dic.append(''.join(item))
    for i in range(len(seq_list)):
        item = seq_list[i]
        seq_list[i] = index_dic[i] + item
    return seq_list


# 去除索引
def remove_index(seq_list, index_length):
    for i in range(len(seq_list)):
        item = seq_list[i]
        seq_list[i] = item[index_length:]
    return seq_list


# 将文件转换为序列对象
def read_seq(file_path):
    with open(file_path, "r", encoding='utf-8') as f:
        fasta_segments = f.read()
        fasta_segments = fasta_segments.strip('\n')
        seq_list = fasta_segments.split('\n')
        f.close()
        # print(len(seq_list))
        return seq_list


def write_seq(file_path, seq_list):
    with open(file_path, 'w', encoding='utf-8') as f:
        for item in seq_list:
            f.write(item)
            f.write('\n')
    f.close()


# 移除子串
def remove_substr(sequence: Seq, position, str_length):
    return sequence[:position] + sequence[position + str_length:]


# 对均聚物进行替换
def homo_replace(sequence: Seq, position, homo_length):
    cons = Constants()
    pos_str = cons.dic1[position]
    homo_char = sequence[position]
    sequence = 'A' + homo_char + pos_str + remove_substr(sequence, position, homo_length + 1)
    return sequence


def homo_unreplace(sequence: Seq, homo_length):
    cons = Constants()
    while sequence[0] == 'A':
        homo_char = sequence[1]
        pos_str = sequence[2:6]
        position = cons.dic1.index(pos_str)
        sequence = remove_substr(sequence, 0, homo_length + 1)
        sequence = sequence[:position] + homo_char * (homo_length + 1) + sequence[position:]
    return sequence


# 对互补子序列进行替换
def revecopl_replace(sequence: Seq, position, str_length):
    return 0


# 对重复子序列进行替换
def repeat_replace(seq, ref_seq, str_length):
    i = 0
    length = len(seq)
    while i + 8 < length:
        subseq = str(seq[i:i + 8])
        position1 = ref_seq.find(subseq)
        position2 = int(i / 8)
        if 0 <= position1:
            pos_ref = cons.dic1[position1]
            pos_seq = cons.dic2[position2]
            new_seq = remove_substr(seq, i, str_length)
            seq = 'C' + pos_ref + pos_seq + new_seq
            # length -= 8
        i += 8

    return seq


# 对重复子序列进行逆替换
def repeat_unreplace(seq, ref_seq, str_length):
    while seq[0] == 'C':
        position1 = seq[1:5]
        position2 = seq[5:8]
        ref_pos = cons.dic1.index(position1)
        seq_pos = cons.dic2.index(position2) * 8
        tmpseq = ref_seq[ref_pos:ref_pos + str_length]
        seq = remove_substr(seq, 0, str_length)
        seq = seq[:seq_pos] + tmpseq + seq[seq_pos:]
    return seq


# 对反向子序列进行替换
def recompl_replace(seq, ref_seq, str_length):
    new_seq = ''
    i = 0
    length = len(seq)
    while i + 8 < length:
        tmpseq = Seq(seq[i:i + 8])
        subseq = str(tmpseq.reverse_complement())
        position1 = ref_seq.find(subseq)
        position2 = int(i / 8)
        if 0 <= position1:
            pos_ref = cons.dic1[position1]
            pos_seq = cons.dic2[position2]
            new_seq = remove_substr(seq, i, str_length)
            seq = 'G' + pos_ref + pos_seq + new_seq
            # length -= 8
        i += 8

    return seq


# 反向互补逆替换
def recompl_unreplace(seq, ref_seq, str_length):
    while seq[0] == 'G':
        position1 = seq[1:5]
        position2 = seq[5:8]
        ref_pos = cons.dic1.index(position1)
        seq_pos = cons.dic2.index(position2) * 8
        tmpseq = Seq(ref_seq[ref_pos:ref_pos + str_length])
        tmpseq = tmpseq.reverse_complement()
        seq = remove_substr(seq, 0, str_length)
        seq = seq[:seq_pos] + tmpseq + seq[seq_pos:]
    return seq


# 预处理
def preprocess(seq_list):
    for i in range(len(seq_list)):
        item = 'TT' + seq_list[i]
        seq_list[i] = item
    return seq_list


# 去除预处理
def unpreprocess(seq_list):
    for i in range(len(seq_list)):
        item = seq_list[i][2:]
        seq_list[i] = item
    return seq_list


# 对重复子序列进行处理
def repeat_process(seq_list, str_length):
    ref_seq = seq_list[0][10:] + seq_list[1][10:]
    ref_seq = ref_seq[:256]
    for i in range(2, len(seq_list)):
        item = seq_list[i]
        seq_list[i] = ''.join(repeat_replace(item, ref_seq, str_length))

    return seq_list


# 对重复子序列进行恢复处理
def repeat_unprocess(seq_list, str_length):
    ref_seq = seq_list[0][10:] + seq_list[1][10:]
    ref_seq = ref_seq[:256]
    for i in range(2, len(seq_list)):
        item = seq_list[i]
        seq_list[i] = ''.join(repeat_unreplace(item, ref_seq, str_length))
    return seq_list


# 对反向互补子序列进行处理
def recompl_process(seq_list, str_length):
    ref_seq = seq_list[0][10:] + seq_list[1][10:]
    ref_seq = ref_seq[:256]
    for i in range(2, len(seq_list)):
        item = seq_list[i]
        seq_list[i] = ''.join(recompl_replace(item, ref_seq, str_length))
    return seq_list


# 对反向互补子序列进行恢复处理
def recompl_unprocess(seq_list, str_length):
    ref_seq = seq_list[0][10:] + seq_list[1][10:]
    ref_seq = ref_seq[:256]
    for i in range(2, len(seq_list)):
        item = seq_list[i]
        seq_list[i] = ''.join(recompl_unreplace(item, ref_seq, str_length))
    return seq_list


# 处理均聚物
def homo_process(seq_list, homo_length):
    for i in range(len(seq_list)):
        sequence = seq_list[i]
        for j in range(0, len(sequence) - homo_length):
            str_tmp = sequence[j: j + homo_length + 1]
            if str_tmp == 'A' * (homo_length + 1) or str_tmp == 'T' * (homo_length + 1) or str_tmp == 'C' * (
                    homo_length + 1) or str_tmp == 'G' * (homo_length + 1):
                sequence = homo_replace(sequence, j, homo_length)
                j += homo_length + 1
        seq_list[i] = sequence
    return seq_list


# 对均聚物进行逆处理
def homo_unprocess(seq_list, homo_length):
    for i in range(len(seq_list)):
        item = seq_list[i]
        seq_list[i] = ''.join(homo_unreplace(item, homo_length))
    return seq_list


# 对文件中情况不同子序列关系计数
def file_count(seq_list, str_len):
    sum = {'repeat': 0, 'reverse_comp': 0, 'homo': 0}
    for item in seq_list:
        item = Seq(item)
        sum['repeat'] += seq_substr_count(item, str_len)['repeat']
        sum['reverse_comp'] += seq_substr_count(item, str_len)['reverse_comp']
        sum['homo'] += seq_substr_count(item, str_len)['homo']
    return sum


# 修改序列情况
def seq_count(seq_list):
    seq_sum = len(seq_list)
    seq_repeat = 0
    seq_recompl = 0
    seq_homo = 0
    for i in range(len(seq_list)):
        if seq_list[i][0] == 'C':
            seq_repeat += 1
        if seq_list[i][0] == 'G':
            seq_recompl += 1
        if seq_list[i][0] == 'A':
            seq_homo += 1
    print("seq_sum:", seq_sum)
    print("seq_repeat:", seq_repeat)
    print("seq_recompl:", seq_recompl)
    print("seq_homo:", seq_homo)
    return 0


#
def file_encode(seq_list, str_len, homo_length, file_path):
    # seq_list = add_index(seq_list, index_length=8)
    seq_list = preprocess(seq_list)
    # seq_list = repeat_process(seq_list, str_len)
    # seq_list = recompl_process(seq_list, str_len)
    seq_list = homo_process(seq_list, homo_length)
    # seq_list = homo_unprocess(seq_list, homo_length)
    seq_list = homo_process(seq_list, homo_length)
    # seq_count(seq_list)
    with open(file_path, 'w', encoding='utf-8') as f:
        for item in seq_list:
            f.write(item)
            f.write('\n')
    f.close()
    # print(seq_list)
    # write_seq(file_path, seq_list)


def homo_encode(seq_list, homo_length):
    seq_list = preprocess(seq_list)
    seq_list = homo_process(seq_list, homo_length)
    seq_list = homo_process(seq_list, homo_length)
    return seq_list


# 解码
def file_decode(seq_list, str_len, homo_length, file_path):
    seq_list = homo_unprocess(seq_list, homo_length)
    # seq_list = homo_unprocess(seq_list, homo_length)
    seq_list = homo_unprocess(seq_list, homo_length)
    # seq_list = recompl_unprocess(seq_list, str_len)
    # seq_list = repeat_unprocess(seq_list, str_len)
    seq_list = unpreprocess(seq_list)
    # seq_count(seq_list)
    with open(file_path, 'w', encoding='utf-8') as f:
        for item in seq_list:
            f.write(item)
            f.write('\n')
    f.close()


def homo_decode(seq_list, homo_length):
    seq_list = homo_unprocess(seq_list, homo_length)
    seq_list = homo_unprocess(seq_list, homo_length)
    seq_list = unpreprocess(seq_list)
    return seq_list


def find_homo(seq_list, homo_length):
    for i in range(len(seq_list)):
        sequence = seq_list[i]
        for j in range(0, len(sequence) - homo_length):
            str_tmp = sequence[j: j + homo_length + 1]
            if str_tmp == 'A' * (homo_length + 1) or str_tmp == 'T' * (homo_length + 1) or str_tmp == 'C' * (
                    homo_length + 1) or str_tmp == 'G' * (homo_length + 1):
                print("行:", i + 1, "  ", "列:", j + 1)
                j += homo_length + 1
    return 0


def seq_compare(base_data, unprocessed_data, processed_data):
    for i in range(len(base_data)):
        if base_data[i] != unprocessed_data[i]:
            print(i + 1)
            print(processed_data[i][0])
    return False


def gc_compute(fasta_segments):
    gc_count = 0.00
    base_count = 0.00

    for seq in fasta_segments:
        for base in seq:
            if base in {'A', 'T', 'C', 'G'}:
                base_count += 1.0
                if base in {'C', 'G'}:
                    gc_count += 1.0
    return round(gc_count / base_count, 6)


def run_length_compute(seq_list):
    max_run_length = 0
    count = 1

    for fasta_segments in seq_list:
        for i in range(len(fasta_segments) - 1):
            if fasta_segments[i] == fasta_segments[i + 1]:
                count += 1
            else:
                if count > max_run_length:
                    max_run_length = count
                count = 1
    return max_run_length


if __name__ == "__main__":
    cons = Constants()
    current_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "secure")

    cipher_path = os.path.join(current_path, "input", "time_memo_files", "197.txt")
    cipher_data = read_seq(cipher_path)
    processed_path = os.path.join(current_path, "process", "processed.fasta")

    # mealy_machine = Matter()

    # print(cipher_data[1])
    # print(mealied_data[1])

    unprocessed_path = os.path.join(current_path, "process", "unprocessed.fasta")
    # aes_cipher_data = read_seq(aes_cipher_path)
    # rsa_cipher_data = read_seq(rsa_cipher_path)

    # file_encode(cipher_data, 8, 5, processed_path)
    processed_data = read_seq(processed_path)
    # file_decode(processed_data, 8, 5, unprocessed_path)
    # # unprocessed_data = read_seq(unprocessed_path)
    find_homo(processed_data, 5)
    # st.count_similarity_file(cipher_path, unprocessed_path)
    # print("plain and cipher file content is same?", trans.compare_files(cipher_path, unprocessed_path))
    # seq_compare(rsa_cipher_data, unprocessed_data, processed_data)
    # print(cons.dic1)
    # print(cons.dic2)
    # dic2 = seq_substr_count(aes_cipher_data, str_len=8)
    # print(dic2)
