import os
import random

import editdistance
import seq_methods

dic1 = ['A', 'T', 'C', 'G']


def count_bases_changes(str1, str2):
    changes = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            changes += 1

    return changes


def count_bases_intensity(str1, str2):
    intensity = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            if str1[i] == 'Z' or str2[i] == 'Z':
                i += 1
            else:
                intensity += abs(dic1.index(str1[i]) - dic1.index(str2[i]))

    return intensity


def count_edit_distance(target, string):
    distance = editdistance.distance(target, string)
    return distance


# 读取序列
def read_seq(file_path):
    with open(file_path, "r", encoding='utf-8') as f:
        fasta_segments = f.read()
        fasta_segments = fasta_segments.strip('\n')
        seq_list = fasta_segments.split('\n')
        f.close()
        return seq_list


#
def count_changes_file(file1_path, file2_path):
    # 读取两个文件中的序列
    seqlist_1 = read_seq(file1_path)
    seqlist_2 = read_seq(file2_path)

    return count_changes_list(seqlist_1, seqlist_2)


def count_intensity_file(file1_path, file2_path):
    # 读取两个文件中的序列
    seqlist_1 = read_seq(file1_path)
    seqlist_2 = read_seq(file2_path)

    return count_intensity_list(seqlist_1, seqlist_2)

def count_changes_list(seqlist_1, seqlist_2):
    changes = 0.00
    seq_num = float(len(seqlist_1))
    base_num = seq_num * len(seqlist_1[0])
    #
    for i in range(len(seqlist_1)):
        changes += count_bases_changes(seqlist_1[i], seqlist_2[i])

    changes_rate = round(changes/base_num, 4)
    print('NBCR: ', changes_rate)

    return changes_rate


def count_intensity_list(seqlist_1, seqlist_2):
    intensity = 0.00
    seq_num = float(len(seqlist_1))
    base_num = seq_num * len(seqlist_1[0])
    #
    for i in range(len(seqlist_1)):
        intensity += count_bases_intensity(seqlist_1[i], seqlist_2[i])

    intensity_rate = round(intensity/(base_num*3), 4)
    print('BACI: ', intensity_rate)

    return intensity_rate


def caculate_similarity_list(seqlist_1, seqlist_2):
    distance_sum = 0.00
    same_seq = 0.00
    seq_num = float(len(seqlist_1))
    base_num = seq_num * len(seqlist_1[0])
    # 计算总编辑距离
    for i in range(len(seqlist_1)):
        seq_distance = count_edit_distance(seqlist_1[i], seqlist_2[i])
        distance_sum += seq_distance
        if seq_distance == 0:
            same_seq += 1
    base_similarity = round((base_num - distance_sum) / base_num, 4)
    seq_similarity = round(same_seq / seq_num, 4)
    print('similariy ', 'base:', base_similarity, 'seq:', seq_similarity)

    return base_similarity, seq_similarity


def caculate_similarity_file(file1_path, file2_path):
    # 读取两个文件中的序列
    seqlist_1 = read_seq(file1_path)
    seqlist_2 = read_seq(file2_path)

    return caculate_similarity_list(seqlist_1, seqlist_2)


def caculate_difference_file(file1_path, file2_path):
    radio_base, radio_seq = caculate_similarity_file(file1_path, file2_path)
    base_difference = round(1 - radio_base, 4)
    seq_difference = round(1 - radio_seq, 4)
    print("error ", 'base:', base_difference, 'seq:', seq_difference)
    return base_difference, seq_difference


def base_rotate(info_base, state_base):
    state_num = dic1.index(state_base)
    info_index = (dic1.index(info_base) + state_num) % 4
    new_base = dic1[info_index]
    return new_base


def ran_err1base(seq_list):
    seq_index = random.randint(0, len(seq_list) - 1)
    base_pos = random.randint(0, len(seq_list[0]) - 1)
    seq = list(seq_list[seq_index])
    err_base = base_rotate(seq[base_pos], 'C')
    seq[base_pos] = err_base
    seq = ''.join(seq)
    seq_list[seq_index] = seq
    return seq_list


def caculate_baci(file1_path, file2_path):
    # return caculate_difference_file(file1_path, file2_path)
    return count_intensity_file(file1_path, file2_path)

def caculate_nbcr(file1_path, file2_path):
    return count_changes_file(file1_path, file2_path)

if __name__ == "__main__":
    current_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "secure")
    fileone_path = os.path.join(current_path, "input", "x_ray_base.txt")
    filetwo_path = os.path.join(current_path, "input", "err1base.txt")
    path1 = r'C:\Users\Administrator\Desktop\papers\codes\IM-Codec-master\output\output_vitro.txt'
    path2 = r'C:\Users\Administrator\Desktop\papers\codes\IM-Codec-master\output\output_vitro_4.txt'
    # ran_err1base(fileone_path, filetwo_path)
    caculate_nbcr(path1, path2)
    caculate_baci(path1, path2)
