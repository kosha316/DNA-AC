import os
from Bio.Seq import Seq
import itertools
from operator import *

dic = ['A', 'T', 'C', 'G']


# 定义有限状态机
class Matter(object):
    # 状态表
    state_table = [
        [0, 3, 1, 2],
        [2, 1, 3, 0],
        [3, 0, 2, 1],
        [2, 1, 0, 3]
    ]
    # transposed_state = list(zip(*state_table))

    # 输出表
    output_table = [
        ['T', 'C', 'G', 'A'],
        ['T', 'C', 'T', 'G'],
        ['G', 'A', 'G', 'T'],
        ['C', 'A', 'C', 'A'],
    ]
    # transposed_output = list(zip(*output_table))

    # 逆状态表
    state_table_r = [
        [2, 0, 3, 1],
        [3, 2, 1, 0],
        [0, 1, 3, 2],
        [3, 1, 0, 2]
    ]

    # 逆输出表
    output_table_r = [
        ['T', 'A', 'C', 'G'],
        ['T', 'G', 'T', 'C'],
        ['G', 'A', 'A', 'C'],
        ['G', 'C', 'T', 'A'],
    ]

    xor_table = [
        ['A', 'T', 'C', 'G'],
        ['T', 'A', 'G', 'C'],
        ['C', 'G', 'A', 'T'],
        ['G', 'C', 'T', 'A']
    ]

    initiate_state = []
    end_state = []

    def base_xor2(self, base1, base2):
        num1 = dic.index(base1)
        num2 = dic.index(base2)
        result_base = self.xor_table[num1][num2]
        return result_base

    def work2(self, input_seq, index_length):
        seq_len = len(input_seq)
        state = dic.index(input_seq[seq_len - 1])
        out_list = list(input_seq)
        for i in range(index_length, seq_len):
            output_char = self.output_table[state][dic.index(input_seq[i])]
            out_list[i] = output_char
            state = self.state_table[state][dic.index(input_seq[i])]
            # state = dic.index(self.base_xor2(dic[state], output_char))

            # print(state, end='')
        # out_list.append('T')
        out_list.append(dic[state])
        out_list.append('C')
        output_seq = ''.join(out_list)
        output_seq = output_seq[::-1]

        # print(self.end_state)
        return output_seq

    # def work(self, input_seq, index_length, state_index):
    #     state = self.initiate_state[state_index]
    #     out_list = list(input_seq)
    #     for i in range(index_length, len(input_seq)):
    #         output_char = self.output_table[state][dic.index(input_seq[i])]
    #         out_list[i] = output_char
    #         state = self.state_table[state][dic.index(input_seq[i])]
    #         # state = dic.index(self.base_xor2(dic[state], output_char))
    #
    #         # print(state, end='')
    #     output_seq = ''.join(out_list)
    #     output_seq = output_seq[::-1]
    #
    #     self.end_state.append(state)
    #     # print(self.end_state)
    #     return output_seq
    #
    # def rework(self, input_seq, seq_num, index_length):
    #     state = self.end_state[seq_num]
    #     out_list = list(input_seq)
    #     for i in range(len(input_seq) - index_length):
    #         # state = dic.index(self.base_xor2(dic[state], input_seq[i]))
    #         output_char = self.output_table_r[state][dic.index(input_seq[i])]
    #         out_list[i] = output_char
    #         state = self.state_table_r[state][dic.index(input_seq[i])]
    #
    #     output_seq = ''.join(out_list)
    #     output_seq = output_seq[::-1]
    #     return output_seq

    def rework2(self, input_seq, index_length):
        input_seq = input_seq[1:]
        state = dic.index(input_seq[0])
        seq_len = len(input_seq)
        out_list = list(input_seq)
        for i in range(1, seq_len - index_length):
            # state = dic.index(self.base_xor2(dic[state], input_seq[i]))
            output_char = self.output_table_r[state][dic.index(input_seq[i])]
            out_list[i] = output_char
            state = self.state_table_r[state][dic.index(input_seq[i])]

        output_seq = ''.join(out_list[1:])
        output_seq = output_seq[::-1]
        return output_seq


def int2atcg(numlist):
    baselist = []
    for num in numlist:
        base = dic[num]
        baselist.append(base)
    baselist = ''.join(baselist)
    return baselist


def atcg2int(seq):
    numlist = []
    for base in seq:
        num = dic.index(base)
        numlist.append(num)
    return numlist


if __name__ == "__main__":
    machine = Matter()
    infoseq = "GGGGAAGCCGATTGACAGCGCGCGCACTGTGAGCCCCATCGGGGAAGTGTCCGCTCTCCGCGCGACTAACTACTACACAATAGTGACTAGCACCGAGCTCTAATTTCCCATTTGAGGCAGACTGTGCCTGGT"
    machine.initiate_state = atcg2int(
        "AAAAAAGCCGATTGACAGCGCGCGCACTGTGAGCCCCATCGGGGAAGTGTCCGCTCTCCGCGCGACTAACTACTACACAATAGTGACTAGCACCGAGCTCTAATTTCCCATTTGAGGCAGACTGTGCCTGGT")
    encry_seq = machine.work2(infoseq, 8)
    decry_seq = machine.rework2(encry_seq, 8)
    print(encry_seq)
    print(decry_seq)
    print(eq(decry_seq, infoseq))
