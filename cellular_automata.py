import os
import Bio
from Bio.Seq import Seq
import itertools
from operator import *
import pandas as pd

dic1 = ['A', 'T', 'C', 'G']


class CellularAutomata(object):
    xor_table = [
        ['A', 'T', 'C', 'G'],
        ['T', 'A', 'G', 'C'],
        ['C', 'G', 'A', 'T'],
        ['G', 'C', 'T', 'A']
    ]

    list_cyclenum = []

    def base_xor2(self, base1, base2):
        num1 = dic1.index(base1)
        num2 = dic1.index(base2)
        result_base = self.xor_table[num1][num2]
        return result_base

    def base_xor3(self, base1, base2, base3):
        mid_base = self.base_xor2(base1, base2)
        result_base = self.base_xor2(mid_base, base3)
        return result_base

    def rule_90(self, state_seq):
        new_state = list(state_seq)
        seq_length = len(state_seq)
        for i in range(0, seq_length):
            if i == 0:
                new_state[i] = self.base_xor2(state_seq[seq_length - 1], state_seq[i + 1])
            if i == seq_length - 1:
                new_state[i] = self.base_xor2(state_seq[i - 1], state_seq[0])
            else:
                new_state[i] = self.base_xor2(state_seq[i - 1], state_seq[i + 1])
        new_state = ''.join(new_state)
        return new_state

    def rule_150(self, state_seq):
        new_state = list(state_seq)
        seq_length = len(state_seq)
        for i in range(0, seq_length):
            if i == 0:
                new_state[i] = self.base_xor3(state_seq[seq_length - 1], state_seq[i], state_seq[i + 1])
            if i == seq_length - 1:
                new_state[i] = self.base_xor3(state_seq[i - 1], state_seq[i], state_seq[0])
            else:
                new_state[i] = self.base_xor3(state_seq[i - 1], state_seq[i], state_seq[i + 1])
        new_state = ''.join(new_state)
        return new_state

    def not_satisfy_gc(self, seq):
        # gc_percent = 0.00
        gc_count = 0.00
        for base in seq:
            if base == 'C' or base == 'G':
                gc_count += 1.0
        gc_percent = gc_count / float(len(seq))
        if 0.6 >= gc_percent >= 0.4:
            return 0
        else:
            return 1

    def not_satify_random(self, seq):
        g_count = 0.00
        c_count = 0.00
        a_count = 0.00
        t_count = 0.00
        for base in seq:
            if base == 'A':
                a_count += 1.0
            if base == 'T':
                t_count += 1.0
            if base == 'C':
                c_count += 1.0
            if base == 'G':
                g_count += 1.0

        a_percent = a_count / (a_count + t_count)
        c_percent = c_count / (c_count + g_count)
        if 0.6 >= a_percent >= 0.4 or 0.6 >= c_percent >= 0.4:
            return 0
        else:
            return 1

    def not_satisfy_homo(self, seq, homo_length):
        for j in range(0, len(seq) - homo_length):
            str_tmp = seq[j: j + homo_length + 1]
            if str_tmp == 'A' * (homo_length + 1) or str_tmp == 'T' * (homo_length + 1) or str_tmp == 'C' * (
                    homo_length + 1) or str_tmp == 'G' * (homo_length + 1):
                return 1
        else:
            return 0

    def not_satisfy_uncorrelated_address(self, seq, index_length):
        address = seq[: index_length]
        if address in seq[index_length:]:
            return 1
        else:
            return 0

    def seq_xor(self, info_seq, state_seq, index_length):
        result_seq = []
        for i in range(index_length):
            result_seq.append(info_seq[i])
        for i in range(index_length, len(info_seq)):
            result_base = self.base_xor2(info_seq[i], state_seq[i])
            result_seq.append(result_base)
        result_seq = ''.join(result_seq)
        return result_seq

    def base_rotate(self, info_base, state_base):
        state_num = dic1.index(state_base)
        info_index = (dic1.index(info_base) + state_num) % 4
        new_base = dic1[info_index]
        return new_base

    def seq_rotate(self, info_seq, state_seq, index_length):
        result_seq = []
        # for i in range(index_length):
        #     result_seq.append(info_seq[i])
        for i in range(0, len(info_seq)):
            result_base = self.base_rotate(info_seq[i], state_seq[i])
            result_seq.append(result_base)
        result_seq = ''.join(result_seq)
        return result_seq

    def base_rerotate(self, info_base, state_base):
        state_num = dic1.index(state_base)
        info_index = (dic1.index(info_base) - state_num) % 4
        new_base = dic1[info_index]
        return new_base

    def seq_rerotate(self, info_seq, state_seq, index_length):
        result_seq = []
        # for i in range(index_length):
        #     result_seq.append(info_seq[i])
        for i in range(0, len(info_seq)):
            result_base = self.base_rerotate(info_seq[i], state_seq[i])
            result_seq.append(result_base)
        result_seq = ''.join(result_seq)
        return result_seq

    def mix_rule(self, state_seq):
        seq_90 = self.rule_90(state_seq)
        seq_150 = self.rule_150(state_seq)
        new_state = list(seq_90)
        for i in range(0, len(state_seq), 2):
            new_state[i] = seq_150[i]
        new_state = ''.join(new_state)
        return new_state

    def seq_encry(self, info_seq, state_seq, homolength, index_length):
        cycle_dic = []
        for item in itertools.product('GATC', repeat=3):
            cycle_dic.append(''.join(item))
        cycle_num = 0

        mid_state_seq = state_seq
        precycle_base = info_seq[0]
        precycle_index = dic1.index(precycle_base)
        # print(precycle_index)
        result_seq = info_seq

        for i in range(0, precycle_index + 1):
            new_state = self.mix_rule(mid_state_seq)
            # result_seq = self.seq_xor(result_seq, new_state, index_length)
            result_seq = self.seq_rotate(result_seq, new_state, index_length)
            mid_state_seq = new_state
            # print(result_seq)
            cycle_num += 1

        while (
                self.not_satisfy_homo(result_seq, homolength)
                or self.not_satisfy_gc(result_seq)
                or self.not_satify_random(result_seq)
                # or self.not_satisfy_uncorrelated_address(result_seq, index_length)
        ):
            new_state = self.mix_rule(mid_state_seq)
            # result_seq = self.seq_xor(result_seq, new_state, index_length)
            result_seq = self.seq_rotate(result_seq, new_state, index_length)
            mid_state_seq = new_state
            cycle_num += 1

        cycle_base = cycle_dic[cycle_num]
        # result_seq = cycle_base + result_seq
        return result_seq, cycle_base

    def seq_decry(self, result_seq, state_seq, index_length, cycle_str):
        # print(seq_num)
        mid_state_seq = state_seq
        cycle_dic = []
        for item in itertools.product('GATC', repeat=3):
            cycle_dic.append(''.join(item))
        cycle_num = cycle_dic.index(cycle_str)

        for i in range(0, cycle_num):
            new_state = self.mix_rule(mid_state_seq)
            # result_seq = self.seq_xor(result_seq, new_state, index_length)
            result_seq = self.seq_rerotate(result_seq, new_state, index_length)
            # print(result_seq)
            mid_state_seq = new_state
        return result_seq

    def ranseq_get(self, stateseq, seqnum, file_path):
        with open(file_path, 'w') as f:
            for i in range(0, seqnum):
                stateseq = self.mix_rule(stateseq)
                f.write(stateseq)


def caculate_rule150table(CA):
    list_base = ['A', 'T', 'C', 'G']
    list_data = []
    list_value = []
    for base1 in list_base:
        for base2 in list_base:
            for base3 in list_base:
                seq = base1 + base2 + base3
                list_data.append(seq)
                list_value.append(CA.base_xor3(base1, base2, base3))
    data = {
        "input": list_data,
        "output": list_value
    }

    pd.set_option("display.max_colwidth", None)  # 设置列宽为无限制
    pd.set_option("display.max_rows", None)  # 设置显示所有行
    pd.set_option("display.max_columns", None)  # 设置显示所有列
    pd.set_option("display.width", None)  # 设置终端宽度为无限制

    df = pd.DataFrame(data)
    df.to_csv("output.csv", index=False)

if __name__ == "__main__":
    CA = CellularAutomata()
    infoseq = "AATCGATTAAACCCGACGAATTACTATTAGATCAGTTTGGAGCTCTCCACCATGCTGTTATACCTCGAGAACGACTTCGCGAAGGAAGGCTTCGCCGTATCGCGGCTTAGACTCCCAGGCCTGGGTGGCGCCG"
    stateseq = "AAAAAAGCCGATTGACAGCGCGCGCACTGTGAGCCCCATCGGGGAAGTGTCCGCTCTCCGCGCGACTAACTACTACACAATAGTGACTAGCACCGAGCTCTAATTTCCCATTTGAGGCAGACTGTGCCTGGTG"
    current_path = r"C:\Users\Administrator\Desktop\papers\codes\Chamaeleo-master"
    generated_file_path = os.path.join(current_path, "examples", "generated_files")
    ranseq_path = os.path.join(generated_file_path, "ranseq", "ranseq.txt")

    # CA.ranseq_get(stateseq, 20000, ranseq_path)

    # encry_seq, cycle_str = CA.seq_encry(infoseq, stateseq, 5, 8)
    # decry_seq = CA.seq_decry(encry_seq, stateseq, 8, cycle_str)
    # print(decry_seq)
    # print(encry_seq[::-1])
    # print(eq(decry_seq, infoseq))
    caculate_rule150table(CA)
