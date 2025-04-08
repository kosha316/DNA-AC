import random
import editdistance
import numpy as np
from PIL import Image
import itertools
import trans_type
import os


def ran1bit(plian_bytes):
    bytelen = len(plian_bytes)
    bytes_array = bytearray(plian_bytes)
    byte_position = random.randint(0, bytelen)
    byte = plian_bytes[byte_position]
    print(byte)
    bit_position = random.randint(0, 8)
    byte ^= (1 << bit_position)
    byte = byte % 256
    bytes_array[byte_position] = byte
    return bytes(bytes_array)


def get_codebook():
    codebook = []
    for item in itertools.product('ATCG', repeat=4):
        codebook.append(''.join(item))
    return codebook
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


def add_index(seq_list, index_len):
    outlist = []
    list_len = len(seq_list)
    # print(index_len)
    dic_index = []
    for item in itertools.product('ATCG', repeat=index_len):
        dic_index.append(''.join(item))
    for i in range(0, list_len):
        newseq = dic_index[i] + seq_list[i]
        # print(newseq)
        outlist.append(newseq)
    return outlist


def remove_index(seq_list, index_len):
    outlist = []
    dic_index = []
    list_len = len(seq_list)
    for item in itertools.product('ATCG', repeat=index_len):
        dic_index.append(''.join(item))
    for i in range(0, list_len):
        index_str = dic_index[i]
        matched_seq = next((s for s in seq_list if s.startswith(index_str)), None)
        # print(matched_seq)
        if matched_seq is None:
            matched_seq = find_most_similar(index_str, seq_list)
        outlist.append(matched_seq[index_len:])
    return outlist


def seq_pad(seq_list, seq_length, pad_str):
    endseq = seq_list[len(seq_list) - 1]
    while len(endseq) < seq_length:
        endseq += pad_str
    seq_list[len(seq_list) - 1] = endseq
    return seq_list


def remove_pad(seq_list, pad_str):
    s = seq_list[len(seq_list) - 1]
    while s.endswith(pad_str) and s != pad_str:
        s = s[:-len(pad_str)]
    seq_list[len(seq_list) - 1] = s
    return seq_list


def encode(file_path, seq_length, index_length, pad_str):
    with open(file_path, 'rb') as f:
        data = f.read()
    #改变一比特雪崩测试
    # data = ran1bit(data)
    seq_list = trans_type.bytes_to_fasta(data, seq_length)
    seq_list = seq_pad(seq_list, seq_length, pad_str)
    # print(seq_list)
    seq_list = add_index(seq_list, index_length)
    return seq_list


def decode(seq_list, index_length, pad_str, file_path):
    # print(seq_list)
    seq_list = remove_index(seq_list, index_length)
    # print(seq_list)
    seq_list = remove_pad(seq_list, pad_str)
    # print(seq_list, len(seq_list))
    data = trans_type.fasta_to_bytes(seq_list)
    with open(file_path, 'wb') as f:
        f.write(data)
    # write_seq(file_path, seq_list)


def img_encode(file_path, seq_length, index_length, pad_str):
    # 构建码本
    codebook = get_codebook()
    list_base = []
    seq_list = []
    # 打开图像
    image = Image.open(file_path)  # 替换为你的图像路径

    # 获取图像的尺寸
    width, height = image.size
    print(f"Image dimensions: Width = {width}, Height = {height}")

    # 遍历图像的所有像素值
    for y in range(height):
        for x in range(width):
            pixel = image.getpixel((x, y))  # 获取像素值
            base = codebook[pixel]
            # 得到对应DNA序列
            list_base.append(base)
    #
    for i in range(0, len(list_base), int(seq_length/4)):
        seq = ''.join(list_base[i:i+int(seq_length/4)])
        seq_list.append(seq)

    seq_list = seq_pad(seq_list, seq_length, pad_str)
    seq_list = add_index(seq_list, index_length)
    print(len(seq_list), len(seq_list[0]))
    return seq_list


def img_decode(seq_list, index_length, pad_str, file_path, width, height):
    # 构建码本
    codebook = get_codebook()
    list_pixel = []
    seq_list = remove_index(seq_list, index_length)
    # print(seq_list)
    seq_list = remove_pad(seq_list, pad_str)
    for seq in seq_list:
        for i in range(0, len(seq), 4):
            base = ''.join(seq[i:i+4])
            pixel = codebook.index(base)
            list_pixel.append(pixel)
    pixel_array = np.array(list_pixel)
    # print(pixel_array)
    pixel_reshape = np.reshape(pixel_array, (width, height))
    # print(pixel_reshape)
    # # 使用Pillow将像素值矩阵转换为图像
    pixel_reshape = pixel_reshape.astype(np.uint8)
    image = Image.fromarray(pixel_reshape)
    image.save(file_path)

if __name__ == "__main__":
    current_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "secure")
    file_path = os.path.join(current_path, "input", "x_ray.png")
    encode_path = os.path.join(current_path, "input", "x_ray_base.txt")
    decode_path = os.path.join(current_path, "input", "decode.png")

    seq_list = img_encode(file_path, 128, 8, 'GCAT')
    img_decode(seq_list, 8, 'GCAT', decode_path, 256, 256)
    # seq_list = encode(file_path, 128, 8, 'ACGT')
    # print(seq_list)
    # write_seq(encode_path, seq_list)
    # seq_list = read_seq(encode_path)
    # decode(seq_list, 8, 'ACGT', decode_path)
