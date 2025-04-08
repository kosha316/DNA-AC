import os
import Bio
from Bio.Seq import Seq
from complexcgr import FCGR
from collections import Counter
import math
import numpy as np
import cv2
from PIL import Image, ImageFile
from skimage.metrics import structural_similarity as ssim


def open_fasta(file_path):
    with open(file_path, "r", encoding='utf-8') as f:
        long_seq = f.read().replace('\n', '')
        # print(long_seq)
        f.close()
        return long_seq


def read_seq(file_path):
    with open(file_path, "r", encoding='utf-8') as f:
        fasta_segments = f.read()
        fasta_segments = fasta_segments.strip('\n')
        seq_list = fasta_segments.split('\n')
        f.close()
        return seq_list


def gc_compute(fasta_segments):
    gc_count = 0.00
    base_count = 0.00

    for base in fasta_segments:
        if base in {'A', 'T', 'C', 'G'}:
            base_count += 1.0
            if base in {'C', 'G'}:
                gc_count += 1.0
    return round(gc_count / base_count, 4)


def vitro_base_count(seq_list):
    base_num = 0
    for seq in seq_list:
        base_num += len(seq)
    return base_num


def vivo_base_count(seq_list, index_length):
    base_num = 0
    for seq in seq_list:
        base_num += (len(seq) - index_length)
    return base_num


def gc_seq_compute(fasta_segments):
    seq_gc_satis = 0.00

    fasta_segments = fasta_segments.strip('\n')
    seq_list = fasta_segments.split('\n')

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


def run_length_compute(fasta_segments):
    max_run_length = 0

    fasta_segments = fasta_segments.strip('\n')
    seq_list = fasta_segments.split('\n')
    for seq in seq_list:
        seq_run_length = compute_seq_runlength(seq)
        if seq_run_length > max_run_length:
            # print(seq)
            max_run_length = seq_run_length

    return max_run_length


def get_bits(file_path):
    with open(file_path, "rb") as f:
        file_bits = f.seek(0, os.SEEK_END) * 8
        f.close()
        return file_bits


def count_density(file_bit, fasta_path, index_length, ec_code):
    seq_list = read_seq(fasta_path)
    if ec_code == "ReedSolomon":
        ec_length = 16
        net_info_length = len(seq_list[0]) - index_length - ec_length
        net_base_num = net_info_length * (len(seq_list) + 1)
        net_info_dens = round(float(file_bit) / float(net_base_num), 4)

    return net_info_dens


def count_base_radio(plain_data, cipher_data):
    plain_base = vitro_base_count(plain_data)
    cipher_base = vitro_base_count(cipher_data)
    return round(float(plain_base) / float(cipher_base), 4)


def count_constrains_files(fasta_path):
    fasta_files = os.listdir(fasta_path)
    tplt = "{0:^20}\t{1:^20}\t{2:^20}\t{3:^20}\t{4:^20}\t{5:^20}"
    print(
        tplt.format("model", "file gc", "seq gc", "homo", "mean_tm", "variance_tm"))
    for file in fasta_files:
        file_path = fasta_path + "\\" + file
        f = open(file_path, 'r', encoding='utf-8')
        data = f.read()
        gc_percent = gc_compute(data)
        gc_seq = gc_seq_compute(data)
        homo_length = run_length_compute(data)
        mean_tm, variance_tm = calculate_tm_statistics(file_path)
        print(tplt.format(file.replace('.txt', ''), gc_percent, gc_seq, homo_length, mean_tm, variance_tm))
    return 0


def calculate_entropy(fasta_path):
    fasta_files = os.listdir(fasta_path)
    tplt = "{0:^20}\t{1:^20}\t{2:^20}\t{3:^20}\t{4:^20}\t{5:^20}"
    print(tplt.format("model", "A", "T", "C", "G", "entropy"))
    for file in fasta_files:
        file_path = fasta_path + "\\" + file
        f = open(file_path, 'r', encoding='utf-8')
        data = f.read()
        data = data.replace('\n', '')

        # 统计每个字符的出现次数
        char_counts = Counter(data)
        total_chars = len(data)

        # 计算每个字符的概率
        probabilities = {char: round(count / total_chars, 4) for char, count in char_counts.items()}
        per_a = probabilities['A']
        per_t = probabilities['T']
        per_c = probabilities['C']
        per_g = probabilities['G']

        # 计算信息熵
        entropy = round(-sum(p * math.log2(p) for p in probabilities.values() if p > 0), 4)
        print(file.replace('.txt', ''), probabilities, 'entropy:', entropy)
        # print(tplt.format(file.replace('.txt', ''), per_a, per_t, per_c, per_g, entropy))

    return 0


def chaos_games_repretation(file_path, img_path):
    fcgr = FCGR(k=7, bits=16)
    seq = ''.join(open_fasta(file_path))
    chaos = fcgr(seq)
    fcgr.plot(chaos)
    fcgr.save_img(chaos, path=img_path)
    return 0


def files_fcgr(fasta_path, img_path):
    fasta_files = os.listdir(fasta_path)
    for file in fasta_files:
        file_path = fasta_path + "\\" + file
        fcgr_path = img_path + "\\" + file + ".png"
        chaos_games_repretation(file_path, fcgr_path)
    return 0


def caculate_mse(img1, img2):
    img1 = img1.astype(np.float64)
    img2 = img2.astype(np.float64)
    mse = np.mean((img1 - img2) ** 2)
    return mse


def calculate_psnr(img1, img2):
    img1 = img1.astype(np.float64)
    img2 = img2.astype(np.float64)
    mse = np.mean((img1 - img2) ** 2)
    if mse == 0:
        return float('inf')
    return 20 * np.log10(255.0 / np.sqrt(mse))


def calculate_ssim(img1, img2, C1=(0.01 * 255) ** 2, C2=(0.03 * 255) ** 2):
    img1 = img1.astype(np.float64)
    img2 = img2.astype(np.float64)
    kernel = cv2.getGaussianKernel(11, 1.5)
    window = np.outer(kernel, kernel.transpose())

    mu1 = cv2.filter2D(img1, -1, window)[5:-5, 5:-5]  # 这里简化了边界处理
    mu2 = cv2.filter2D(img2, -1, window)[5:-5, 5:-5]
    mu1_sq = mu1 ** 2
    mu2_sq = mu2 ** 2
    mu1_mu2 = mu1 * mu2
    sigma1_sq = cv2.filter2D(img1 ** 2, -1, window)[5:-5, 5:-5] - mu1_sq
    sigma2_sq = cv2.filter2D(img2 ** 2, -1, window)[5:-5, 5:-5] - mu2_sq
    sigma12 = cv2.filter2D(img1 * img2, -1, window)[5:-5, 5:-5] - mu1_mu2

    ssim_map = ((2 * mu1_mu2 + C1) * (2 * sigma12 + C2)) / ((mu1_sq + mu2_sq + C1) *
                                                            (sigma1_sq + sigma2_sq + C2))
    return ssim_map.mean()


def calculate_psnr_ssim(oriimg_fath, recover_path):
    img_files = os.listdir(recover_path)

    tplt = "{0:^20}\t{1:^20}\t{2:^20}\t{3:^20}\t{4:^20}\t{5:^20}\t{6:^20}"
    print(tplt.format("error", "PSNR", "SSIM", 'MSE', 'FSIM', 'MS-SSIM', 'UQI'))
    ori_img = cv2.imread(oriimg_fath)

    for file in img_files:
        file_path = recover_path + "\\" + file
        recover_img = cv2.imread(file_path)
        psnr = calculate_psnr(ori_img, recover_img)
        ssim = calculate_ssim(ori_img, recover_img)
        mse = caculate_mse(ori_img, recover_img)
        fsim = calculate_fsim(ori_img, recover_img)
        msssim = caculate_msssim(ori_img, recover_img)
        uqi = calculate_uqi(ori_img, recover_img)
        print(tplt.format(file.replace('.png', ''), round(psnr, 4), round(ssim, 4), round(mse, 4),
                          round(fsim, 4), round(msssim, 4), round(uqi, 4)))

    return 0


def calculate_phase_congruency(image):
    """
    计算图像的相位一致性。
    """
    # 使用傅里叶变换计算相位信息
    F = np.fft.fft2(image)
    phase = np.angle(F)
    return phase


def calculate_gradient_magnitude(image):
    """
    计算图像的梯度幅度。
    """
    grad_x = cv2.Sobel(image, cv2.CV_64F, 1, 0, ksize=3)
    grad_y = cv2.Sobel(image, cv2.CV_64F, 0, 1, ksize=3)
    magnitude = np.sqrt(grad_x ** 2 + grad_y ** 2)
    return magnitude


def calculate_fsim(image1, image2):
    """
    计算两幅图像的FSIM值。
    """
    # 将图像转换为灰度图
    gray1 = cv2.cvtColor(image1, cv2.COLOR_BGR2GRAY).astype(np.float64)
    gray2 = cv2.cvtColor(image2, cv2.COLOR_BGR2GRAY).astype(np.float64)

    # 计算相位一致性
    phase1 = calculate_phase_congruency(gray1)
    phase2 = calculate_phase_congruency(gray2)
    phase_similarity = (2 * phase1 * phase2 + 1) / (phase1 ** 2 + phase2 ** 2 + 1)

    # 计算梯度幅度
    grad1 = calculate_gradient_magnitude(gray1)
    grad2 = calculate_gradient_magnitude(gray2)
    grad_similarity = (2 * grad1 * grad2 + 1) / (grad1 ** 2 + grad2 ** 2 + 1)

    # 计算FSIM
    fsim_map = phase_similarity * grad_similarity
    fsim_score = np.mean(fsim_map)
    return fsim_score


def gaussian_pyramid(image, levels):
    """构建高斯金字塔"""
    pyramid = [image]
    for _ in range(levels - 1):
        image = cv2.pyrDown(image)
        pyramid.append(image)
    return pyramid


def caculate_msssim(img1, img2, levels=5):
    """计算MS-SSIM"""
    pyramid1 = gaussian_pyramid(img1, levels)
    pyramid2 = gaussian_pyramid(img2, levels)

    win_size = min(img1.shape)  # 窗口大小取图像较小边的值
    if win_size % 2 == 0:  # 确保窗口大小为奇数
        win_size -= 1

    msssim_value = 1.0
    for i in range(levels):
        ssim_value = ssim(pyramid1[i], pyramid2[i], win_size=win_size, multichannel=True)
        msssim_value *= ssim_value

    return abs(msssim_value) ** (1.0 / levels)


def calculate_uqi(img1, img2):


    # 检查图像尺寸是否相同
    if img1.shape != img2.shape:
        raise ValueError("输入的两幅图像必须具有相同的尺寸和通道数")

    # 将图像转换为浮点型，并归一化到[0,1]范围
    img1 = img1.astype(np.float64) / 255.0
    img2 = img2.astype(np.float64) / 255.0

    # 计算平均值
    mu_x = np.mean(img1, axis=(0, 1))
    mu_y = np.mean(img2, axis=(0, 1))

    # 计算方差
    sigma_x_sq = np.var(img1, axis=(0, 1))
    sigma_y_sq = np.var(img2, axis=(0, 1))

    # 计算协方差
    sigma_xy = np.mean((img1 - mu_x) * (img2 - mu_y), axis=(0, 1))

    # 计算UQI每个通道
    numerator = 4 * sigma_xy * mu_x * mu_y
    denominator = (sigma_x_sq + sigma_y_sq) * (mu_x ** 2 + mu_y ** 2)

    # 防止除以零
    denominator = np.where(denominator == 0, 1e-10, denominator)

    uqi_channels = numerator / denominator

    # 取平均UQI
    uqi = np.mean(uqi_channels)

    return uqi


def fix_img(img_path, fix_img_path):
    ImageFile.LOAD_TRUNCATED_IMAGES = True
    img_files = os.listdir(img_path)
    for file in img_files:
        file_path = img_path + "\\" + file
        img = Image.open(file_path).convert("L")  # 转换为灰度图
        img.save(fix_img_path + "\\" + file)


def caculate_code_rate(file_path, fasta_path, key_path, index_length):
    file_bits = float(get_bits(file_path))
    seq_list = read_seq(fasta_path)
    key_seq = read_seq(key_path)
    vitro_base_num = float(vitro_base_count(seq_list))
    vivo_base_num = float(vivo_base_count(seq_list, index_length))
    key_num = float(vitro_base_count(key_seq))
    vitro_code_rate = round(file_bits / (vitro_base_num + key_num), 4)
    vivo_code_rate = round(file_bits / (vivo_base_num + key_num), 4)
    print('vitro_code_rate:', vitro_code_rate, 'vivo_code_rate:', vivo_code_rate)
    return vitro_code_rate, vivo_code_rate


def calculate_tm(sequence, method="nearest_neighbor"):
    """
    计算单条DNA序列的解链温度 (Tm)。
    支持两种方法：
    1. 简化公式（适用于短序列 <14nt）
    2. 最邻近法（适用于长序列 ≥14nt）
    """
    # 计算碱基数量
    A = sequence.count("A")
    T = sequence.count("T")
    G = sequence.count("G")
    C = sequence.count("C")
    length = len(sequence)

    if method == "simple" or length < 14:
        # 简化公式
        tm = 2 * (A + T) + 4 * (G + C)
    elif method == "nearest_neighbor":
        # 最邻近法（简化版）
        # 热力学参数（单位：kcal/mol）
        delta_H = {
            "AA": -7.9, "AT": -7.2, "AC": -8.4, "AG": -7.8,
            "TA": -7.2, "TT": -7.9, "TC": -8.2, "TG": -8.5,
            "CA": -8.5, "CT": -7.8, "CC": -8.0, "CG": -10.6,
            "GA": -8.2, "GT": -8.4, "GC": -9.8, "GG": -11.0
        }
        delta_S = {
            "AA": -22.2, "AT": -20.4, "AC": -22.4, "AG": -21.0,
            "TA": -21.3, "TT": -22.2, "TC": -22.2, "TG": -22.7,
            "CA": -22.7, "CT": -21.0, "CC": -19.9, "CG": -27.2,
            "GA": -20.4, "GT": -22.4, "GC": -24.4, "GG": -26.7
        }
        # 计算焓变和熵变
        total_delta_H = 0
        total_delta_S = 0
        for i in range(length - 1):
            pair = sequence[i:i + 2]
            if pair[0] == 'Z':
                pair = 'A' + pair[1]
            if pair[1] == 'Z':
                pair = pair[0] + 'A'
            total_delta_H += delta_H[pair]
            total_delta_S += delta_S[pair]

        # 计算Tm
        tm = (total_delta_H * 1000) / (total_delta_S + 1.987 * np.log(0.5)) - 273.15
    else:
        raise ValueError("Invalid method. Choose 'simple' or 'nearest_neighbor'.")

    return tm


def calculate_tm_statistics(file_path, method="nearest_neighbor"):
    """
    计算多条DNA序列的Tm均值和方差。
    """
    sequences = read_seq(file_path)
    tm_values = [calculate_tm(seq, method) for seq in sequences]
    mean_tm = round(np.mean(tm_values), 4)
    variance_tm = round(np.var(tm_values), 4)
    return mean_tm, variance_tm


if __name__ == "__main__":
    current_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    oriimg_file_path = os.path.join(current_path, "secure", "input", "x_ray.png")
    plain_path = os.path.join(current_path, "secure", "input", "x_ray_base.txt")
    fasta_path = os.path.join(current_path, "secure", "biotest", 'rate_files', 'info_file', 'ARC4.txt')
    fre_path = os.path.join(current_path, "secure", "biotest", 'fre_files')
    key_path = os.path.join(current_path, "secure", "biotest", 'rate_files', 'key_file', 'ARC4.key')
    recovery_img_path = os.path.join(current_path, "secure", "biotest", "images")
    fixed_img_path = os.path.join(current_path, "secure", "biotest", "recon_imgs")

    # file_bits = get_bits(oriimg_file_path)
    # count_constrains_files(fre_path)

    calculate_entropy(fre_path)
    # caculate_code_rate(oriimg_file_path, fasta_path, key_path, index_length=8)
    # fix_img(recovery_img_path, fixed_img_path)
    # calculate_psnr_ssim(oriimg_file_path, fixed_img_path)
    # files_fcgr(fasta_files, fasta_path, img_path)
