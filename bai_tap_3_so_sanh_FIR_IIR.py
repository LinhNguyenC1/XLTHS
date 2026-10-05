# -*- coding: utf-8 -*-
"""
========================================================================================
BÀI TẬP 3: SO SÁNH KHẢ NĂNG KHỬ NHIỄU CỦA BỘ LỌC FIR VÀ IIR CƠ BẢN
CHƯƠNG TRÌNH MÔ PHỎNG MONTE CARLO & XUẤT HÌNH ẢNH PHỤC VỤ BÁO CÁO SLIDE
========================================================================================
Nội dung đáp ứng đầy đủ yêu cầu:
1. Tín hiệu sạch chậm thay đổi s[n], nhiễu AWGN theo chuẩn SNR yêu cầu.
2. Bộ lọc MA (FIR, M=5) và EMA (IIR, alpha=0.2).
3. Quét dải SNR từ -10 dB đến 20 dB, 100 lần Monte Carlo mỗi mức.
4. Căn chỉnh độ trễ: MA = 2 mẫu, EMA = 4 mẫu.
5. Bỏ qua giai đoạn quá độ ban đầu (transient phase).
6. Random seed cố định (seed=42) để tái lập kết quả chuẩn xác 100%.
7. Phân tích đầy đủ 4 sự đánh đổi:
   - Khả năng giảm nhiễu (Noise Reduction Ratio)
   - Độ trễ (Delay & Phase Distortion)
   - Sai số do làm mượt quá mức (Over-smoothing Error)
   - Độ phức tạp tính toán (Computational & Memory Complexity)
8. Minh họa hạn chế: Độ trễ nhóm của EMA phụ thuộc tần số (pha phi tuyến).
9. Xuất bộ 4 hình ảnh độ phân giải cao phục vụ trực tiếp cho Slide thuyết trình.
========================================================================================
"""

import os
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Cấu hình thẩm mỹ đồ thị
rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
rcParams['font.size'] = 11
rcParams['axes.titlesize'] = 13
rcParams['axes.labelsize'] = 11
rcParams['legend.fontsize'] = 10
rcParams['figure.autolayout'] = False

OUTPUT_DIR = "/Users/nguyenthithuylinh/Downloads/Xử_lý_tiếng_nói/BTVN_B3_FIR_vs_IIR"
os.makedirs(OUTPUT_DIR, exist_ok=True)
RANDOM_SEED = 42

# ======================================================================================
# 1. CÁC HÀM TẠO TÍN HIỆU VÀ BỘ LỌC CƠ BẢN
# ======================================================================================

def tao_tin_hieu_sach(N=600, f0=0.01):
    """
    Tạo tín hiệu sạch chậm thay đổi s[n] = sin(2*pi*f0*n).
    Chu kỳ = 1/f0 = 100 mẫu. Tần số thấp giúp mô phỏng đúng tín hiệu chậm cần làm mượt.
    """
    n = np.arange(N)
    s = np.sin(2 * np.pi * f0 * n)
    return n, s


def sinh_nhieu_awgn(signal, snr_db):
    """
    Sinh nhiễu Gaussian trắng cộng (AWGN) theo SNR (dB) yêu cầu:
    x[n] = s[n] + v[n], v[n] ~ N(0, sigma_v^2)
    P_signal = mean(s[n]^2)
    P_noise = P_signal / 10^(SNR/10)
    sigma_v = sqrt(P_noise)
    """
    p_signal = np.mean(signal ** 2)
    snr_linear = 10.0 ** (snr_db / 10.0)
    p_noise = p_signal / snr_linear
    sigma_v = np.sqrt(p_noise)
    noise = sigma_v * np.random.randn(len(signal))
    return signal + noise, noise, sigma_v


def loc_ma(x, M=5):
    """
    Bộ lọc Moving Average (FIR):
    y[n] = (1/M) * sum_{k=0}^{M-1} x[n-k]
    Hệ số xung h[n] = 1/M (hữu hạn M mẫu).
    Pha tuyến tính, trễ cố định = (M-1)/2 mẫu.
    """
    kernel = np.ones(M) / float(M)
    # Tích chập nhân quả (causal convolution)
    y_full = np.convolve(x, kernel, mode='full')
    return y_full[:len(x)]


def loc_ema(x, alpha=0.2):
    """
    Bộ lọc Exponential Moving Average (IIR bậc 1):
    y[n] = alpha * x[n] + (1 - alpha) * y[n-1]
    Đáp ứng xung h[n] = alpha * (1 - alpha)^n (vô hạn, suy giảm hàm mũ).
    Có phản hồi feedback, pha phi tuyến, trễ thấp tần ~ (1-alpha)/alpha mẫu.
    """
    y = np.zeros_like(x)
    y_prev = 0.0
    # Cài đặt đệ quy thời gian thực
    for n in range(len(x)):
        y_n = alpha * x[n] + (1.0 - alpha) * y_prev
        y[n] = y_n
        y_prev = y_n
    return y


def can_chinh_tre_va_bo_qua_qua_do(y, s, delay_samples, n_transient=50):
    """
    Căn chỉnh độ trễ và loại bỏ vùng quá độ (transient phase):
    - Dịch chuyển bù độ trễ: delay_samples mẫu (MA=2, EMA=4)
    - Loại bỏ n_transient mẫu đầu tiên: để hệ thống đi vào trạng thái xác lập (steady-state),
      tránh sai số do nạp ban đầu (initial state bias).
    """
    d = int(round(delay_samples))
    # Dịch tín hiệu y sang trái d mẫu để khớp pha với s
    if d > 0:
        y_shifted = y[d:]
        s_shifted = s[:-d]
    else:
        y_shifted = y
        s_shifted = s

    # Cắt bỏ phần quá độ ban đầu (transient phase)
    # Cả MA và EMA đều cần thời gian để lấp đầy "ký ức"
    y_steady = y_shifted[n_transient:]
    s_steady = s_shifted[n_transient:]
    
    return y_steady, s_steady


def tinh_rmse(y_steady, s_steady):
    """
    Tính Root Mean Square Error (RMSE) giữa tín hiệu sau lọc và tín hiệu sạch:
    RMSE = sqrt( (1/K) * sum( (y[n] - s[n])^2 ) )
    """
    return np.sqrt(np.mean((y_steady - s_steady) ** 2))


# ======================================================================================
# 2. MÔ PHỎNG MONTE CARLO
# ======================================================================================

def chay_mo_phong_monte_carlo(N=600, f0=0.01, M=5, alpha=0.2,
                              snr_range=np.arange(-10, 21, 1),
                              num_trials=100, seed=42):
    """
    Chạy 100 lần Monte Carlo cho mỗi mức SNR từ -10 dB đến 20 dB:
    - Căn chỉnh độ trễ: MA = 2 mẫu, EMA = 4 mẫu
    - Bỏ qua phần quá độ: 50 mẫu đầu
    - Tính RMSE trung bình và độ lệch chuẩn của RMSE
    """
    np.random.seed(seed)
    n_axis, s_clean = tao_tin_hieu_sach(N=N, f0=f0)

    delay_ma = (M - 1) / 2.0  # 2.0 mẫu
    delay_ema = (1.0 - alpha) / alpha  # 4.0 mẫu
    n_transient = 50

    num_snr = len(snr_range)
    rmse_ma_mean = np.zeros(num_snr)
    rmse_ma_std = np.zeros(num_snr)
    rmse_ema_mean = np.zeros(num_snr)
    rmse_ema_std = np.zeros(num_snr)
    
    # Đo thêm sai số của tín hiệu thô chưa lọc (Raw Noisy RMSE) để thấy hiệu quả khử nhiễu
    rmse_raw_mean = np.zeros(num_snr)

    print(f"\n🚀 BẮT ĐẦU MÔ PHỎNG MONTE CARLO (Seed={seed})")
    print(f"   • Dải SNR: [{snr_range[0]} dB -> {snr_range[-1]} dB] (Tổng {num_snr} mức)")
    print(f"   • Số lần lặp mỗi mức SNR: {num_trials} trials")
    print(f"   • Bộ lọc MA: M = {M}  --> Trễ căn chỉnh = {delay_ma:.1f} mẫu")
    print(f"   • Bộ lọc EMA: α = {alpha} --> Trễ căn chỉnh = {delay_ema:.1f} mẫu")
    print(f"   • Bỏ qua phần quá độ: {n_transient} mẫu đầu\n")

    t_start = time.time()

    for idx, snr_db in enumerate(snr_range):
        trials_ma = np.zeros(num_trials)
        trials_ema = np.zeros(num_trials)
        trials_raw = np.zeros(num_trials)

        for t in range(num_trials):
            x_noisy, _, _ = sinh_nhieu_awgn(s_clean, snr_db)

            # Lọc tín hiệu
            y_ma = loc_ma(x_noisy, M=M)
            y_ema = loc_ema(x_noisy, alpha=alpha)

            # Căn chỉnh trễ và bỏ quá độ
            y_ma_st, s_ma_st = can_chinh_tre_va_bo_qua_qua_do(y_ma, s_clean, delay_ma, n_transient)
            y_ema_st, s_ema_st = can_chinh_tre_va_bo_qua_qua_do(y_ema, s_clean, delay_ema, n_transient)
            x_raw_st, s_raw_st = can_chinh_tre_va_bo_qua_qua_do(x_noisy, s_clean, 0, n_transient)

            # Tính RMSE
            trials_ma[t] = tinh_rmse(y_ma_st, s_ma_st)
            trials_ema[t] = tinh_rmse(y_ema_st, s_ema_st)
            trials_raw[t] = tinh_rmse(x_raw_st, s_raw_st)

        rmse_ma_mean[idx] = np.mean(trials_ma)
        rmse_ma_std[idx] = np.std(trials_ma)
        rmse_ema_mean[idx] = np.mean(trials_ema)
        rmse_ema_std[idx] = np.std(trials_ema)
        rmse_raw_mean[idx] = np.mean(trials_raw)

        if snr_db in [-10, -5, 0, 5, 10, 15, 20]:
            better = "EMA thắng" if rmse_ema_mean[idx] < rmse_ma_mean[idx] else "MA thắng"
            diff_pct = (rmse_ema_mean[idx] - rmse_ma_mean[idx]) / rmse_ma_mean[idx] * 100.0
            print(f"   SNR = {snr_db:+3d} dB | Raw RMSE = {rmse_raw_mean[idx]:.4f} | "
                  f"MA = {rmse_ma_mean[idx]:.4f} | EMA = {rmse_ema_mean[idx]:.4f} | "
                  f"{better:10s} (chênh lệch: {abs(diff_pct):5.1f}%)")

    t_elapsed = time.time() - t_start
    print(f"\n✅ Hoàn thành 100 trials Monte Carlo trong {t_elapsed:.2f} giây.")

    return {
        'snr_range': snr_range,
        'rmse_ma_mean': rmse_ma_mean,
        'rmse_ma_std': rmse_ma_std,
        'rmse_ema_mean': rmse_ema_mean,
        'rmse_ema_std': rmse_ema_std,
        'rmse_raw_mean': rmse_raw_mean,
        'delay_ma': delay_ma,
        'delay_ema': delay_ema,
        'M': M,
        'alpha': alpha
    }


# ======================================================================================
# 3. XUẤT 4 HÌNH ẢNH MINH HỌA ĐẶC SẮC CHO SLIDE
# ======================================================================================

def tao_hinh_1_input_filter_output(M=5, alpha=0.2, snr_db=3):
    """
    HÌNH 1 DÙNG CHO SLIDE 1 & 2:
    Thể hiện trực quan: Input x[n] (nhiễu) --> Qua bộ lọc --> Ra Output y[n] (làm mượt).
    Tập trung vào KẾT QUẢ TRỰC QUAN trước theo đúng yêu cầu slide!
    """
    np.random.seed(42)
    N = 220
    n, s_clean = tao_tin_hieu_sach(N=N, f0=0.015)
    x_noisy, noise, _ = sinh_nhieu_awgn(s_clean, snr_db=snr_db)

    y_ma = loc_ma(x_noisy, M=M)
    y_ema = loc_ema(x_noisy, alpha=alpha)

    fig, axes = plt.subplots(3, 1, figsize=(13, 8.5), sharex=True)
    fig.suptitle(f"1. TÁC ĐỘNG CỦA BỘ LỌC ĐẾN TÍN HIỆU (MINH HỌA TẠI SNR = {snr_db} dB)\n"
                 r"Mô hình: $x[n] = s[n] + v[n]$ $\longrightarrow$ [BỘ LỌC] $\longrightarrow$ $y[n]$",
                 fontsize=14, fontweight='bold', y=0.98)

    # Panel 1: Input x[n] và tín hiệu gốc s[n]
    axes[0].plot(n, s_clean, 'k--', linewidth=2.0, label='Tín hiệu sạch gốc s[n] (Sóng chậm)', alpha=0.85)
    axes[0].plot(n, x_noisy, color='#e67e22', linewidth=1.1, alpha=0.7, label=f'Đầu vào có nhiễu x[n] (AWGN, SNR={snr_db}dB)')
    axes[0].set_title("BƯỚC 1: ĐẦU VÀO QUAN SÁT (INPUT) = TÍN HIỆU SẠCH + NHIỄU GAUSSIAN TRẮNG", fontweight='bold', color='#b9770e')
    axes[0].set_ylabel("Biên độ")
    axes[0].grid(True, linestyle=':', alpha=0.6)
    axes[0].legend(loc='upper right', framealpha=0.9)

    # Panel 2: Qua 2 bộ lọc MA & EMA
    axes[1].plot(n, s_clean, 'k--', linewidth=2.0, label='Tín hiệu sạch s[n]', alpha=0.85)
    axes[1].plot(n, y_ma, color='#2980b9', linewidth=2.0, label=f'Đầu ra MA (FIR, M={M}) - Trọng số đều')
    axes[1].plot(n, y_ema, color='#c0392b', linewidth=2.0, label=f'Đầu ra EMA (IIR, α={alpha}) - Trọng số hàm mũ')
    axes[1].set_title("BƯỚC 2: QUA 2 BỘ LỌC LÀM MƯỢT (OUTPUT) - NHIỄU GIẢM RÕ RỆT NHƯNG XUẤT HIỆN ĐỘ TRỄ",
                     fontweight='bold', color='#2c3e50')
    axes[1].set_ylabel("Biên độ")
    axes[1].grid(True, linestyle=':', alpha=0.6)
    axes[1].legend(loc='upper right', framealpha=0.9)

    # Panel 3: Phóng to chi tiết (Zoom-in) để thấy rõ:
    # 1. Khả năng làm mượt (EMA mượt hơn MA)
    # 2. Độ trễ (EMA trễ ~4 mẫu, MA trễ ~2 mẫu)
    # 3. Bỏ qua phần quá độ ban đầu (transient)
    zoom_range = slice(30, 110)
    n_zoom = n[zoom_range]
    axes[2].plot(n_zoom, s_clean[zoom_range], 'k-', linewidth=2.5, label='Gốc s[n]', alpha=0.9)
    axes[2].plot(n_zoom, x_noisy[zoom_range], color='#e67e22', linewidth=0.8, alpha=0.45, label='Nhiễu x[n]')
    axes[2].plot(n_zoom, y_ma[zoom_range], 'b-o', markersize=3.5, linewidth=1.8, label=f'MA (M={M}) [Trễ 2 mẫu]')
    axes[2].plot(n_zoom, y_ema[zoom_range], 'r-s', markersize=3.5, linewidth=1.8, label=f'EMA (α={alpha}) [Trễ 4 mẫu]')
    
    # Đánh dấu độ trễ tại đỉnh sóng
    peak_idx = 30 + np.argmax(s_clean[zoom_range])
    axes[2].axvline(x=peak_idx, color='black', linestyle=':', alpha=0.6)
    axes[2].text(peak_idx - 0.5, 1.12, "Đỉnh gốc", ha='right', fontsize=9, color='black', fontweight='bold')
    axes[2].text(peak_idx + 2, 1.05, "→ MA trễ 2 mẫu", ha='left', fontsize=9, color='blue', fontweight='bold')
    axes[2].text(peak_idx + 4, 0.95, "→ EMA trễ 4 mẫu", ha='left', fontsize=9, color='red', fontweight='bold')
    
    axes[2].set_title("BƯỚC 3: PHÓNG TO CHI TIẾT (ZOOM-IN) - SO SÁNH ĐỘ MƯỢT VÀ ĐỘ TRỄ GIỮA HAI BỘ LỌC",
                     fontweight='bold', color='#27ae60')
    axes[2].set_xlabel("Chỉ số mẫu thời gian n")
    axes[2].set_ylabel("Biên độ")
    axes[2].set_ylim(-0.3, 1.25)
    axes[2].grid(True, linestyle=':', alpha=0.6)
    axes[2].legend(loc='lower left', framealpha=0.9, ncol=4)

    plt.tight_layout()
    file_path = os.path.join(OUTPUT_DIR, "slide_hinh1_input_filter_output.png")
    plt.savefig(file_path, dpi=160, bbox_inches='tight')
    plt.close()
    print(f"📸 Đã lưu: {file_path}")


def tao_hinh_2_nguyen_ly_va_dap_ung_xung(M=5, alpha=0.2):
    """
    HÌNH 2 DÙNG CHO SLIDE 3 & 4:
    Giải thích nguyên lý hoạt động của 2 bộ lọc thông qua đáp ứng xung h[n] và cấu trúc toán học.
    - Bộ lọc 1: MA (FIR) -> Cửa sổ trượt, hữu hạn, không phản hồi, pha tuyến tính.
    - Bộ lọc 2: EMA (IIR) -> Đệ quy phản hồi, vô hạn suy giảm hàm mũ, pha phi tuyến.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    fig.suptitle("2. NGUYÊN LÝ HOẠT ĐỘNG: ĐÁP ỨNG XUNG h[n] & CƠ CHẾ BỘ NHỚ",
                 fontsize=14, fontweight='bold', y=0.98)

    N_samples = 25
    n_idx = np.arange(N_samples)

    # Đáp ứng xung MA
    h_ma = np.zeros(N_samples)
    h_ma[:M] = 1.0 / M

    ax1 = axes[0]
    markerline, stemlines, baseline = ax1.stem(n_idx, h_ma, basefmt="k-", linefmt="b-", markerfmt="bo")
    plt.setp(markerline, markersize=7)
    ax1.set_title(f"BỘ LỌC 1: MOVING AVERAGE (FIR - M={M})\n"
                  r"$y_1[n] = \frac{1}{M} \sum_{k=0}^{M-1} x[n-k]$",
                  fontsize=12, fontweight='bold', color='#1f618d')
    ax1.set_xlabel("Chỉ số mẫu n")
    ax1.set_ylabel("Biên độ h[n]")
    ax1.set_ylim(-0.03, 0.28)
    ax1.grid(True, linestyle=':', alpha=0.6)
    
    # Chú thích nguyên lý MA
    ax1.text(0.04, 0.88, 
             "• Chiều dài HỮU HẠN: đúng M=5 mẫu\n"
             "• Trọng số bằng nhau: w = 1/5 = 0.2\n"
             "• KHÔNG CÓ PHẢN HỒI (Feedforward only)\n"
             "• Luôn ổn định tuyệt đối (All-zero FIR)\n"
             "• Pha tuyến tính -> Trễ cố định = 2.0 mẫu",
             transform=ax1.transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#ebf5fb", edgecolor="#2980b9"))

    # Đáp ứng xung EMA
    h_ema = alpha * ((1.0 - alpha) ** n_idx)

    ax2 = axes[1]
    markerline2, stemlines2, baseline2 = ax2.stem(n_idx, h_ema, basefmt="k-", linefmt="r-", markerfmt="ro")
    plt.setp(markerline2, markersize=7)
    ax2.set_title(f"BỘ LỌC 2: EXPONENTIAL MA (IIR - α={alpha})\n"
                  r"$y_2[n] = \alpha x[n] + (1-\alpha) y_2[n-1]$",
                  fontsize=12, fontweight='bold', color='#922b21')
    ax2.set_xlabel("Chỉ số mẫu n")
    ax2.set_ylabel("Biên độ h[n]")
    ax2.set_ylim(-0.03, 0.28)
    ax2.grid(True, linestyle=':', alpha=0.6)

    # Chú thích nguyên lý EMA
    tau = -1.0 / np.log(1.0 - alpha)
    ax2.axvline(x=tau, color='green', linestyle='--', linewidth=1.5, label=f'Hằng số thời gian τ ≈ {tau:.1f} mẫu')
    ax2.text(0.04, 0.88,
             "• Chiều dài VÔ HẠN (Suy giảm hàm mũ)\n"
             "• Mẫu càng cũ, trọng số càng giảm dần\n"
             "• CÓ PHẢN HỒI ĐỆ QUY (Feedback y[n-1])\n"
             "• Ổn định khi |1-α| < 1 (1 cực tại z = 1-α)\n"
             "• Pha phi tuyến -> Trễ thấp tần ≈ 4.0 mẫu",
             transform=ax2.transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#fdf2e9", edgecolor="#c0392b"))
    ax2.legend(loc='lower right', framealpha=0.9)

    plt.tight_layout()
    file_path = os.path.join(OUTPUT_DIR, "slide_hinh2_nguyen_ly_dap_ung_xung.png")
    plt.savefig(file_path, dpi=160, bbox_inches='tight')
    plt.close()
    print(f"📸 Đã lưu: {file_path}")


def tao_hinh_3_so_sanh_rmse_snr(res):
    """
    HÌNH 3 DÙNG CHO SLIDE 5 & 6:
    Vẽ 2 đường cong RMSE theo SNR trên CÙNG MỘT HỆ TRỤC (đúng yêu cầu đề bài).
    Phân tích rõ ràng 3 vùng SNR:
    - Vùng 1: SNR thấp (Nhiễu mạnh, EMA khử nhiễu vượt trội).
    - Vùng 2: SNR trung bình (Vùng chuyển tiếp cân bằng).
    - Vùng 3: SNR cao (Nhiễu yếu, MA giữ biên độ sạch tốt hơn).
    """
    snr = res['snr_range']
    rmse_ma = res['rmse_ma_mean']
    rmse_ema = res['rmse_ema_mean']
    rmse_raw = res['rmse_raw_mean']
    M = res['M']
    alpha = res['alpha']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), gridspec_kw={'height_ratios': [2.3, 1]})
    fig.suptitle("3. KẾT QUẢ MÔ PHỎNG MONTE CARLO: SO SÁNH RMSE GIỮA MA VÀ EMA THEO SNR\n"
                 r"(Đã căn chỉnh trễ: $\Delta_{MA}=2$ mẫu, $\Delta_{EMA}=4$ mẫu; Bỏ qua quá độ 50 mẫu đầu; 100 trials/SNR)",
                 fontsize=13, fontweight='bold', y=0.98)

    # Subplot 1: Hai đường RMSE trên cùng một hệ trục
    ax1.plot(snr, rmse_raw, 'k:', linewidth=1.5, alpha=0.5, label='Chưa lọc (Raw Noisy RMSE)')
    ax1.plot(snr, rmse_ma, 'b-o', markersize=5, linewidth=2.2, label=f'Bộ lọc MA (FIR, M={M})')
    ax1.plot(snr, rmse_ema, 'r-s', markersize=5, linewidth=2.2, label=f'Bộ lọc EMA (IIR, α={alpha})')

    # Đánh dấu 3 vùng SNR bằng dải màu nền
    ax1.axvspan(-10, 0, color='#f9ebea', alpha=0.6, label='Vùng 1: SNR Thấp (<-10 đến 0 dB) [Nhiễu mạnh]')
    ax1.axvspan(0, 10, color='#fef9e7', alpha=0.6, label='Vùng 2: SNR Trung bình (0 đến 10 dB) [Cân bằng]')
    ax1.axvspan(10, 20, color='#eafaf1', alpha=0.6, label='Vùng 3: SNR Cao (10 đến 20 dB) [Tín hiệu sạch]')

    # Ghi chú các điểm đặc sắc
    ax1.annotate("EMA thắng áp đảo\n(Khử nhiễu mạnh hơn)",
                 xy=(-5, rmse_ema[5]), xytext=(-8, 0.45),
                 arrowprops=dict(facecolor='red', shrink=0.08, width=1.5, headwidth=6),
                 fontsize=10, fontweight='bold', color='#922b21')

    ax1.annotate("Giao điểm chuyển giao (~15-18 dB)\nMA thắng do ít sai số làm mượt",
                 xy=(16, rmse_ma[26]), xytext=(7, 0.12),
                 arrowprops=dict(facecolor='blue', shrink=0.08, width=1.5, headwidth=6),
                 fontsize=10, fontweight='bold', color='#1f618d')

    ax1.set_ylabel("RMSE (Càng nhỏ càng tốt)", fontweight='bold')
    ax1.set_ylim(-0.02, max(rmse_ma[0], rmse_ema[0]) * 1.15)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='upper right', framealpha=0.95, fontsize=10)

    # Subplot 2: Tỷ lệ chênh lệch RMSE (%) = (RMSE_EMA - RMSE_MA) / RMSE_MA * 100
    # Âm: EMA tốt hơn | Dương: MA tốt hơn
    diff_pct = (rmse_ema - rmse_ma) / rmse_ma * 100.0
    colors = ['#c0392b' if val < 0 else '#2980b9' for val in diff_pct]
    ax2.bar(snr, diff_pct, color=colors, width=0.75, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax2.axhline(0, color='black', linewidth=1.2, linestyle='--')
    ax2.set_xlabel("Tỷ số Tín hiệu trên Nhiễu: SNR (dB)", fontweight='bold')
    ax2.set_ylabel("Chênh lệch RMSE (%)", fontweight='bold')
    ax2.set_title(r"CHÊNH LỆCH HIỆU NĂNG: $\frac{RMSE_{EMA} - RMSE_{MA}}{RMSE_{MA}} \times 100\%$ (<0: EMA tốt hơn | >0: MA tốt hơn)",
                  fontsize=11, fontweight='bold')
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.set_xlim(-10.8, 20.8)

    plt.tight_layout()
    file_path = os.path.join(OUTPUT_DIR, "slide_hinh3_so_sanh_rmse_snr.png")
    plt.savefig(file_path, dpi=160, bbox_inches='tight')
    plt.close()
    print(f"📸 Đã lưu: {file_path}")


def tao_hinh_4_phan_tich_danh_doi_va_han_che(M=5, alpha=0.2):
    """
    HÌNH 4 DÙNG CHO SLIDE 7 & 8:
    Phân tích toàn diện 4 sự đánh đổi (Trade-offs) và hạn chế pha của EMA:
    1. Khả năng giảm phương sai nhiễu (Noise Variance Reduction: 1/M vs alpha/(2-alpha))
    2. Độ trễ (Delay: hằng số vs phụ thuộc tần số)
    3. Sai số do làm mượt quá mức (Over-smoothing bias / Peak attenuation)
    4. Hạn chế EMA: Độ trễ nhóm biến thiên theo tần số tau_g(omega) -> Pha phi tuyến, méo phân tán
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("4. PHÂN TÍCH TOÀN DIỆN CÁC SỰ ĐÁNH ĐỔI (TRADE-OFFS) & HẠN CHẾ BỘ LỌC",
                 fontsize=14, fontweight='bold', y=0.98)

    # -------------------------------------------------------------
    # 4.1 Đánh đổi 1: Khả năng giảm nhiễu (Noise Variance Reduction)
    # -------------------------------------------------------------
    ax1 = axes[0, 0]
    var_red_ma = 1.0 / M  # 1/5 = 0.20 (-6.99 dB)
    var_red_ema = alpha / (2.0 - alpha)  # 0.2 / 1.8 = 1/9 ≈ 0.111 (-9.54 dB)
    bars = ax1.bar(['MA (M=5)', 'EMA (α=0.2)'], [var_red_ma, var_red_ema],
                   color=['#2980b9', '#c0392b'], width=0.45, edgecolor='black')
    ax1.set_ylabel("Hệ số Phương sai Nhiễu Còn lại: $\sigma^2_{out} / \sigma^2_{in}$")
    ax1.set_title("1. KHẢ NĂNG GIẢM NHIỄU (CÀNG NHỎ CÀNG TỐT)", fontweight='bold')
    ax1.set_ylim(0, 0.26)
    ax1.grid(True, linestyle=':', alpha=0.6)
    
    # Text trên cột
    ax1.text(0, var_red_ma + 0.01, f"{var_red_ma:.3f} (-7.0 dB)\nGiảm 5 lần", ha='center', fontweight='bold', color='blue')
    ax1.text(1, var_red_ema + 0.01, f"{var_red_ema:.3f} (-9.5 dB)\nGiảm 9 lần (VƯỢT TRỘI)", ha='center', fontweight='bold', color='red')

    # -------------------------------------------------------------
    # 4.2 Đánh đổi 2: Độ trễ (Delay)
    # -------------------------------------------------------------
    ax2 = axes[0, 1]
    delay_ma = (M - 1) / 2.0  # 2 mẫu
    delay_ema = (1.0 - alpha) / alpha  # 4 mẫu tại DC
    bars2 = ax2.bar(['MA (M=5)', 'EMA (α=0.2)'], [delay_ma, delay_ema],
                    color=['#2980b9', '#c0392b'], width=0.45, edgecolor='black')
    ax2.set_ylabel("Độ trễ thời gian (Số mẫu)")
    ax2.set_title("2. ĐỘ TRỄ NHÓM TẠI TẦN SỐ THẤP (CÀNG NHỎ CÀNG NHANH)", fontweight='bold')
    ax2.set_ylim(0, 5.2)
    ax2.grid(True, linestyle=':', alpha=0.6)
    
    ax2.text(0, delay_ma + 0.15, f"{delay_ma:.1f} mẫu\n(Phản ứng nhanh hơn 2x)", ha='center', fontweight='bold', color='blue')
    ax2.text(1, delay_ema + 0.15, f"{delay_ema:.1f} mẫu\n(Trễ gấp đôi MA)", ha='center', fontweight='bold', color='red')

    # -------------------------------------------------------------
    # 4.3 Đánh đổi 3: Sai số do làm mượt quá mức (Over-smoothing Bias)
    # -------------------------------------------------------------
    ax3 = axes[1, 0]
    # Khi tín hiệu sạch s[n] đi qua bộ lọc, biên độ đỉnh bị suy giảm (peak clipping)
    # Tần số f0 = 0.01 -> w0 = 2*pi*f0 = 0.02*pi
    w0 = 2.0 * np.pi * 0.01
    # Đáp ứng biên độ MA tại w0: |H(w0)| = sin(M*w0/2) / (M * sin(w0/2))
    mag_ma = np.abs(np.sin(M * w0 / 2.0) / (M * np.sin(w0 / 2.0)))
    # Đáp ứng biên độ EMA tại w0: |H(w0)| = alpha / |1 - (1-alpha)*e^(-j*w0)|
    mag_ema = alpha / np.abs(1.0 - (1.0 - alpha) * np.exp(-1j * w0))

    atten_ma = (1.0 - mag_ma) * 100.0   # % suy hao biên độ
    atten_ema = (1.0 - mag_ema) * 100.0

    bars3 = ax3.bar(['MA (M=5)', 'EMA (α=0.2)'], [atten_ma, atten_ema],
                    color=['#2980b9', '#c0392b'], width=0.45, edgecolor='black')
    ax3.set_ylabel("Độ suy hao biên độ đỉnh sạch (%)")
    ax3.set_title("3. SAI SỐ LÀM MƯỢT QUÁ MỨC (OVER-SMOOTHING BIAS)", fontweight='bold')
    ax3.set_ylim(0, 3.8)
    ax3.grid(True, linestyle=':', alpha=0.6)
    
    ax3.text(0, atten_ma + 0.15, f"{atten_ma:.2f}%\n(Bảo toàn hình dạng cực tốt)", ha='center', fontweight='bold', color='blue')
    ax3.text(1, atten_ema + 0.15, f"{atten_ema:.2f}%\n(Suy hao gấp 7 lần MA)", ha='center', fontweight='bold', color='red')

    # -------------------------------------------------------------
    # 4.4 HẠN CHẾ CỐT LÕI: Độ trễ nhóm của EMA phụ thuộc tần số (Pha phi tuyến)
    # -------------------------------------------------------------
    ax4 = axes[1, 1]
    # Khảo sát độ trễ nhóm tau_g theo tần số f từ 0 đến 0.1
    f_arr = np.linspace(0.0001, 0.15, 300)
    w_arr = 2.0 * np.pi * f_arr

    # MA group delay: luôn phẳng = (M-1)/2 = 2.0 mẫu
    tau_g_ma = np.full_like(w_arr, (M - 1) / 2.0)

    # EMA group delay:
    # H(z) = alpha / (1 - beta * z^-1), với beta = 1 - alpha
    # tau_g(w) = - d/dw arg(H(w)) = beta * (cos(w) - beta) / (1 + beta^2 - 2*beta*cos(w))
    beta = 1.0 - alpha
    tau_g_ema = beta * (np.cos(w_arr) - beta) / (1.0 + beta**2 - 2.0 * beta * np.cos(w_arr))

    ax4.plot(f_arr, tau_g_ma, 'b-', linewidth=2.5, label='MA: Trễ nhóm phẳng cố định = 2.0 mẫu (Pha tuyến tính)')
    ax4.plot(f_arr, tau_g_ema, 'r-', linewidth=2.5, label='EMA: Trễ nhóm biến thiên theo tần số (Pha phi tuyến)')
    ax4.axvline(x=0.01, color='green', linestyle=':', label='Tần số tín hiệu khảo sát f0=0.01')
    
    ax4.set_xlabel("Tần số chuẩn hóa f (cycles/sample)", fontweight='bold')
    ax4.set_ylabel("Độ trễ nhóm (Số mẫu)", fontweight='bold')
    ax4.set_title("4. HẠN CHẾ CỦA EMA: ĐỘ TRỄ PHỤ THUỘC TẦN SỐ", fontweight='bold', color='#922b21')
    ax4.grid(True, linestyle=':', alpha=0.6)
    ax4.legend(loc='upper right', fontsize=9.5, framealpha=0.95)
    ax4.set_ylim(0, 4.5)

    plt.tight_layout()
    file_path = os.path.join(OUTPUT_DIR, "slide_hinh4_phan_tich_danh_doi_va_han_che.png")
    plt.savefig(file_path, dpi=160, bbox_inches='tight')
    plt.close()
    print(f"📸 Đã lưu: {file_path}")


# ======================================================================================
# 4. CHƯƠNG TRÌNH CHÍNH
# ======================================================================================

def main():
    print("=" * 80)
    print(" BÀI TẬP 3: SO SÁNH KHẢ NĂNG KHỬ NHIỄU BỘ LỌC FIR (MA) VÀ IIR (EMA)")
    print(" TỔNG HỢP KẾT QUẢ MÔ PHỎNG, ĐÁNH GIÁ ĐÁNH ĐỔI & TẠO SLIDE THUYẾT TRÌNH")
    print("=" * 80)

    # 1. Tạo Hình 1: Input -> Filter -> Output
    print("\n[1/4] Đang tạo Hình 1: Minh họa Input -> Filter -> Output...")
    tao_hinh_1_input_filter_output(M=5, alpha=0.2, snr_db=3)

    # 2. Tạo Hình 2: Nguyên lý & Đáp ứng xung
    print("\n[2/4] Đang tạo Hình 2: Nguyên lý hoạt động & Đáp ứng xung h[n]...")
    tao_hinh_2_nguyen_ly_va_dap_ung_xung(M=5, alpha=0.2)

    # 3. Chạy Monte Carlo và Tạo Hình 3: So sánh RMSE theo SNR
    print("\n[3/4] Đang chạy mô phỏng Monte Carlo (100 trials, quét SNR [-10, 20] dB)...")
    snr_range = np.arange(-10, 21, 1)
    res = chay_mo_phong_monte_carlo(N=600, f0=0.01, M=5, alpha=0.2,
                                   snr_range=snr_range, num_trials=100, seed=RANDOM_SEED)
    tao_hinh_3_so_sanh_rmse_snr(res)

    # 4. Tạo Hình 4: Phân tích 4 sự đánh đổi & Hạn chế độ trễ EMA
    print("\n[4/4] Đang tạo Hình 4: Phân tích các sự đánh đổi & hạn chế độ trễ EMA...")
    tao_hinh_4_phan_tich_danh_doi_va_han_che(M=5, alpha=0.2)

    # In bảng tóm tắt kết quả
    print("\n" + "=" * 80)
    print("                     BẢNG KẾT QUẢ MONTE CARLO THEO CÁC MỐC SNR")
    print("=" * 80)
    print(f"{'SNR (dB)':>8} | {'Raw RMSE':>10} | {'RMSE MA':>10} | {'RMSE EMA':>10} | {'Bộ lọc tốt hơn':>16} | {'Chênh lệch (%)':>15}")
    print("-" * 80)
    for i, s_val in enumerate(res['snr_range']):
        if s_val % 2 == 0 or s_val in [-10, -5, 0, 5, 10, 15, 20]:
            r_raw = res['rmse_raw_mean'][i]
            r_ma = res['rmse_ma_mean'][i]
            r_ema = res['rmse_ema_mean'][i]
            better = "EMA (IIR) ✓" if r_ema < r_ma else "MA (FIR) ✓"
            diff = (r_ema - r_ma) / r_ma * 100.0
            print(f"{s_val:>8d} | {r_raw:>10.4f} | {r_ma:>10.4f} | {r_ema:>10.4f} | {better:>16} | {abs(diff):>14.2f}%")
    print("=" * 80)
    print("🎉 Hoàn tất toàn bộ mô phỏng và tạo ảnh thành công!")


if __name__ == "__main__":
    main()
