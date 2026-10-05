# Bài tập 3: Giải thích chi tiết (FIR/MA vs IIR/EMA)

## 1) Thiết lập mô phỏng đúng yêu cầu
- Tín hiệu sạch chậm thay đổi: `s[n] = sin(2πf0n)`, mặc định `f0=0.01`.
- Nhiễu AWGN theo SNR đầu vào `-10..20 dB`:
  - `x[n] = s[n] + v[n]`
  - `P_noise = P_signal / 10^(SNR/10)`.
- Bộ lọc MA (FIR):
  - `y_MA[n] = (1/M) * Σ_{k=0..M-1} x[n-k]`, mặc định `M=5`.
- Bộ lọc EMA (IIR):
  - `y_EMA[n] = αx[n] + (1-α)y_EMA[n-1]`, mặc định `α=0.2`.
- Monte Carlo: `100 trials` cho mỗi mức SNR.

## 2) Căn chỉnh RMSE (metric chính)
Script dùng **aligned RMSE** là metric chính:
1. MA: trễ nhóm chính xác `Δ_MA=(M-1)/2=2 mẫu`.
2. EMA: trễ nhóm **xấp xỉ thấp tần** `Δ_EMA≈(1-α)/α=4 mẫu` (không cố định với mọi tần số).
3. Cắt bỏ quá độ khởi tạo (`n_transient=50`, với mặc định EMA khởi tạo `y[-1]=0`).
4. So sánh trên **cùng vùng chỉ số sạch** cho cả MA/EMA để công bằng.

Ngoài ra script cũng xuất thêm RMSE chưa căn chỉnh (`unaligned`) để minh họa ảnh hưởng của delay.

## 3) Công thức lý thuyết cần phân biệt
Với nhiễu trắng độc lập:
- MA: `σ²_out/σ²_in = 1/M`.
- EMA ổn định: `σ²_out/σ²_in = α/(2-α)`.

Đây là hệ số giảm phương sai **noise-only**. RMSE tổng trong mô phỏng còn phụ thuộc:
- bias do làm mượt quá mức,
- đáp ứng biên độ tại tần số tín hiệu,
- điều kiện biên/khởi tạo.

## 4) Over-smoothing/Bias đúng ngữ cảnh
Over-smoothing được tính từ đáp ứng biên độ của chính bộ lọc tại `f0`:
- MA: `|H_MA(e^{jω0})|`.
- EMA: `|H_EMA(e^{jω0})| = α / |1-(1-α)e^{-jω0}|`.

Vì vậy số liệu bias phụ thuộc `f0`, `M`, `α`; không dùng hằng số cố định ngoài ngữ cảnh.

## 5) Độ trễ và pha
- MA (FIR đối xứng) có pha tuyến tính, trễ cố định.
- EMA (IIR bậc 1) có pha phi tuyến, trễ nhóm phụ thuộc tần số.
  - Công thức trễ nhóm EMA:
  `τg(ω)=((1-α)(cosω-(1-α)))/(1+(1-α)^2-2(1-α)cosω)`.

## 6) Cách chạy
```bash
cd /home/runner/work/XLTHS/XLTHS
python bai_tap_3_so_sanh_FIR_IIR.py
```
Output mặc định trong `outputs/bai_tap_3/` gồm PNG + CSV + JSON.

Self-check:
```bash
python bai_tap_3_so_sanh_FIR_IIR.py --self-check
```
