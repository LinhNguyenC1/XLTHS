# XLTHS - Bài tập 3: So sánh MA (FIR) và EMA (IIR)

## Chạy nhanh (1 lệnh)
```bash
cd /home/runner/work/XLTHS/XLTHS
python bai_tap_3_so_sanh_FIR_IIR.py
```

## Cài đặt
```bash
python -m pip install numpy matplotlib
```

## Kết quả sinh ra
Mặc định lưu trong thư mục tương đối `outputs/bai_tap_3/`:
- `slide_hinh1_input_filter_output.png` (a) input->filter->output
- `slide_hinh2_nguyen_ly_dap_ung_xung.png` (b) MA, (c) EMA
- `slide_hinh3_so_sanh_rmse_snr.png` (e) RMSE theo SNR
- `slide_hinh4_phan_tich_danh_doi_va_han_che.png` (f) trade-off
- `ket_qua_rmse_theo_snr.csv`
- `tong_hop_ket_qua.json`

## Tham số mặc định đúng đề
- MA: `M=5`
- EMA: `alpha=0.2`
- SNR: `-10..20 dB` (bước 1 dB)
- Monte Carlo: `100 trials/SNR`
- Seed tái lập: `42`
- Metric chính: **aligned RMSE** (có căn chỉnh trễ + bỏ quá độ)

## Lưu ý học thuật quan trọng
- MA có group delay chính xác: `(M-1)/2 = 2` mẫu.
- EMA có group delay **xấp xỉ thấp tần**: `(1-alpha)/alpha = 4` mẫu khi `alpha=0.2`.
- Delay EMA phụ thuộc tần số (pha phi tuyến), nên không phải delay cố định cho mọi thành phần tần số.
- Công thức giảm phương sai nhiễu trắng (noise-only):
  - MA: `1/M`
  - EMA: `alpha/(2-alpha)`
  Các công thức này **không đồng nhất** với RMSE tổng (vì RMSE còn gồm bias/over-smoothing + hiệu ứng biên độ).

## Self-check tối thiểu
```bash
python bai_tap_3_so_sanh_FIR_IIR.py --self-check
```
Kiểm tra:
1. Tổng hệ số MA bằng 1.
2. Tổng đáp ứng xung EMA xấp xỉ 1.
3. SNR sinh ra gần SNR mục tiêu.
4. Shape input/output nhất quán.

## Tùy chọn thường dùng
```bash
# đổi thư mục output
python bai_tap_3_so_sanh_FIR_IIR.py --output-dir outputs/custom

# chỉ chạy số liệu (không vẽ hình)
python bai_tap_3_so_sanh_FIR_IIR.py --skip-figures

# đổi số trials/seed
python bai_tap_3_so_sanh_FIR_IIR.py --trials 200 --seed 123
```

## Luồng đề xuất cho slide
(a) input->filter->output → (b) MA → (c) EMA → (d) Monte Carlo + alignment → (e) RMSE theo SNR → (f) trade-off → (g) ứng dụng/kết luận.
