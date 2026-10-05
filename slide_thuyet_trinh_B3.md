# Slide thuyết trình Bài 3 (khớp code hiện tại)

## Thứ tự đề xuất
1. **(a) Input -> Filter -> Output**
   - File: `slide_hinh1_input_filter_output.png`
2. **(b) MA/FIR: nguyên lý**
   - File: `slide_hinh2_nguyen_ly_dap_ung_xung.png` (panel MA)
3. **(c) EMA/IIR: nguyên lý**
   - File: `slide_hinh2_nguyen_ly_dap_ung_xung.png` (panel EMA)
4. **(d) Phương pháp Monte Carlo + Alignment**
   - `M=5`, `alpha=0.2`, SNR `-10..20 dB`, `100 trials/SNR`
   - RMSE chính thức: **aligned RMSE** (cắt trễ + transient + cùng vùng chỉ số)
5. **(e) RMSE theo SNR**
   - File: `slide_hinh3_so_sanh_rmse_snr.png`
   - Trình bày theo dữ liệu thực tế từ CSV/JSON, không khẳng định thắng/thua trước khi chạy.
6. **(f) So sánh trade-off**
   - File: `slide_hinh4_phan_tich_danh_doi_va_han_che.png`
   - Nội dung: giảm nhiễu, độ trễ, over-smoothing, độ phức tạp.
7. **(g) Ứng dụng/kết luận**
   - Kết luận có điều kiện theo dải SNR quan sát được trong lần chạy.

## Ghi chú học thuật bắt buộc khi thuyết trình
- MA: trễ nhóm cố định `2` mẫu (với `M=5`).
- EMA: trễ nhóm thấp tần xấp xỉ `4` mẫu (với `alpha=0.2`), **phụ thuộc tần số**.
- Công thức `1/M` và `alpha/(2-alpha)` chỉ mô tả giảm phương sai nhiễu trắng, không đồng nhất với RMSE tổng.
- EMA khởi tạo mặc định `y[-1]=0` (có thể đổi bằng CLI `--ema-init first_sample`).

## Dữ liệu dùng cho slide
- `ket_qua_rmse_theo_snr.csv`: số liệu bảng theo từng SNR.
- `tong_hop_ket_qua.json`: config + kết quả + phân tích giao cắt MA/EMA.

## Lệnh chạy tái lập
```bash
python bai_tap_3_so_sanh_FIR_IIR.py --seed 42 --trials 100
```
