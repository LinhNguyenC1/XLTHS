# Bài Tập 3: So sánh khả năng khử nhiễu của bộ lọc FIR và IIR cơ bản
## Tài liệu giải thích chi tiết các thuật ngữ và ý nghĩa

---

## 1. TỔNG QUAN BÀI TẬP

### Mục tiêu
So sánh hiệu năng khử nhiễu của **hai bộ lọc làm mượt** (smoothing filter) cơ bản nhất trong xử lý tín hiệu:

| Bộ lọc | Loại | Công thức | Bản chất |
|--------|------|-----------|----------|
| **Moving Average (MA)** | FIR | `y1[n] = (1/M) Σ x[n-k]` | Lấy trung bình M mẫu gần nhất |
| **Exponential MA (EMA)** | IIR | `y2[n] = α·x[n] + (1-α)·y2[n-1]` | Trọng số giảm dần theo hàm mũ |

### Phương pháp
Sử dụng **mô phỏng Monte Carlo** (100 trials) để đánh giá thống kê chính xác hiệu năng ở nhiều mức **SNR** khác nhau (-10 dB đến 20 dB).

---

## 2. CÁC THUẬT NGỮ QUAN TRỌNG

### 2.1 FIR (Finite Impulse Response - Đáp ứng xung hữu hạn)

> **Định nghĩa**: Bộ lọc mà đáp ứng xung h[n] có chiều dài **hữu hạn** (kết thúc sau M mẫu).

**Đặc điểm chính:**
- Đầu ra `y[n]` chỉ phụ thuộc vào các mẫu **đầu vào** `x[n], x[n-1], ..., x[n-M+1]`
- **KHÔNG** có phản hồi (feedback) từ đầu ra
- Dạng tổng quát: `y[n] = b₀·x[n] + b₁·x[n-1] + ... + b_{M-1}·x[n-M+1]`
- Luôn **ổn định** (stable) vì không có cực (pole) nào
- Có thể thiết kế với **pha tuyến tính** → không méo dạng tín hiệu

**Ví dụ cụ thể - Bộ lọc MA (M=5):**
```
y[n] = (1/5) * (x[n] + x[n-1] + x[n-2] + x[n-3] + x[n-4])

Đáp ứng xung: h = [0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, ...]
                    ← 5 giá trị khác 0 → HỮU HẠN
```

### 2.2 IIR (Infinite Impulse Response - Đáp ứng xung vô hạn)

> **Định nghĩa**: Bộ lọc mà đáp ứng xung h[n] kéo dài **vô hạn** (không bao giờ chấm dứt hoàn toàn).

**Đặc điểm chính:**
- Đầu ra `y[n]` phụ thuộc vào cả mẫu **đầu vào** `x[n]` VÀ mẫu **đầu ra trước đó** `y[n-1]`
- **CÓ** phản hồi (feedback/recursion) → đây là đặc trưng cốt lõi
- Dạng tổng quát: `y[n] = b₀·x[n] + ... + b_M·x[n-M] - a₁·y[n-1] - ... - a_N·y[n-N]`
- Có thể **không ổn định** nếu thiết kế sai
- Pha **phi tuyến** → có thể gây méo dạng nhẹ

**Ví dụ cụ thể - Bộ lọc EMA (α=0.2):**
```
y[n] = 0.2 * x[n] + 0.8 * y[n-1]
       ↑ mẫu mới      ↑ "ký ức" cũ

Đáp ứng xung: h[n] = 0.2 * (0.8)^n
h = [0.200, 0.160, 0.128, 0.102, 0.082, 0.066, 0.052, ...]
    → Giảm dần theo hàm mũ nhưng KHÔNG BAO GIỜ = 0 → VÔ HẠN
```

### 2.3 AWGN (Additive White Gaussian Noise)

> **Định nghĩa**: Mô hình nhiễu tiêu chuẩn trong xử lý tín hiệu, với 3 tính chất:

| Tính chất | Ý nghĩa | Giải thích trực quan |
|-----------|---------|---------------------|
| **Additive** (Cộng) | `x[n] = s[n] + w[n]` | Nhiễu "cộng thêm" vào tín hiệu, không nhân hay biến đổi |
| **White** (Trắng) | Phổ phẳng ở mọi tần số | Như ánh sáng trắng chứa mọi màu sắc → nhiễu chứa mọi tần số |
| **Gaussian** (Gauss) | Phân bố chuẩn N(0, σ²) | Đa số giá trị nhiễu gần 0, ít khi có giá trị cực lớn |

**Tại sao dùng AWGN?**
- Mô hình hóa nhiễu nhiệt (thermal noise) trong mạch điện tử
- Là chuẩn benchmark quốc tế: mọi thuật toán phải chứng minh hoạt động tốt với AWGN trước
- Đơn giản về mặt toán học nhưng phản ánh đúng thực tế

### 2.4 SNR (Signal-to-Noise Ratio)

> **Định nghĩa**: Tỷ số giữa công suất tín hiệu và công suất nhiễu.

**Công thức:**
```
SNR_linear = P_signal / P_noise

SNR_dB = 10 * log₁₀(P_signal / P_noise)
```

**Bảng quy đổi và ý nghĩa vật lý:**

| SNR (dB) | Tỷ lệ tuyến tính | Ý nghĩa thực tế |
|----------|-------------------|------------------|
| +20 dB | Signal gấp 100× Noise | Rất sạch, gần như không nhiễu |
| +10 dB | Signal gấp 10× Noise | Có nhiễu nhẹ, vẫn nhận diện được |
| 0 dB | Signal = Noise | Ngang nhau, bắt đầu khó phân biệt |
| -10 dB | Noise gấp 10× Signal | Tín hiệu "chìm" trong nhiễu |
| -20 dB | Noise gấp 100× Signal | Chỉ thấy nhiễu, tín hiệu gần như biến mất |

### 2.5 Mô phỏng Monte Carlo

> **Định nghĩa**: Phương pháp thống kê dùng lấy mẫu ngẫu nhiên lặp lại nhiều lần để ước lượng kết quả.

**Tại sao cần Monte Carlo trong bài tập này?**

```
Lần 1: Thêm nhiễu ngẫu nhiên #1 → RMSE = 0.312
Lần 2: Thêm nhiễu ngẫu nhiên #2 → RMSE = 0.287
Lần 3: Thêm nhiễu ngẫu nhiên #3 → RMSE = 0.324
...
Lần 100: Thêm nhiễu ngẫu nhiên #100 → RMSE = 0.301

RMSE trung bình = (0.312 + 0.287 + ... + 0.301) / 100 = 0.306 ± 0.015
```

→ Nếu chỉ chạy **1 lần**, kết quả có thể là 0.312 (may mắn) hoặc 0.350 (xui xẻo)
→ Chạy **100 lần** rồi lấy trung bình → kết quả **đáng tin cậy**, phản ánh đúng bản chất thuật toán

### 2.6 RMSE (Root Mean Square Error)

> **Định nghĩa**: Sai số toàn phương trung bình, đo khoảng cách giữa tín hiệu lọc và tín hiệu gốc.

**Công thức:**
```
RMSE = √[ (1/N) * Σ(y[n] - s[n])² ]
```

- RMSE = 0: Khôi phục hoàn hảo (không bao giờ đạt được trong thực tế)
- RMSE càng nhỏ → bộ lọc càng tốt

### 2.7 Căn chỉnh pha (Phase Alignment)

> **Vấn đề**: Bộ lọc gây ra **trễ nhóm** (group delay), làm tín hiệu đầu ra bị dịch về phía phải trên trục thời gian.

```
Tín hiệu gốc s[n]:        ╱╲    ╱╲    ╱╲
                          ╱  ╲  ╱  ╲  ╱  ╲
Sau lọc MA (trễ 2 mẫu):      ╱╲    ╱╲    ╱╲
                             ╱  ╲  ╱  ╲  ╱  ╲
                          ← 2 mẫu →
```

**Nếu KHÔNG căn chỉnh:** So sánh `y[n]` với `s[n]` → sai số bao gồm cả lỗi do lệch pha → RMSE bị phồng lên giả tạo

**Nếu CÓ căn chỉnh:** Dịch `y` sang trái (hoặc `s` sang phải) → RMSE chỉ phản ánh lỗi do nhiễu → đúng bản chất

**Độ trễ:**
- MA(M=5): delay = (M-1)/2 = **2 mẫu** (chính xác, vì pha tuyến tính)
- EMA(α=0.2): delay ≈ (1-α)/α = **4 mẫu** (xấp xỉ, vì pha phi tuyến)

---

## 3. Ý NGHĨA CỦA BÀI TẬP

### 3.1 Ý nghĩa học thuật
- Hiểu sự khác biệt cơ bản giữa bộ lọc FIR và IIR
- Nắm được phương pháp mô phỏng Monte Carlo
- Biết cách đánh giá hiệu năng bộ lọc bằng RMSE theo SNR
- Hiểu tầm quan trọng của căn chỉnh pha trong so sánh

### 3.2 Ý nghĩa thực tiễn
- **MA** được dùng phổ biến trong: phân tích tài chính (đường trung bình 50 ngày, 200 ngày), xử lý ảnh (blur filter), đo lường vật lý
- **EMA** được dùng trong: điều khiển tự động (PID controller), streaming data, IoT sensor, hệ thống thời gian thực
- Việc chọn bộ lọc phụ thuộc vào **yêu cầu cụ thể**: cần pha tuyến tính → MA, cần tiết kiệm bộ nhớ → EMA

### 3.3 Kết quả mong đợi
1. Khi SNR cao: Cả hai bộ lọc đều có RMSE thấp (lọc nhiễu dễ)
2. Khi SNR thấp: RMSE tăng lên (nhiễu quá mạnh, bộ lọc không thể khử hết)
3. Với α=0.2 và M=5 (cùng "ký ức"): Hiệu năng **tương đương** nhau
4. MA có ưu thế nhẹ về pha tuyến tính, EMA có ưu thế về tính toán

---

## 4. CÁCH CHẠY CHƯƠNG TRÌNH

```bash
cd /Users/nguyenthithuylinh/Downloads/Xử_lý_tiếng_nói/BTVN_B3_FIR_vs_IIR
python bai_tap_3_so_sanh_FIR_IIR.py
```

### Đầu ra:
1. **3 cửa sổ đồ thị**:
   - Đáp ứng xung của MA và EMA
   - So sánh RMSE theo SNR (4 subplot)
   - Quá trình lọc ở 4 mức SNR: -10, 0, 10, 20 dB
2. **Bảng kết quả** in ra terminal
3. **Nhận xét tự động** phân tích chi tiết

### Hình ảnh được lưu:
- `dap_ung_xung.png` - Đáp ứng xung h[n]
- `ket_qua_so_sanh_MA_EMA.png` - RMSE theo SNR
- `chi_tiet_loc_nhieu.png` - Ví dụ trực quan ở từng mức SNR
