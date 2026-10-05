# BÁO CÁO THUYẾT TRÌNH & DÀN Ý SLIDE: BÀI TẬP 3
## SO SÁNH KHẢ NĂNG KHỬ NHIỄU CỦA BỘ LỌC FIR (MA) VÀ IIR (EMA)

---

## 📑 BỐ CỤC SLIDE TRÌNH BÀY (DỰ KIẾN 8 SLIDE)

| Slide | Tiêu đề nội dung | Trọng tâm trình bày |
| :---: | :--- | :--- |
| **Slide 1** | **Tổng quan & Tác động của 2 bộ lọc lên tín hiệu** | Tập trung vào **KẾT QUẢ TRỰC QUAN**: $Input \to [Bộ\ lọc] \to Output$ |
| **Slide 2** | **Nguyên lý hoạt động Bộ lọc 1: Moving Average (MA - FIR)** | Cửa sổ trượt, trung bình cộng, $h[n]$ hữu hạn, không phản hồi |
| **Slide 3** | **Nguyên lý hoạt động Bộ lọc 2: Exponential MA (EMA - IIR)** | Phương trình sai phân đệ quy, $h[n]$ vô hạn suy giảm mũ, phản hồi $y[n-1]$ |
| **Slide 4** | **Phương pháp thực nghiệm: Mô phỏng Monte Carlo & Căn chỉnh trễ** | 100 trials, quét SNR [-10, 20] dB, căn chỉnh trễ 2 & 4 mẫu, bỏ quá độ |
| **Slide 5** | **Kết quả thực nghiệm: Đường cong RMSE theo SNR** | 2 đường trên cùng hệ trục; Phân tích 3 vùng SNR (Thấp, Trung bình, Cao) |
| **Slide 6** | **So sánh toàn diện 6 tiêu chí & Phân tích các sự đánh đổi (Trade-offs)** | Giảm nhiễu, Độ trễ, Sai số làm mượt quá mức, Độ phức tạp tính toán |
| **Slide 7** | **Hạn chế cốt lõi của EMA: Độ trễ nhóm phụ thuộc tần số** | Pha phi tuyến của IIR $\to$ méo phân tán (dispersion) |
| **Slide 8** | **Tổng kết & Ứng dụng thực tế** | Khi nào nên dùng MA? Khi nào nên dùng EMA? |

---

## 🖥️ CHI TIẾT TỪNG SLIDE THUYẾT TRÌNH

### SLIDE 1: TÁC ĐỘNG CỦA 2 BỘ LỌC LÊN TÍN HIỆU (KẾT QUẢ TRỰC QUAN)

#### 1. Sơ đồ chuỗi xử lý tín hiệu
$$\text{Tín hiệu sạch } s[n] \xrightarrow{+ \text{Nhiễu AWGN } v[n]} \text{Input } x[n] \xrightarrow{\text{Bộ lọc MA / EMA}} \text{Output } y[n]$$

#### 2. Kết quả quan sát trực quan từ hình ảnh (`slide_hinh1_input_filter_output.png`)
- **Đầu vào (Input $x[n]$):** Tín hiệu quan sát bị nhiễu dao động mạnh, che lấp dạng sóng chậm gốc $s[n]$.
- **Sau khi qua 2 bộ lọc (Output $y[n]$):**
  - Cả 2 bộ lọc đều **triệt tiêu dao động tần số cao** của nhiễu Gaussian trắng, khôi phục lại dáng điệu mượt mà của sóng sin gốc.
  - **EMA ($\alpha=0.2$):** Làm mượt "sâu" hơn, đường cong ra rất êm nhưng bị dịch trễ về bên phải rõ rệt hơn.
  - **MA ($M=5$):** Bám nhanh hơn (ít trễ hơn), nhưng đường cong vẫn còn gợn sóng nhẹ (nhiễu sót lại nhiều hơn một chút).
- **Phần quá độ ban đầu (Transient):** Trong những mẫu đầu tiên, tín hiệu ra phải mất một thời gian nạp bộ nhớ mới đạt trạng thái bám ổn định.

---

### SLIDE 2: NGUYÊN LÝ HOẠT ĐỘNG BỘ LỌC 1 — MOVING AVERAGE (MA / FIR)

#### 1. Công thức toán học
$$y_1[n] = \frac{1}{M} \sum_{k=0}^{M-1} x[n-k] = \frac{x[n] + x[n-1] + \dots + x[n-M+1]}{M}$$

#### 2. Bản chất kỹ thuật FIR (Finite Impulse Response)
- **Cơ chế:** Dùng một "cửa sổ trượt" kích thước $M$ mẫu thời gian. Mỗi bước thời gian lấy trung bình cộng đều nhau (trọng số $1/M = 0.2$ với $M=5$).
- **Không có phản hồi (Non-recursive / Feedforward only):** Đầu ra $y[n]$ chỉ tính từ các mẫu đầu vào quá khứ và hiện tại $x[n-k]$, **không dùng** $y[n-1]$.
- **Đáp ứng xung hữu hạn $h[n]$:**
  $$h[n] = \begin{cases} \frac{1}{M} & \text{với } 0 \le n \le M-1 \\ 0 & \text{với } n \ge M \end{cases}$$
  Chiều dài xung đúng $M=5$ mẫu rồi tắt hẳn.
- **Tính ổn định & Pha:**
  - Luôn luôn ổn định tuyệt đối (không có cực ngoài gốc tọa độ, hệ thống All-Zero).
  - **Pha tuyến tính tuyệt đối (Linear Phase):** Vì hệ số $h[n]$ đối xứng, mọi tần số đều trễ một khoảng thời gian cố định:
    $$\tau = \frac{M - 1}{2} = \frac{5 - 1}{2} = 2.0 \text{ mẫu}$$
    $\implies$ **Không gây méo dạng sóng do trễ pha!**

---

### SLIDE 3: NGUYÊN LÝ HOẠT ĐỘNG BỘ LỌC 2 — EXPONENTIAL MOVING AVERAGE (EMA / IIR)

#### 1. Công thức sai phân bậc 1
$$y_2[n] = \alpha \cdot x[n] + (1 - \alpha) \cdot y_2[n-1]$$
- Với $\alpha = 0.2$: $y_2[n] = 0.2 \cdot x[n] + 0.8 \cdot y_2[n-1]$.
- Giá trị mới bằng: **20% mẫu hiện tại** cộng với **80% "ký ức" đầu ra trước đó**.

#### 2. Bản chất kỹ thuật IIR (Infinite Impulse Response)
- **Cơ chế đệ quy (Recursive / Feedback):** Tín hiệu đầu ra trước đó $y[n-1]$ được đưa ngược trở lại làm đầu vào của bước tiếp theo.
- **Khai triển chuỗi thời gian:**
  $$y_2[n] = \alpha x[n] + \alpha(1-\alpha)x[n-1] + \alpha(1-\alpha)^2 x[n-2] + \dots = \sum_{k=0}^{\infty} \alpha(1-\alpha)^k x[n-k]$$
- **Đáp ứng xung vô hạn $h[n]$:**
  $$h[n] = \alpha (1 - \alpha)^n u[n]$$
  Các mẫu quá khứ không bị cắt đột ngột như MA, mà trọng số suy giảm theo hàm mũ, kéo dài vô hạn về mặt lý thuyết.
- **Độ trễ và Pha:**
  - Hằng số thời gian: $\tau = -\frac{1}{\ln(1-\alpha)} \approx 4.48$ mẫu.
  - Độ trễ nhóm ở tần số rất thấp:
    $$\tau_g \approx \frac{1 - \alpha}{\alpha} = \frac{1 - 0.2}{0.2} = 4.0 \text{ mẫu}$$
  - **Pha phi tuyến (Non-linear Phase):** Độ trễ nhóm thay đổi theo tần số $\implies$ gây méo phân tán (sẽ phân tích kỹ ở Slide 7).

---

### SLIDE 4: PHƯƠNG PHÁP THỰC NGHIỆM — MÔ PHỎNG MONTE CARLO & CĂN CHỈNH TRỄ

#### 1. Thiết lập mô phỏng
- **Tín hiệu sạch:** $s[n] = \sin(2\pi f_0 n)$ với $f_0 = 0.01$ (chu kỳ $T = 100$ mẫu), tổng độ dài $N = 600$ mẫu.
- **Nhiễu AWGN:** Nhiễu Gauss cộng $v[n] \sim \mathcal{N}(0, \sigma_v^2)$ với công suất tính theo SNR:
  $$\sigma_v = \sqrt{\frac{P_s}{10^{\text{SNR}/10}}}$$
- **Quét SNR:** Từ $-10\text{ dB}$ đến $+20\text{ dB}$ (31 mức, bước 1 dB).
- **Số lần lặp Monte Carlo:** 100 trials độc lập cho mỗi mức SNR.
- **Seed ngẫu nhiên:** Cố định `np.random.seed(42)` để kết quả độc lập, khách quan và hoàn toàn tái lập được.

#### 2. Hai bước tiền xử lý bắt buộc trước khi tính RMSE
1. **Căn chỉnh độ trễ (Delay Alignment):**
   - MA: Dịch sang trái $\Delta = 2$ mẫu.
   - EMA: Dịch sang trái $\Delta = 4$ mẫu.
   - *Nếu không bù trễ, RMSE sẽ đo sai lệch pha thay vì đo hiệu quả khử nhiễu!*
2. **Bỏ qua phần quá độ ban đầu (Transient Phase):**
   - Cắt bỏ 50 mẫu đầu tiên ($n_{\text{transient}} = 50$).
   - Đảm bảo hệ thống đạt trạng thái xác lập (steady-state), triệt tiêu sai số nạp ban đầu $y[-1]=0$.

---

### SLIDE 5: KẾT QUẢ MÔ PHỎNG — ĐƯỜNG CONG RMSE THEO SNR

*(Minh họa bằng đồ thị `slide_hinh3_so_sanh_rmse_snr.png`)*

#### Bảng đối chiếu số liệu thực nghiệm (Trích xuất từ 100 trials Monte Carlo)

| SNR (dB) | Chưa lọc (Raw RMSE) | RMSE MA ($M=5$) | RMSE EMA ($\alpha=0.2$) | Bộ lọc tốt hơn | Chênh lệch (%) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **-10 dB** | 2.2386 | 1.0010 | **0.7482** | **EMA thắng áp đảo** | EMA thấp hơn **25.26%** |
| **-5 dB** | 1.2599 | 0.5586 | **0.4149** | **EMA thắng áp đảo** | EMA thấp hơn **25.72%** |
| **0 dB** | 0.7079 | 0.3184 | **0.2379** | **EMA thắng áp đảo** | EMA thấp hơn **25.30%** |
| **+5 dB** | 0.3971 | 0.1786 | **0.1361** | **EMA thắng** | EMA thấp hơn **23.78%** |
| **+10 dB** | 0.2239 | 0.1001 | **0.0796** | **EMA thắng** | EMA thấp hơn **20.49%** |
| **+15 dB** | 0.1256 | 0.0562 | **0.0497** | **EMA thắng** | EMA thấp hơn **11.72%** |
| **+18 dB** | 0.0891 | **0.0398** | 0.0401 | **MA vượt lên** | MA thấp hơn **0.71%** |
| **+20 dB** | 0.0706 | **0.0317** | 0.0357 | **MA thắng rõ rệt** | MA thấp hơn **12.70%** |

#### Nhận xét chi tiết tại 3 vùng SNR:
1. **Vùng SNR thấp (<-10 dB đến 0 dB) — "Vùng nhiễu thống trị":**
   - Nhiễu rất mạnh, tín hiệu bị chôn vùi.
   - **EMA thắng áp đảo** với RMSE thấp hơn MA xấp xỉ **25.3%** tại mọi điểm.
   - *Lý do:* Hệ số triệt tiêu phương sai nhiễu của EMA ($1/9 \approx 0.111$) tốt hơn nhiều so với MA ($1/5 = 0.200$). Tỷ số lý thuyết: $\sqrt{5/9} \approx 0.745 \implies$ sai số giảm đúng 25.5%!
2. **Vùng SNR trung bình (0 dB đến 15 dB) — "Vùng cân bằng":**
   - Tín hiệu và nhiễu cạnh tranh nhau.
   - EMA vẫn duy trì ưu thế giảm nhiễu tốt hơn MA, nhưng khoảng cách dần thu hẹp từ 25% xuống 11%.
3. **Vùng SNR cao (16 dB đến 20 dB) — "Vùng tín hiệu sạch":**
   - Nhiễu ngẫu nhiên đã rất nhỏ. Sai số lúc này không còn do nhiễu, mà do **sai số làm mượt quá mức (Over-smoothing bias)** của bộ lọc gây ra trên tín hiệu sạch.
   - **MA vượt lên thắng EMA:** Tại 20 dB, MA tốt hơn EMA **12.7%** vì MA bảo toàn biên độ sóng tốt hơn nhiều (MA chỉ làm suy hao đỉnh 0.41%, trong khi EMA làm suy hao tới 2.95%).

---

### SLIDE 6: SO SÁNH TOÀN DIỆN 6 TIÊU CHÍ & SỰ ĐÁNH ĐỔI (TRADE-OFFS)

*(Minh họa bằng đồ thị `slide_hinh4_phan_tich_danh_doi_va_han_che.png`)*

| Tiêu chí so sánh | Bộ lọc MA ($M=5$) | Bộ lọc EMA ($\alpha=0.2$) | Ý nghĩa vật lý & Sự đánh đổi |
| :--- | :--- | :--- | :--- |
| **1. Khả năng giảm nhiễu** | $\sigma^2_{\text{out}}/\sigma^2_{\text{in}} = \frac{1}{M} = \mathbf{0.20}$ (-7.0 dB) | $\sigma^2_{\text{out}}/\sigma^2_{\text{in}} = \frac{\alpha}{2-\alpha} = \mathbf{0.111}$ (-9.5 dB) | **EMA vượt trội:** Triệt tiêu phương sai nhiễu gấp 1.8 lần so với MA |
| **2. Độ trễ thời gian (Delay)** | Cố định **2.0 mẫu** | Xấp xỉ **4.0 mẫu** (tại tần số thấp) | **MA nhanh hơn gấp đôi:** EMA phản ứng chậm chạp hơn, trễ nhiều hơn |
| **3. Sai số làm mượt quá mức** | Suy hao đỉnh: **0.41%** (Rất nhỏ) | Suy hao đỉnh: **2.95%** (Gấp 7 lần MA) | **MA vượt trội ở tín hiệu sạch:** Ít làm tù đỉnh, ít méo biên độ |
| **4. Độ phức tạp tính toán** | $M$ phép cộng, 1 phép chia mỗi mẫu (hoặc đệm vòng) | **1 nhân, 1 trừ, 1 cộng** (hoặc 2 nhân, 1 cộng) | **EMA cực nhẹ:** Số phép tính không phụ thuộc vào độ dài cửa sổ |
| **5. Chi phí bộ nhớ (RAM)** | Cần mảng nhớ lưu $M=5$ mẫu cũ $x[n-k]$ | **Chỉ cần 1 biến duy nhất** lưu $y[n-1]$ | **EMA tối ưu phần cứng:** Thích hợp nhất cho vi điều khiển (MCU, IoT) |
| **6. Ứng dụng tiêu biểu** | Xử lý ảnh, xử lý tín hiệu sau thu thập (Offline), DSP pha tuyến tính | Cảm biến thời gian thực, IoT, điều khiển nhúng, tài chính (Trading) | Phụ thuộc vào ưu tiên: Cần pha chuẩn $\to$ MA; Cần nhẹ & mượt $\to$ EMA |

---

### SLIDE 7: HẠN CHẾ CỐT LÕI CỦA EMA — ĐỘ TRỄ NHÓM PHỤ THUỘC TẦN SỐ

#### 1. Nguyên nhân lý thuyết
- Hàm truyền của EMA:
  $$H_{\text{EMA}}(e^{j\omega}) = \frac{\alpha}{1 - (1-\alpha)e^{-j\omega}}$$
- Độ trễ nhóm (Group Delay) của EMA:
  $$\tau_g(\omega) = -\frac{d}{d\omega} \arg[H(\omega)] = \frac{(1-\alpha)(\cos\omega - (1-\alpha))}{1 + (1-\alpha)^2 - 2(1-\alpha)\cos\omega}$$
- Tại $\omega \to 0$ (DC): $\tau_g(0) = \frac{1-\alpha}{\alpha} = \frac{0.8}{0.2} = 4.0$ mẫu.
- Nhưng khi tần số tăng lên: **Độ trễ nhóm $\tau_g(\omega)$ GIẢM DẦN** (xem Subplot 4 trên `slide_hinh4_phan_tich_danh_doi_va_han_che.png`).

#### 2. Hậu quả thực tế (Pha phi tuyến)
- Các thành phần hài có tần số khác nhau trong tín hiệu sẽ bị **trễ những khoảng thời gian khác nhau**.
- Gây ra **hiện tượng tán sắc / méo dạng pha (Phase Dispersion)**: Nếu tín hiệu đầu vào là âm thanh hoặc xung nhịp tim phức tạp gồm nhiều tần số, EMA sẽ làm tín hiệu bị biến dạng dáng sóng (shape distortion).
- Ngược lại, **MA (FIR)** có đáp ứng pha tuyến tính tuyệt đối, độ trễ nhóm phẳng cố định 2.0 mẫu ở mọi tần số, hoàn toàn không gây méo tán sắc!

---

### SLIDE 8: TỔNG KẾT & KẾT LUẬN

1. **Hiểu bản chất hai bộ lọc:**
   - Cả MA và EMA đều là bộ lọc thông thấp (Low-pass smoothing filters) hiệu quả trên miền thời gian.
   - MA là đại diện tiêu biểu của **FIR** (cửa sổ hữu hạn, feedforward).
   - EMA là đại diện tiêu biểu của **IIR** (đáp ứng vô hạn, feedback đệ quy).
2. **Quy luật đánh đổi (Trade-off Law):**
   - Không có bộ lọc nào hoàn hảo toàn diện: Muốn giảm nhiễu mạnh (EMA) thì phải chấp nhận độ trễ lớn hơn và sai số làm mượt ở vùng sạch.
   - Muốn pha tuyến tính không méo hình (MA) thì phải tốn bộ nhớ lưu đệm mẫu và chấp nhận nhiễu sót lại nhiều hơn.
3. **Ý nghĩa của việc căn chỉnh trễ:**
   - Khi đánh giá thuật toán khử nhiễu, bắt buộc phải bù trễ pha (MA trễ 2 mẫu, EMA trễ 4 mẫu) và loại bỏ quá độ; nếu không kết quả RMSE sẽ bị bóp méo hoàn toàn bởi sai số pha thay vì sai số biên độ nhiễu.
4. **Khuyến nghị lựa chọn:**
   - **Chọn EMA khi:** Thiết bị tài nguyên hạn chế (MCU, IoT sensor, chip nhúng), cần thuật toán chạy online/real-time cực nhanh, môi trường có mức nhiễu cao ($SNR < 15\text{ dB}$).
   - **Chọn MA khi:** Tín hiệu nhạy cảm với pha (âm thanh, y sinh ECG/EEG), cần pha tuyến tính chuẩn, hoặc tín hiệu có mức SNR cao cần bảo toàn biên độ đỉnh.

---
*Tài liệu được chuẩn bị phục vụ báo cáo bài tập lớn Xử lý Tiếng nói & Tín hiệu số.*
