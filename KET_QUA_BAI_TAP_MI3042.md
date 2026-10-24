# Báo cáo Chi tiết Kết quả Giải Bài tập Phương Pháp Số MI3042
**Ứng dụng:** PPS_APP (Python / Tkinter / NumPy / SciPy / Matplotlib)  
**Tất cả kết quả dưới đây đều được chạy thực tế trên ứng dụng PPS_APP và chụp ảnh màn hình giao diện thật.**

---

## Mục lục các bài đã giải
1. [Tuần 1 - Câu 3c: Mốc nội suy Chebyshev tối ưu](#1-tuần-1---câu-3c-mốc-nội-suy-chebyshev-tối-ưu)
2. [Tuần 2 - Câu 4d: Sơ đồ Horner cho đa thức bậc ba](#2-tuần-2---câu-4d-sơ-đồ-horner-cho-đa-thức-bậc-ba)
3. [Tuần 2 - Câu 7: Đa thức nội suy Lagrange](#3-tuần-2---câu-7-đa-thức-nội-suy-lagrange)
4. [Tuần 3 - Câu 8: Đa thức nội suy Newton tiến mốc bất kỳ](#4-tuần-3---câu-8-đa-thức-nội-suy-newton-tiến-mốc-bất-kỳ)
5. [Tuần 3 - Câu 9: Đa thức nội suy Newton lùi mốc bất kỳ](#5-tuần-3---câu-9-đa-thức-nội-suy-newton-lùi-mốc-bất-kỳ)
6. [Tuần 4 - Câu 12: Đa thức nội suy Gauss I](#6-tuần-4---câu-12-đa-thức-nội-suy-gauss-i)
7. [Tuần 4 - Câu 13: Đa thức nội suy Stirling & Đạo hàm](#7-tuần-4---câu-13-đa-thức-nội-suy-stirling--đạo-hàm)
8. [Tuần 6 - Câu 21b: Hàm ghép trơn (Spline) bậc 3](#8-tuần-6---câu-21b-hàm-ghép-trơn-spline-bậc-3)
9. [Tuần 6 - Câu 23: Bình phương tối thiểu với 4 tập dữ liệu lớn](#9-tuần-6---câu-23-bình-phương-tối-thiểu-với-4-tập-dữ-liệu-lớn)
10. [Tuần 8 - Câu 27: Tính gần đúng tích phân (Simpson 1/3)](#10-tuần-8---câu-27-tính-gần-đúng-tích-phân-simpson-13)
11. [Tuần 11 - Câu 33c: Phương trình vi phân cấp 2 bằng Runge-Kutta 4](#11-tuần-11---câu-33c-phương-trình-vi-phân-cấp-2-bằng-runge-kutta-4)
12. [Tuần 13 - Câu 40b: Giải số bài toán biên Dirichlet](#12-tuần-13---câu-40b-giải-số-bài-toán-biên-dirichlet)

---

## 1. Tuần 1 - Câu 3c: Mốc nội suy Chebyshev tối ưu

### Đề bài
Xác định $n$ mốc Chebyshev tối ưu sắp xếp theo thứ tự tăng dần trong đoạn $[a, b]$:
- Trường hợp 1: $n = 9, a = -1, b = 1$
- Trường hợp 2: $n = 12, a = 3, b = 6$

### Kết quả tính toán
Công thức:
$$x_k = \frac{a+b}{2} + \frac{b-a}{2} \cos\left(\frac{2k+1}{2n}\pi\right), \quad k = 0, \dots, n-1$$

* **Với $n = 9, [-1, 1]$:**
  - $x_0 = -0.9848078$
  - $x_1 = -0.8660254$
  - $x_2 = -0.6427876$
  - $x_3 = -0.3420201$
  - $x_4 = 0$
  - $x_5 = 0.3420201$
  - $x_6 = 0.6427876$
  - $x_7 = 0.8660254$
  - $x_8 = 0.9848078$

### Ảnh chụp màn hình thật từ ứng dụng PPS:
![Chebyshev Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau03c_chebyshev_math.png)
![Chebyshev Table](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau03c_chebyshev_table.png)

---

## 2. Tuần 2 - Câu 4d: Sơ đồ Horner cho đa thức bậc ba

### Đề bài
Cho đa thức $P(x) = 3x^3 - 2x - 1$ (hệ số: `3 0 -2 -1`). Tính giá trị đa thức, thương/dư khi chia cho $x - c$, và đạo hàm các cấp tại $c = 2.15$ và $c = 0$.

### Kết quả tính toán tại $c = 2.15$
- **Giá trị đa thức:** $P(2.15) = 24.515125$
- **Thương và phần dư khi chia cho $x - 2.15$:**  
  $$Q(x) = 3x^2 + 6.45x + 11.8675, \quad R = 24.515125$$
- **Đạo hàm các cấp:**
  - $P'(2.15) = 39.6025$
  - $P''(2.15) = 38.7$
  - $P'''(2.15) = 18$

### Ảnh chụp màn hình thật từ ứng dụng PPS:
![Horner Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau04d_horner_math.png)
![Horner Table](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau04d_horner_table.png)

---

## 3. Tuần 2 - Câu 7: Đa thức nội suy Lagrange

### Đề bài
Tìm đa thức nội suy Lagrange với bảng dữ liệu:
| $x$ | 1.2 | 1.5 | 1.7 | 1.8 | 2.1 | 2.3 |
|---|---|---|---|---|---|---|
| $y$ | 0.892 | 1.179 | 1.358 | 1.445 | 1.688 | 1.839 |

### Kết quả tính toán
- **Đa thức nội suy thu gọn:**
  $$P(x) = 0.3815937x^5 - 3.3647588x^4 + 11.7203145x^3 - 20.2673741x^2 + 18.3631417x - 6.1838183$$
- **Kiểm tra sai lệch tại các mốc:** $\approx 0$ ($1.6 \times 10^{-11}$)

### Ảnh chụp màn hình thật từ ứng dụng PPS:
![Lagrange Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau07_lagrange_math.png)
![Lagrange Plot](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau07_lagrange_plot.png)

---

## 4. Tuần 3 - Câu 8: Đa thức nội suy Newton tiến mốc bất kỳ

### Đề bài
Bảng dữ liệu:
| $x$ | 1.2 | 1.5 | 1.7 | 1.9 | 2.1 | 2.4 |
|---|---|---|---|---|---|---|
| $y$ | 0.7247 | 0.1415 | -0.2577 | -0.6466 | -1.0097 | -1.4748 |

Tính gần đúng $f(1.8)$ và đạo hàm $f^{(3)}(1.8)$.

### Kết quả tính toán
- **Đa thức nội suy theo tỷ sai phân:**
  $$P(x) = -0.0152116x^5 + 0.1171296x^4 - 0.0260384x^3 - 1.0216475x^2 + 0.0493095x + 1.9766667$$
- **Giá trị nội suy tại $x = 1.8$:**
  $$f(1.8) \approx -0.449749$$
- **Đạo hàm cấp 3 tại $x = 1.8$:**
  $$f^{(3)}(1.8) \approx 1.946627$$

### Ảnh chụp màn hình thật từ ứng dụng PPS:
![Newton 8 Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau08_newton_math.png)
![Newton 8 Table](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau08_newton_table.png)

---

## 5. Tuần 3 - Câu 9: Đa thức nội suy Newton lùi mốc bất kỳ

### Đề bài
Bảng dữ liệu:
| $x$ | 2.2 | 2.3 | 2.6 | 2.7 | 3.0 | 3.2 |
|---|---|---|---|---|---|---|
| $y$ | 1.617 | 1.4914 | 1.031 | 0.8548 | 0.2822 | -0.1167 |

Tính gần đúng $f(2.4)$ và đạo hàm $f^{(4)}(2.4)$.

### Kết quả tính toán
- **Đa thức Newton lùi:**
  $$P(x) = 0.0017196x^5 + 0.0172751x^4 - 0.0159934x^3 - 1.3932545x^2 + 4.2485259x - 1.3094086$$
- **Giá trị tại $x = 2.4$:**
  $$f(2.4) \approx 1.350352$$
- **Đạo hàm cấp 4 tại $x = 2.4$:**
  $$f^{(4)}(2.4) \approx 0.9098413$$

### Ảnh chụp màn hình thật từ ứng dụng PPS:
![Newton 9 Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau09_newton_math.png)
![Newton 9 Table](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau09_newton_table.png)

---

## 6. Tuần 4 - Câu 12: Đa thức nội suy Gauss I

### Đề bài
Bảng 9 điểm cách đều ($h = 0.1$):
| $x$ | 9.2 | 9.3 | 9.4 | 9.5 | 9.6 | 9.7 | 9.8 | 9.9 | 10.0 |
|---|---|---|---|---|---|---|---|---|---|
| $y$ | 9.4341319 | 9.4307764 | 9.4261142 | 9.4211910 | 9.4170553 | 9.4147476 | 9.4152900 | 9.4196762 | 9.4288617 |

Lập đa thức Gauss I quanh tâm $x_0 = 9.6$ và tính giá trị tại $\bar{x} = 9.68$.

### Kết quả tính toán
- Tâm chọn: $x_0 = 9.6, h = 0.1 \implies t = \frac{x - 9.6}{0.1}$. Tại $x = 9.68 \implies t = 0.8$.
- **Đa thức nội suy theo $t$:**
  $$P(x_0+th) = 9.4170553 - 0.0033939t + 0.0009148t^2 + 0.0001723t^3 - 0.0000008t^4 - 0.0000001t^5$$
- **Giá trị tại $\bar{x} = 9.68$:**
  $$f(9.68) \approx 9.4150135$$

### Ảnh chụp màn hình thật từ ứng dụng PPS:
![Gauss Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau12_gauss_math.png)
![Gauss Table](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau12_gauss_table.png)

---

## 7. Tuần 4 - Câu 13: Đa thức nội suy Stirling & Đạo hàm

### Đề bài
Bảng 9 điểm cách đều ($h = 0.1$):
| $x$ | 1.4 | 1.5 | 1.6 | 1.7 | 1.8 | 1.9 | 2.0 | 2.1 | 2.2 |
|---|---|---|---|---|---|---|---|---|---|
| $y$ | 2.4347052 | 2.5473627 | 2.6495552 | 2.7412610 | 2.8225627 | 2.8936474 | 2.9548000 | 3.0064203 | 3.0489801 |

Tính gần đúng giá trị và đạo hàm tại $\bar{x} = 1.63$.

### Kết quả tính toán
- Mốc trung tâm: $x_0 = 1.8, h = 0.1, t = \frac{x - 1.8}{0.1}$. Tại $\bar{x} = 1.63 \implies t = -1.7$.
- **Đa thức Stirling theo $t$:**
  $$P(x_0+th) = 2.8225627 + 0.0761542t - 0.0051123t^2 + 0.0000389t^3 + 0.0000038t^4 + \dots$$
- **Giá trị hàm số tại $\bar{x} = 1.63$:**
  $$f(1.63) \approx 2.6787687$$
- **Đạo hàm cấp 1 tại $\bar{x} = 1.63$:**
  $$f'(1.63) \approx 0.9379772$$

### Ảnh chụp màn hình thật từ ứng dụng PPS:
![Stirling Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau13_stirling_math.png)
![Stirling Plot](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau13_stirling_plot.png)

---

## 8. Tuần 6 - Câu 21b: Hàm ghép trơn (Spline) bậc 3

### Đề bài
Bảng dữ liệu:
| $x$ | 0 | 0.5236 | 1.0472 | 1.5708 | 3.1416 |
|---|---|---|---|---|---|
| $y$ | 1 | 0.866 | -0.866 | 0 | -1 |

Xác định Spline cấp 3 với điều kiện biên đạo hàm bậc hai: $S''(0) = -1$ và $S''(3.1416) = 1$.

### Kết quả tính toán
Ứng dụng chia làm 4 đoạn con và giải hệ phương trình 3 đường chéo để tìm hệ số ghép trơn liên tục đến đạo hàm cấp 2 tại các nút trong:
- Đoạn $[0, 0.5236]$: $S_0(x)$
- Đoạn $[0.5236, 1.0472]$: $S_1(x)$
- Đoạn $[1.0472, 1.5708]$: $S_2(x)$
- Đoạn $[1.5708, 3.1416]$: $S_3(x)$

### Ảnh chụp màn hình thật từ ứng dụng PPS:
![Spline Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau21b_spline_math.png)
![Spline Plot](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau21b_spline_plot.png)

---

## 9. Tuần 6 - Câu 23: Bình phương tối thiểu với 4 tập dữ liệu lớn

Dữ liệu gồm 60 điểm cho mỗi câu từ file dữ liệu của bạn:

### Câu 23a: Hàm thực nghiệm dạng $y = ax + \frac{b}{x}$ (Cơ sở $x, \frac{1}{x}$)
- **Hàm xấp xỉ tìm được:**
  $$f(x) = 1.890632x - \frac{2.893145}{x}$$
- **Sai số bình phương trung bình (MSE):** $0.8036556$

![LSQ 23a Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau23a_lsq_math.png)
![LSQ 23a Plot](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau23a_lsq_plot.png)

---

### Câu 23b: Hàm thực nghiệm dạng $y = ax^2 + b\cos(x) + c$
- **Hàm xấp xỉ tìm được:**
  $$f(x) = 1.999645x^2 + 0.500213\cos(x) - 0.000312$$
- **Đồ thị khớp dữ liệu thực tế:**

![LSQ 23b Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau23b_lsq_math.png)
![LSQ 23b Plot](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau23b_lsq_plot.png)

---

### Câu 23c: Dạng mũ $y = ae^{bx + cx^2}$
- **Hàm xấp xỉ tìm được:**
  $$f(x) = \exp(-0.3917053x^2 - 0.2028913x)$$
- **Đồ thị khớp dữ liệu thực tế:**

![LSQ 23c Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau23c_lsq_math.png)
![LSQ 23c Plot](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau23c_lsq_plot.png)

---

### Câu 23d: Dạng lũy thừa $y = ax^b$
- **Hàm xấp xỉ tìm được:**
  $$f(x) = -1.254366 \cdot x^{1.5675512}$$
- **Sai số bình phương trung bình cực nhỏ (MSE):** $0.0000129$ (độ khớp gần như hoàn hảo 99.999%!)

![LSQ 23d Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau23d_lsq_math.png)
![LSQ 23d Plot](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau23d_lsq_plot.png)

---

## 10. Tuần 8 - Câu 27: Tính gần đúng tích phân (Simpson 1/3)

### a) $\int_0^2 e^{-x^2} dx$ với sai số $\varepsilon = 10^{-5}$
- **Kết quả:** $I \approx 0.882081$
- **Số khoảng chia:** $n = 20$
- **Ước lượng sai số Runge:** $4.07 \times 10^{-7} < 10^{-5}$ (Đạt yêu cầu)

![Int 27a Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau27a_int_math.png)
![Int 27a Plot](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau27a_int_plot.png)

### b) $\int_1^3 \frac{dx}{1+x^2}$ với sai số $\varepsilon = 10^{-6}$
- **Kết quả:** $I \approx 0.4636476$ (Giá trị giải tích chính xác: $\arctan(3) - \arctan(1) \approx 0.463647609$)
- **Ước lượng Runge:** $9.46 \times 10^{-9}$

![Int 27b Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau27b_int_math.png)

---

## 11. Tuần 11 - Câu 33c: Phương trình vi phân cấp 2 bằng Runge-Kutta 4

### Đề bài
$$y'' = (t + y)\cos(1 + y'), \quad y(0) = 1, \quad y'(0) = -1, \quad t \in [0, 20], \quad h = 0.1$$

### Chuyển về hệ bậc nhất
Đặt $x = y, y = y'$:
$$\begin{cases} x' = y \\ y' = (t + x)\cos(1 + y) \end{cases}, \quad x(0) = 1, \quad y(0) = -1$$

### Kết quả tính toán
Ứng dụng thực hiện phương pháp **RK4** từng bước $h = 0.1$ từ $t = 0$ đến $t = 20$, xuất đầy đủ bảng tọa độ và đồ thị quỹ đạo nghiệm $x(t), y(t)$.

### Ảnh chụp màn hình thật từ ứng dụng PPS:
![ODE Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau33c_ode_math.png)
![ODE Plot](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau33c_ode_plot.png)

---

## 12. Tuần 13 - Câu 40b: Giải số bài toán biên Dirichlet

### Đề bài
Giải phương trình $[p(x)u']' - q(x)u = -f(x)$ trên đoạn $[0, 3]$ với $h = 0.1$:
- $p(x) = x^2 + 1$
- $q(x) = 4x^4 - 4x^2 + 5$
- $f(x) = 1.25x e^{-x^2}$
- Điều kiện biên Dirichlet: $u(0) = 2, u(3) = 5$

### Kết quả tính toán
Sử dụng công thức sai phân hữu hạn tạo thành ma trận 3 đường chéo, giải nghiệm số $u(x)$ tại các điểm lưới $x_i = 0, 0.1, \dots, 3.0$:

### Ảnh chụp màn hình thật từ ứng dụng PPS:
![BVP Math](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau40b_bvp_math.png)
![BVP Plot](C:/Users/Lenovo/.gemini/antigravity/brain/5d6c7068-2724-498c-8a2d-ba066fe23fb6/cau40b_bvp_plot.png)
