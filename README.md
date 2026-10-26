# PPS_APP — Ứng dụng Tính toán & Phương pháp số

**PPS_APP** là phần mềm tính toán khoa học và giải bài toán số học chuyên sâu viết bằng **Python** kết hợp giao diện đồ họa trực quan (**Tkinter**). Ứng dụng hỗ trợ trình bày lời giải toán học bằng công thức LaTeX sắc nét (MathText), vẽ đồ thị tương tác, lập bảng sai phân / bảng tích chuẩn mực và xuất báo cáo đa định dạng (HTML, Excel, CSV, JSON).

---

## 🌟 Tính năng chính

### 1. Nội suy & Xấp xỉ hàm số
- **Đa thức nội suy**:
  - Lagrange, Newton (mốc bất kỳ qua tỷ sai phân, mốc cách đều tiến / lùi).
  - Nội suy trung tâm: Gauss I, Gauss II, Stirling, Bessel (hỗ trợ ký hiệu chuẩn $t$, $h$, $x_0$).
  - Bảng sai phân so le và bảng tích bậc giảm dần theo đúng chuẩn bài giảng.
- **Nội suy ngược**: Tìm nghiệm theo phương pháp lặp tiến/lùi, tự động kiểm tra điều kiện đơn điệu và hội tụ.
- **Đường cong Spline**: Tuyến tính (Linear), Bậc hai (Quadratic biết đạo hàm biên), Bậc ba (Cubic Natural, Clamped hoặc biết đạo hàm bậc hai).
- **Hồi quy bình phương tối thiểu (Least Squares)**:
  - Đa thức bậc tùy ý, hàm lũy thừa ($y = a x^b$), hàm mũ ($y = a e^{bx}$), hoặc tổ hợp hàm cơ sở tự chọn.
  - Tự động chuẩn hóa/dịch dữ liệu và giải hệ bằng phân tích SVD ổn định.
- **Đa thức Chebyshev & Sơ đồ Horner**:
  - Sinh $n+1$ mốc Chebyshev tối ưu trên đoạn $[a, b]$.
  - Thuật toán Horner tính giá trị đa thức và các cấp đạo hàm.

### 2. Đạo hàm & Tích phân số
- **Đạo hàm số từ bảng mốc**: Lập công thức đạo hàm nhiều cấp (cấp 1 đến $p$) dựa trên hệ số trọng số Lagrange.
- **Tích phân số 1 chiều**:
  - Nhập từ biểu thức giải tích $f(x)$ hoặc bảng số liệu rời rạc $(x, y)$.
  - Các công thức: Hình thang (Trapezoidal), Simpson 1/3, Simpson 3/8, Trung điểm (Midpoint), Newton-Cotes bậc 1 đến 10.
  - Tự động chia lưới theo dung sai mong muốn hoặc ước lượng sai số theo quy tắc Runge / cận đạo hàm $M$.
- **Tích phân kép**: Tính tích phân trên miền chữ nhật $[a, b] \times [c, d]$ bằng phương pháp tích Tensor lưới tự tinh chỉnh.

### 3. Phương trình vi phân (ODE)
- Giải bài toán giá trị ban đầu (Cauchy) cho **một phương trình** hoặc **hệ phương trình vi phân cấp 1**:
  - Phương pháp Euler (Hiện & Ẩn).
  - Phương pháp Hình thang (Trapezoidal ODE).
  - Họ Runge-Kutta: RK2 (Euler cải tiến, Heun), RK3, RK4 cổ điển, RK tham số.
  - Phương pháp đa bước Adams: Adams-Bashforth (dự báo) kết hợp Adams-Moulton (hiệu chỉnh) từ cấp 2 đến cấp 5.
  - Bảng Butcher tùy biến: Cho phép nhập bảng Butcher tổng quát dạng ma trận.

### 4. Bài toán biên (BVP) & Bài toán trị riêng
- **Bài toán biên cấp 2**: Giải phương trình dạng $[p(x)u']' - q(x)u = -f(x)$ với các loại điều kiện biên:
  - Dirichlet ($u$ tại biên)
  - Neumann ($u'$ tại biên)
  - Robin (kết hợp tuyến tính $\alpha u' + \beta u = \gamma$)
  - Sử dụng sai phân hữu hạn đưa về hệ ma trận ba đường chéo (Tridiagonal matrix).
- **Bài toán trị riêng**: Tìm các trị riêng nhỏ nhất và đồ thị hàm riêng tương ứng.

---

## 🚀 Hướng dẫn cài đặt & Khởi chạy

### Yêu cầu hệ thống
- Hệ điều hành: Windows 10/11 (hoặc Linux/macOS có cài đặt môi trường đồ họa Tk).
- Python: Phiên bản **3.11 trở lên**.

### Cách 1: Khởi chạy nhanh trên Windows (Khuyên dùng)
1. **Lần đầu tiên**: Nhấp đúp vào **`CAI_DAT.bat`** để tạo môi trường ảo `.venv` và tự động cài đặt các thư viện cần thiết.
2. **Sử dụng**: Nhấp đúp vào **`CHAY_PPS.bat`** để mở ứng dụng.

### Cách 2: Khởi chạy thủ công bằng dòng lệnh
```bash
# 1. Tạo môi trường ảo (nếu chưa có)
python -m venv .venv

# 2. Kích hoạt môi trường ảo
# Trên Windows:
.venv\Scripts\activate
# Trên Linux/macOS:
source .venv/bin/activate

# 3. Cài đặt các thư viện phụ thuộc
pip install -r requirements.txt

# 4. Khởi chạy ứng dụng
python main.py
```

---

## 🖥️ Hướng dẫn sử dụng nhanh

1. **Chọn phương pháp**: Ở cột danh mục bên trái, chọn nhóm công cụ (Nội suy, Đạo hàm, Tích phân, Vi phân, Biên...) rồi chọn thuật toán cụ thể.
2. **Nhập dữ liệu**:
   - Dán hoặc gõ dãy mốc $x, y$ (cách nhau bởi dấu cách hoặc xuống dòng).
   - Nhập từ file **Excel (.xlsx)** hoặc **CSV** bằng nút **Nhập bảng**.
   - Hỗ trợ nhập hàm số toán học: `sin(x)`, `cos(x)`, `exp(x)`, `ln(x)`, `sqrt(x)`, `x^2`...
3. **Tính toán**: Nhấn nút **Tính kết quả** (hoặc phím tắt **`Ctrl + Enter`**).
4. **Xem & Xuất kết quả**:
   - **Tab Lời giải toán học**: Hiển thị công thức ký hiệu toán LaTeX trực quan, các bước phân tích chi tiết.
   - **Tab Bảng số liệu**: Xem bảng sai phân, bảng tích, bảng bước tính dạng bảng trực quan.
   - **Tab Đồ thị**: Đồ thị trực quan với các công cụ phóng to (zoom), di chuyển (pan), đặt lại góc nhìn và lưu ảnh.
   - **Lưu bài / Mở bài** (`Ctrl + S` / `Ctrl + O`): Lưu toàn bộ thông số đầu vào vào file `.json` để tiếp tục làm việc sau.
   - **Lưu lời giải**: Xuất báo cáo lời giải độc lập ra file HTML (có sẵn công thức KaTeX & hình ảnh) hoặc TXT.
   - **Xuất bảng**: Lưu bảng dữ liệu số ra file Excel (`.xlsx`) hoặc `.csv`.

---

## 📁 Cấu trúc thư mục

```text
PPS_APP/
├── pps/                        # Gói mã nguồn cốt lõi của ứng dụng
│   ├── gui.py                  # Giao diện người dùng Tkinter & điều phối tác vụ
│   ├── interpolation.py        # Thuật toán nội suy, xấp xỉ, spline, bình phương tối thiểu
│   ├── calculus.py             # Tính đạo hàm & tích phân số
│   ├── differential.py         # Phương trình vi phân ODE, bài toán biên BVP & trị riêng
│   ├── common.py               # Xử lý chuỗi biểu thức, bảng sai phân, xác thực dữ liệu
│   └── math_view.py            # Dựng công thức MathText LaTeX và xuất báo cáo HTML
├── examples/                   # Các ví dụ mẫu (.json) cho từng phương pháp
├── tests/                      # Bộ kiểm thử tự động (Unit tests & GUI tests)
├── data/                       # Dữ liệu bài tập đã trích xuất
├── output/                     # Báo cáo, ảnh lời giải và kết quả chạy bài tập
├── CAI_DAT.bat                 # Script cài đặt môi trường tự động
├── CHAY_PPS.bat                # Script chạy nhanh ứng dụng
├── KIEM_TRA.bat                # Script chạy kiểm thử độ chính xác
├── main.py                     # Điểm khởi động ứng dụng
├── requirements.txt            # Danh sách thư viện Python cần thiết
└── README.md                   # Tài liệu giới thiệu dự án
```

---

## 🧪 Kiểm thử (Testing)

Để kiểm tra độ chính xác của toàn bộ các thuật toán và giao diện, chạy script:
```bash
KIEM_TRA.bat
```
Hoặc qua dòng lệnh:
```bash
python -m unittest discover tests
```
Bộ kiểm thử kiểm tra tái lập đa thức chuẩn, các phương pháp tích phân, ODE, cách trình bày và tải/xuất dữ liệu Excel/CSV.
