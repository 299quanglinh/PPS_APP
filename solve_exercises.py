import sys
sys.stdout.reconfigure(encoding='utf-8')
import json
import math
import struct
import ctypes
from ctypes import wintypes
from pathlib import Path
from PIL import Image
import numpy as np

from pps.gui import App, dispatch

user32 = ctypes.windll.user32
gdi32 = ctypes.windll.gdi32

def capture_window(win, filename):
    win.update_idletasks()
    win.update()
    hwnd = user32.GetAncestor(win.winfo_id(), 2)
    rect = wintypes.RECT()
    user32.GetWindowRect(hwnd, ctypes.byref(rect))
    w = rect.right - rect.left
    h = rect.bottom - rect.top

    hwndDC = user32.GetWindowDC(hwnd)
    mfcDC = gdi32.CreateCompatibleDC(hwndDC)
    saveBitMap = gdi32.CreateCompatibleBitmap(hwndDC, w, h)
    gdi32.SelectObject(mfcDC, saveBitMap)

    user32.PrintWindow(hwnd, mfcDC, 2)

    bmpinfo = bytearray(40)
    struct.pack_into('<IiiHHIIIIII', bmpinfo, 0, 40, w, -h, 1, 32, 0, 0, 0, 0, 0, 0)
    buf = bytearray(w * h * 4)
    gdi32.GetDIBits(mfcDC, saveBitMap, 0, h, (ctypes.c_char * len(buf)).from_buffer(buf), (ctypes.c_char * len(bmpinfo)).from_buffer(bmpinfo), 0)

    img = Image.frombuffer('RGBA', (w, h), bytes(buf), 'raw', 'BGRA', 0, 1)
    img.convert('RGB').save(filename)

    gdi32.DeleteObject(saveBitMap)
    gdi32.DeleteDC(mfcDC)
    user32.ReleaseDC(hwnd, hwndDC)

# Load extracted dataset
with open('data/extracted_data.json', 'r', encoding='utf-8') as f:
    datasets = json.load(f)

app = App()
app.geometry('1280x820')
app.update()

results = {}

# Helper to run and capture
def run_and_capture(tool, method, values, x_str, y_str, math_cap=None, table_cap=None, plot_cap=None):
    app.select_tool(tool)
    app.method.set(method)
    app.show_form()
    for k, v in values.items():
        if k in app.form:
            app.form[k].set(str(v))
    if hasattr(app, 'xinput') and app.xinput.winfo_exists():
        app.xinput.delete('1.0', 'end'); app.xinput.insert('1.0', x_str)
    if hasattr(app, 'yinput') and app.yinput.winfo_exists():
        app.yinput.delete('1.0', 'end'); app.yinput.insert('1.0', y_str)
    
    rep = dispatch(tool, method, {**{k: v.get() for k, v in app.form.items()}, **values}, x_str, y_str)
    app.display(rep)
    
    if math_cap:
        app.notebook.select(0)
        capture_window(app, math_cap)
    if table_cap:
        app.notebook.select(2)
        app.show_table()
        capture_window(app, table_cap)
    if plot_cap:
        app.notebook.select(3)
        capture_window(app, plot_cap)
        
    return rep

# ==========================================
# 1. Tuần 1 - Câu 3c: Chebyshev
# ==========================================
print("1. Chạy Câu 3c: Chebyshev...")
r_cheb1 = run_and_capture('Chọn mốc', 'Chebyshev', {'degree': '9', 'a': '-1', 'b': '1'}, '', '',
                          math_cap='output/cau03c_chebyshev_math.png',
                          table_cap='output/cau03c_chebyshev_table.png')
results['cau3c_1'] = r_cheb1.text

r_cheb2 = run_and_capture('Chọn mốc', 'Chebyshev', {'degree': '12', 'a': '3', 'b': '6'}, '', '')
results['cau3c_2'] = r_cheb2.text

# ==========================================
# 2. Tuần 2 - Câu 4d: Horner cho P(x) = 3x^3 - 2x - 1 tại c = 2.15
# ==========================================
print("2. Chạy Câu 4d: Horner...")
r_horner = run_and_capture('Đạo hàm & Horner', 'Horner', {'coef': '3 0 -2 -1', 'target': '2.15', 'order': '3'}, '', '',
                           math_cap='output/cau04d_horner_math.png',
                           table_cap='output/cau04d_horner_table.png')
results['cau4d'] = r_horner.text

# ==========================================
# 3. Tuần 2 - Câu 7: Lagrange
# ==========================================
print("3. Chạy Câu 7: Lagrange...")
x7 = "1.2 1.5 1.7 1.8 2.1 2.3"
y7 = "0.892 1.179 1.358 1.445 1.688 1.839"
r_lagrange = run_and_capture('Nội suy đa thức', 'Lagrange', {'target': '1.8', 'order': '0'}, x7, y7,
                             math_cap='output/cau07_lagrange_math.png',
                             plot_cap='output/cau07_lagrange_plot.png')
results['cau7'] = r_lagrange.text

# ==========================================
# 4. Tuần 3 - Câu 8: Newton tiến mốc bất kỳ
# ==========================================
print("4. Chạy Câu 8: Newton tiến mốc bất kỳ...")
x8 = "1.2 1.5 1.7 1.9 2.1 2.4"
y8 = "0.7247 0.1415 -0.2577 -0.6466 -1.0097 -1.4748"
r_newton8 = run_and_capture('Nội suy đa thức', 'Newton bất kỳ', {'target': '1.8', 'order': '3'}, x8, y8,
                            math_cap='output/cau08_newton_math.png',
                            table_cap='output/cau08_newton_table.png')
results['cau8'] = r_newton8.text

# ==========================================
# 5. Tuần 3 - Câu 9: Newton lùi mốc bất kỳ
# ==========================================
print("5. Chạy Câu 9: Newton lùi mốc bất kỳ...")
x9 = "2.2 2.3 2.6 2.7 3 3.2"
y9 = "1.617 1.4914 1.031 0.8548 0.2822 -0.1167"
r_newton9 = run_and_capture('Newton', 'Mốc bất kỳ', {'target': '2.4', 'order': '4', 'direction': 'Lùi'}, x9, y9,
                            math_cap='output/cau09_newton_math.png',
                            table_cap='output/cau09_newton_table.png')
results['cau9'] = r_newton9.text

# ==========================================
# 6. Tuần 4 - Câu 12: Gauss 1
# ==========================================
print("6. Chạy Câu 12: Gauss 1...")
x12 = "9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10"
y12 = "9.4341319 9.4307764 9.4261142 9.421191 9.4170553 9.4147476 9.41529 9.4196762 9.4288617"
r_gauss12 = run_and_capture('Nội suy đa thức', 'Gauss I', {'target': '9.68', 'order': '0'}, x12, y12,
                            math_cap='output/cau12_gauss_math.png',
                            table_cap='output/cau12_gauss_table.png')
results['cau12'] = r_gauss12.text

# ==========================================
# 7. Tuần 4 - Câu 13: Stirling
# ==========================================
print("7. Chạy Câu 13: Stirling...")
x13 = "1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2"
y13 = "2.4347052 2.5473627 2.6495552 2.741261 2.8225627 2.8936474 2.9548 3.0064203 3.0489801"
r_stirling13 = run_and_capture('Nội suy đa thức', 'Stirling', {'target': '1.63', 'order': '1'}, x13, y13,
                               math_cap='output/cau13_stirling_math.png',
                               plot_cap='output/cau13_stirling_plot.png')
results['cau13'] = r_stirling13.text

# ==========================================
# 8. Tuần 6 - Câu 21: Spline cấp 3
# ==========================================
print("8. Chạy Câu 21: Spline...")
# Sắp xếp các điểm: (0, 1), (0.5236, 0.866), (1.0472, -0.866), (1.5708, 0), (3.1416, -1)
x21 = "0 0.5236 1.0472 1.5708 3.1416"
y21 = "1 0.866 -0.866 0 -1"
# Câu 21b: Cubic đạo hàm bậc 2 với S''(0) = -1, S''(pi) = 1
r_spline21b = run_and_capture('Spline', 'Cubic đạo hàm bậc hai', {'left': '-1', 'right': '1', 'target': '1.0'}, x21, y21,
                              math_cap='output/cau21b_spline_math.png',
                              plot_cap='output/cau21b_spline_plot.png')
results['cau21b'] = r_spline21b.text

# Câu 21c: Cubic kẹp S'(0) = 0, S'(pi) = 0
r_spline21c = run_and_capture('Spline', 'Cubic kẹp', {'left': '0', 'right': '0', 'target': '1.0'}, x21, y21)
results['cau21c'] = r_spline21c.text

# ==========================================
# 9. Tuần 6 - Câu 23: Bình phương tối thiểu
# ==========================================
print("9. Chạy Câu 23: Least Squares...")
# a) y = ax + b/x -> basis: x; 1/x
xa_str = ' '.join(map(str, datasets['cau23_a']['x']))
ya_str = ' '.join(map(str, datasets['cau23_a']['y']))
r_lsq23a = run_and_capture('Bình phương tối thiểu', 'Tuyến tính theo cơ sở', {'basis': 'x; 1/x'}, xa_str, ya_str,
                           math_cap='output/cau23a_lsq_math.png',
                           plot_cap='output/cau23a_lsq_plot.png')
results['cau23a'] = r_lsq23a.text

# b) y = ax^2 + b*cos(x) + c -> basis: x^2; cos(x); 1
xb_str = ' '.join(map(str, datasets['cau23_b']['x']))
yb_str = ' '.join(map(str, datasets['cau23_b']['y']))
r_lsq23b = run_and_capture('Bình phương tối thiểu', 'Tuyến tính theo cơ sở', {'basis': 'x^2; cos(x); 1'}, xb_str, yb_str,
                           math_cap='output/cau23b_lsq_math.png',
                           plot_cap='output/cau23b_lsq_plot.png')
results['cau23b'] = r_lsq23b.text

# c) y = ae^(bx+cx^2) -> Hàm mũ với basis x; x^2
xc_str = ' '.join(map(str, datasets['cau23_c']['x']))
yc_str = ' '.join(map(str, datasets['cau23_c']['y']))
r_lsq23c = run_and_capture('Bình phương tối thiểu', 'Hàm mũ', {'basis': 'x; x^2'}, xc_str, yc_str,
                           math_cap='output/cau23c_lsq_math.png',
                           plot_cap='output/cau23c_lsq_plot.png')
results['cau23c'] = r_lsq23c.text

# d) y = ax^b -> Lũy thừa
xd_str = ' '.join(map(str, datasets['cau23_d']['x']))
yd_str = ' '.join(map(str, datasets['cau23_d']['y']))
r_lsq23d = run_and_capture('Bình phương tối thiểu', 'Lũy thừa', {}, xd_str, yd_str,
                           math_cap='output/cau23d_lsq_math.png',
                           plot_cap='output/cau23d_lsq_plot.png')
results['cau23d'] = r_lsq23d.text

# ==========================================
# 10. Tuần 8 - Câu 27: Tích phân số
# ==========================================
print("10. Chạy Câu 27: Tích phân...")
# a) e^(-x^2) trên [0, 2], eps = 1e-5
r_int27a = run_and_capture('Tích phân', 'Simpson', {'expr': 'exp(-x^2)', 'a': '0', 'b': '2', 'n': '10', 'tol': '0.00001', 'source': 'hàm'}, '', '',
                           math_cap='output/cau27a_int_math.png',
                           plot_cap='output/cau27a_int_plot.png')
results['cau27a'] = r_int27a.text

# b) 1/(1+x^2) trên [1, 3], eps = 1e-6
r_int27b = run_and_capture('Tích phân', 'Simpson', {'expr': '1/(1+x^2)', 'a': '1', 'b': '3', 'n': '10', 'tol': '0.000001', 'source': 'hàm'}, '', '',
                           math_cap='output/cau27b_int_math.png')
results['cau27b'] = r_int27b.text

# c) 1/ln(ln(2x)) trên [2, 5], eps = 1e-7
r_int27c = run_and_capture('Tích phân', 'Simpson', {'expr': '1/ln(ln(2*x))', 'a': '2', 'b': '5', 'n': '10', 'tol': '0.0000001', 'source': 'hàm'}, '', '')
results['cau27c'] = r_int27c.text

# d) 1/sqrt(x^2+3x+1) trên [1, 3], eps = 5e-6
r_int27d = run_and_capture('Tích phân', 'Simpson', {'expr': '1/sqrt(x^2+3*x+1)', 'a': '1', 'b': '3', 'n': '10', 'tol': '0.000005', 'source': 'hàm'}, '', '')
results['cau27d'] = r_int27d.text

# ==========================================
# 11. Tuần 11 - Câu 33c: Phương trình vi phân cấp 2
# y'' = (t+y)*cos(1+y'), y(0)=1, y'(0)=-1, t in [0, 20], h=0.1
# Hệ: y; (t+x)*cos(1+y)
# ==========================================
print("11. Chạy Câu 33c: ODE...")
r_ode33c = run_and_capture('Phương trình vi phân', 'RK4', {'expr': 'y; (t+x)*cos(1+y)', 'initial': '1 -1', 'a': '0', 'b': '20', 'h': '0.1', 'tol': '1e-6'}, '', '',
                           math_cap='output/cau33c_ode_math.png',
                           plot_cap='output/cau33c_ode_plot.png')
results['cau33c'] = r_ode33c.text

# ==========================================
# 12. Tuần 13 - Câu 40b: Bài toán biên Dirichlet
# p = x^2+1, q = 4*x^4-4*x^2+5, f = 1.25*x*exp(-x^2), a=0, b=3, h=0.1, left: 0 1 2, right: 0 1 5
# ==========================================
print("12. Chạy Câu 40b: BVP Dirichlet...")
r_bvp40b = run_and_capture('Bài toán biên', 'Biên tổng quát', {'p': 'x^2+1', 'q': '4*x^4-4*x^2+5', 'f': '1.25*x*exp(-x^2)', 'a': '0', 'b': '3', 'h': '0.1', 'left': '0 1 2', 'right': '0 1 5'}, '', '',
                           math_cap='output/cau40b_bvp_math.png',
                           plot_cap='output/cau40b_bvp_plot.png')
results['cau40b'] = r_bvp40b.text

app.destroy()

with open('output/solved_results.json', 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)

print("\n>>> ALL CALCULATIONS AND SCREENSHOTS CAPTURED SUCCESSFULLY! <<<")
