# Hướng dẫn dùng PPS khi thi ngoại tuyến

## Chuẩn bị trên máy mang vào phòng thi

Chạy `KIEM_TRA.bat` và `CHAY_CONSOLE.bat examples\02_bessel.json` trước khi ngắt mạng; ví dụ Bessel phải in `P^(1)(1.5) = 3`. Các thư viện trong `.venv` phải được cài **trước** khi thi; `CAI_DAT.bat` cần Internet để tải gói nếu máy chưa có sẵn. `CHAY_PPS.bat` mở GUI. Bộ kiểm thử cửa sổ đã chạy qua ngoài sandbox; sandbox kiểm định gây lỗi Tcl `init.tcl`. Giữ `CHAY_CONSOLE.bat` làm đường chạy dự phòng và thử cả hai trên chính máy thi. Terminal gọi một file JSON; xem các file trong `examples/` để biết cấu trúc `tool`, `method`, `values`, `x`, `y`. Đổi **tất cả** dữ liệu và tham số theo đề; các giá trị ví dụ không phải dữ liệu mặc định cho bài thi.

Trong JSON, `x` và `y` là chuỗi số cách nhau bằng dấu cách; dấu thập phân là dấu chấm. Ví dụ:

```json
{"tool":"Nội suy đa thức","method":"Bessel","values":{"target":"0.64","order":"1"},"x":"0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0","y":"0.427 0.264 0.125 0.016 -0.057 -0.088 -0.071 0"}
```

Với `f(x)=x³−2x+1`, tung độ đúng cho `x=0.3,…,1.0` lần lượt là `0.427, 0.264, 0.125, 0.016, -0.057, -0.088, -0.071, 0`. Kết quả kiểm tra: `P(0.64)=-0.017856`, `P′(0.64)=-0.7712`. Nếu đổi ví dụ, hãy tính lại toàn bộ y từ dữ liệu đề.

## Bảng lớn và số mốc của Stirling

Trong **Nội suy đa thức**, ô **Số mốc dùng** điều khiển số mốc thực sự đi vào đa thức:

- **Để trống:** dùng toàn bộ bảng. Giới hạn 100 mốc đã được bỏ. Với bảng lớn, báo cáo giữ hệ số dạng tích và tính bằng Newton lồng với số chữ số trung gian tăng theo bậc và độ chênh thang dữ liệu. Không ép khai triển hệ số theo x khi thao tác đó dễ tràn số.
- **Nhập K:** chọn đúng K mốc gần điểm cần tính. Phải nhập điểm tính; kết quả ghi `K/N` và liệt kê mốc đã chọn. Đây là đa thức của K mốc, có thể khác đa thức toàn bảng.
- **Stirling:** K phải lẻ. Bảng 100 mốc không thể dùng toàn bộ với công thức Stirling đối xứng đang học; nếu đề yêu cầu 9 mốc thì nhập `9`. Bảng 101 mốc có thể để trống để dùng đủ 101 mốc.
- **Bessel:** K phải chẵn và ít nhất 2. **Gauss I/II:** hỗ trợ cả chẵn lẫn lẻ; số chẵn dùng gốc giữa trái cho Gauss I và giữa phải cho Gauss II. Các công thức này đều yêu cầu lưới x cách đều.

Ví dụ ngoại tuyến đã kiểm tra:

```bat
CHAY_CONSOLE.bat examples\08_stirling_chon_9_tu_100_moc.json
CHAY_CONSOLE.bat examples\09_stirling_toan_bo_101_moc.json
```

Ví dụ 08 dùng `y=x²`, chọn 9/100 mốc và trả `P′(49.2)=98.4`. Ví dụ 09 dùng đủ 101 mốc của `y=x³−2x+1`, trả `P(0.25)=0.515625` và `P″(0.25)=1.5`. Trong Python có thể gọi `ip.interpolate(x,y,'Stirling',z,1,count=9)`. Đối tượng `report.value` hỗ trợ gọi giá trị, `.deriv(r)` và `.coef` theo bậc tăng; `.coef` báo lỗi nếu hệ số khai triển vượt miền float, trong khi giá trị tại điểm vẫn có thể tính được.

Dùng bậc cao không tự bảo đảm xấp xỉ tốt hơn hàm thật. Phải giữ số mốc đề yêu cầu; khi đề cho phép tự chọn thì xét khoảng tính, sai số dữ liệu và cận sai số. Tính toàn bảng dùng thời gian bậc hai theo số mốc và cần thêm chữ số trung gian; không có cam kết chạy nhanh với bảng kích thước bất kỳ.

## Tra nhanh theo dạng câu hỏi

| Đề hỏi gì | Gọi hàm nào | Nhập gì | Đọc kết quả thế nào; phần tự làm |
|---|---|---|---|
| Sinh mốc tối ưu | `ip.chebyshev(n,a,b)` | `n` là bậc; hàm trả `n+1` mốc | `report.value` là dãy x. Muốn chứng minh minimax, dùng tính chất Chebyshev trong PDF trang 7–8. |
| Chọn K mốc quanh z | `ip.select_points(x,y,z,K).value` | Toàn bộ bảng và z | Nhận `(xx,yy)`; lưới đều cho nhóm liên tiếp, không đều dùng K điểm gần nhất. Ghi nếu ngoại suy. |
| Lập P và tính P(z) | `ip.interpolate(xx,yy,'Lagrange',z,0)` hoặc `'Newton bất kỳ'` | Mốc x phân biệt, bất kỳ | `report.value` là đa thức theo x; `report.text()` có bảng và hệ số. Mốc đều có thể dùng `'Newton tiến'`/`'Newton lùi'`. |
| Stirling/Gauss/Bessel và đạo hàm | `ip.interpolate(xx,yy,'Bessel',z,1)` | Bessel: số mốc chẵn; Stirling: số mốc lẻ; Gauss hỗ trợ cả hai; x đều | Đọc `P^(r)(z)` và biến `t` trong báo cáo; đạo hàm đã chia `h^r`. Chép đúng gốc t và công thức sai phân. |
| Horner cho đa thức | `ip.horner([a_n,...,a_0],c,r)` | **Horner giao diện dùng bậc giảm**, kể cả hệ số 0 | `report.text()` in P và đạo hàm tới cấp r. Mã giả PDF A02–A04 dùng bậc tăng: đảo danh sách khi chuyển. |
| Cận sai số nội suy tại z | `ch1.interpolation_error_bound(xx,z,M_N)` | M_N là cận hợp lệ của `|f^(N)|`; N là số mốc | `report.value` là cận tại z, không tự là cận toàn đoạn. Với lưới đều dùng `interpolation_global_bound`; với mốc Chebyshev trên `[a,b]` dùng `chebyshev_global_bound(N,a,b,M_N)`. |
| Điền ô y trống | `ch1.fill_missing(x,[y0,None,y2,...],degree)` | Đề **phải cho** bậc ≤ degree và có ít nhất degree+1 điểm đã biết | `report.value` là toàn bộ y; hàm báo mâu thuẫn nếu các điểm đã biết khác không khớp. Không tự suy y từ bảng thiếu giả thiết. |
| Tìm x từ mức y | `ch1.scan_level`, `ch1.monotone_branches`, rồi `ch1.inverse_by_swap` hoặc `ch1.bisect_model` | Chọn nhánh; đảo bảng cần y đơn điệu mẫu; chia đôi cần đa thức thuận và khoảng trái dấu | Đảo bảng trả nghiệm **mô hình ngược**; chia đôi trả nghiệm và khoảng của **mô hình thuận**. Hai kết quả có thể khác. Không kết luận duy nhất cho hàm thật chỉ từ dãy mẫu. |
| Lặp nội suy ngược | `ch1.inverse_controlled(xx,yy,y_target,(a,b),backward,tolerance)` hoặc `ip.inverse(...)` | x đều; chọn J=[a,b] trong nhánh. Với hàm kiểm soát, nhập thêm `residual_tolerance` nếu cần | Nhánh kiểm soát chỉ trả khi chứng nhận `φ(J)⊆J`, `q<1` và cả hai dung sai; cận chỉ cho nghiệm mô hình. Nhánh `ip.inverse` cũ chỉ dừng thực hành. Nếu kiểm tra co thất bại, dùng chia đôi khi có trái dấu. |
| Spline bậc 1–3 | `ip.spline(x,y,kind,left,right,z)` | Chọn đúng kiểu biên theo đề | `report.value` là PPoly. Các kiểu: tuyến tính, bậc hai biên trái/phải, cubic tự nhiên/kẹp/cấp hai/tuần hoàn/not-a-knot. Dùng `ch1.spline_at(pp,z,r)` để kiểm tra đạo hàm ở nút. |
| Spline bậc 4 | `ch1.spline_degree(x,y,4,[s1,s2,s3])` | Ba đạo hàm đầu trái do đề cho | Bảng hệ số theo `u=x-x_i` bậc tăng; nếu biên nằm ở cả hai đầu, gọi `ch1.spline_boundary(x,y,4,conditions)`. Mỗi điều kiện có dạng `('trái',1,giá_trị)` hoặc `([(phía,cấp,hệ_số),...],giá_trị)` cho tổ hợp tuyến tính. Không tự chọn ba điều kiện. |
| Hồi quy cơ sở tuyến tính | `ip.least_squares(x,y,basis=['1','x','cos(x)'],weights=w)` | Cơ sở đúng mô hình; w dương hoặc bỏ qua | `report.text()` có hệ số, ma trận thiết kế, ma trận chuẩn có trọng số, SSE/MSE/RMSE. Dữ liệu cùng x có thể lặp lại vì đó là nhiều quan sát. |
| Mũ/lũy thừa | `ip.least_squares(x,y,'Hàm mũ',basis=['x'])` hoặc `'Lũy thừa'` | y khác 0 và cùng dấu; lũy thừa còn cần x>0 | Đọc tham số A,b đã khôi phục và phân biệt SSE thang log với SSE thang y. |
| Mũ/lũy thừa với mục tiêu SSE thang y | `ch1.original_scale_fit(x,y,'Hàm mũ',basis=['x'])` hoặc `'Lũy thừa'` | Cơ sở đúng mô hình; có thể nhập `initial=[A,b1,...]` | Trả cực tiểu cục bộ tìm được bằng lặp phi tuyến; xem phần dư và thử tham số đầu khác nếu nhạy. Trong GUI chọn `Hàm mũ — SSE thang y` hoặc `Lũy thừa — SSE thang y`. |
| Các mô hình biến đổi khác | `ch1.transformed_fit(x,y,model)` | `model` là `1/(a+b*x)`, `x/(a+b*x)`, `a+b*ln(x)`, `A*b^x` | Kiểm tra mẫu/log; đọc bảng tham số và sai số ở cả hai thang. Với hai dạng nghịch đảo, nếu đề đòi tối thiểu SSE thang y thì cách biến đổi này chưa đáp ứng. |
| Đạo hàm từ bảng | `ca.differentiate(x,y,z,degree=p,order=r)` | `p+1` mốc; p≥r≥1 | Bảng có trọng số và từng đóng góp. Cận sai số đạo hàm cần chứng minh riêng theo chương 2. |
| Tích phân từ bảng/hàm | `ca.integrate(...)` | Chọn phương pháp, a,b,n; nếu bảng thì x,y phải khớp lưới | `report.value` là tích phân xấp xỉ. Runge là ước lượng, không phải cận bảo đảm; cận lý thuyết chỉ khi M hợp lệ. |
| Tích phân đảo cận hoặc cận bằng nhau | `ca.integrate(...,a=b_cũ,b=a_cũ)`; `ca.double_integral(...)` | Có thể đảo một hoặc hai chiều; đoạn bằng 0 trả 0 | Bảng trọng số dùng đoạn tăng dần; dấu kết quả và cận trong công thức theo đúng đầu vào. Khi tự chia lưới, chương trình còn đối chiếu lưới không lồng nhau để giảm bỏ sót dao động. |
| Euler/RK/Adams cho ODE | `de.ode(['y','-x'],[1,0],0,1,0.01,'RK4')` | Biến độc lập t; một phương trình dùng x, hai phương trình dùng x/y, ba dùng x/y/z. h cùng chiều a→b | `t,Y=report.value`, mỗi cột Y là một thành phần. Dung sai chỉ dành cho phép giải bước ẩn; cần đối chiếu h và h/2 để đánh giá sai số toàn cục. |
| RK2/RK3 tham số | `tab=de.custom_rk2('alpha2',c2)` hoặc `de.custom_rk3(c2,c3)`; truyền `custom=tab` vào `de.ode(...,method='RK tùy chỉnh')` | Nhập tham số đúng ký hiệu; các họ RK3 đặc biệt có quy ước hệ số được in trong bảng Butcher | Nếu đề cho hệ số tự do khác quy ước, nhập nguyên bảng A,b,c bằng **RK tùy chỉnh**. |
| Adams s | `de.ode(...,method='Adams 3')` | s từ 2 đến 5; h cố định | Dùng RK4 khởi tạo, AB s giá trị cũ dự báo, AM có thêm giá trị mới hiệu chỉnh. Bậc toàn cục còn phụ thuộc phần khởi tạo. |
| Bài toán biên | `de.bvp('1+x','0','-2-4*x',0,1,0.05,(1,0,0),(1,1,3))` | Ứng dụng giải **(p·u′)′−q·u=−f**. Mỗi biên là `(alpha,beta,gamma)` của **alpha·u′+beta·u=gamma**, đạo hàm theo chiều x tăng ở cả hai đầu | `x,u=report.value`. Ví dụ này có nghiệm u=x². Biên dùng pháp tuyến ngoài phải đổi dấu đạo hàm trước khi nhập. Phần dư nhỏ của hệ sai phân không phải cận sai số so với nghiệm liên tục. |
| Sturm–Liouville | `de.eigen('1','0','1',0,1,0.01,3)` | Giải **−(p·u′)′+q·u=λ·r·u**, hai đầu Dirichlet 0; p tại nửa bước và r tại nút phải dương | `lam,x,V=report.value`. Trả các trị riêng nhỏ nhất, kể cả âm nếu có; cột V là hàm riêng chuẩn hóa max tuyệt đối bằng 1. Xét thêm h/2 khi so với bài toán liên tục. |
| Gauss/Thomas | `ch1.gaussian_elimination(A,b)`; `ch1.thomas(lower,diag,upper,b)` | Thomas yêu cầu `lower[0]=upper[-1]=0`; pivot gần 0 cần bộ giải có chọn trụ | Các bước khử hiển thị sau chuẩn hóa từng hàng; phần dư tính trên hệ gốc. Khi chép đúng mã giả của đề phải ghi thêm thao tác đổi thang nếu dùng bảng của ứng dụng. |

Với hồi quy, tham số rất lớn có thể được giữ dưới dạng `exp(1000)` thay vì `inf`; phải giữ nguyên dạng này khi trình bày. Với spline, gọi `ch1.spline_at(pp,z,r)` nếu cần kiểm tra đạo hàm chung tại nút; đối tượng PPoly thô trả NaN khi gọi ngoài miền. Với các mô hình cần lặp phi tuyến, thông báo hội tụ không chứng minh đó là cực tiểu toàn cục.

Simpson 3/8 được gọi bằng **Newton–Cotes bậc 3**, với n chia hết cho 3. Các ví dụ bổ sung chạy được bằng terminal:

```bat
CHAY_CONSOLE.bat examples\10_bien_robin_p_bien_thien.json
CHAY_CONSOLE.bat examples\11_rk3_moc_trung_hop_le.json
CHAY_CONSOLE.bat examples\12_tich_phan_ham_dao_dong.json
```

Ví dụ 10 có nghiệm liên tục u=x²; sai số lớn nhất của lưới h=0.05 khoảng 0.00086634. Ví dụ 11 dùng họ RK3 có c₂=c₃=2/3, x(1) gần e với sai số dưới 2·10⁻⁶. Ví dụ 12 có tích phân thật bằng 0; phép đối chiếu lưới bổ sung giúp phát hiện trường hợp lưới ban đầu trùng chu kỳ. Không lấy các tham số ví dụ làm mặc định cho dữ liệu đề thi.

Ví dụ Python dùng ngay tại thư mục dự án:

```python
from pps import interpolation as ip
from pps import chapter1_tools as ch1
from pps import calculus as ca, differential as de

x = [0, 1, 2]
y = [0, 1, 4]
nhanh = ch1.inverse_by_swap(x, y, 2).value       # 5/3: mô hình đảo
P = ip.interpolate(x, y, 'Newton bất kỳ').value
nghiem, khoang = ch1.bisect_model(P, 2, 1, 2).value  # gần sqrt(2): mô hình thuận
print(nhanh, nghiem, khoang)
```

## Đối chiếu từng đề tham khảo

| Đề/câu | Chuỗi thao tác và dữ liệu cần nhập | Đầu ra đủ hay còn phải làm |
|---|---|---|
| 20191 PDF đề thống kê câu 1 | Từ các khoảng điểm và tần số, tự lập bảng cộng dồn; chọn mô hình cho hàm phân bố rồi nội suy tại tổng điểm 12 và tìm ngưỡng sao cho số trúng tuyển kỳ vọng `6.5/0.9` nghìn | Code **không có** công cụ riêng cho tần số nhóm, điểm trung bình, hay số học sinh; phải tự đổi trung bình 4 thành tổng 12 và xử lý nhóm điểm/đơn vị. Nội suy phân bố là giả thiết mô hình, không phải giá trị chính xác từ tần số nhóm. |
| 20191 PDF đề thống kê câu 2 | File `testdata_ab.text` → chọn dạng hàm từ đồ thị/phần dư → `least_squares` hoặc `transformed_fit` → MSE thang đề yêu cầu | File dữ liệu chưa có. Công thức hệ chuẩn, lựa chọn mô hình và mã giả phải tự viết. |
| 20191 PDF đề tầm nhìn câu 1 | File `Kip2Cau1.txt` → khớp mô hình thích hợp hoặc nội suy gần 17 m → quét nhánh và đảo bảng/chia đôi để tìm x tại 20 km | File chưa có. Thuật toán tìm độ cao y bất kỳ và điều kiện một nghiệm phải trình bày; dữ liệu thực nghiệm không tự bảo đảm mô hình. |
| 20191 PDF đề tầm nhìn câu 2 | Ba mốc `(0,1),(π/2,0),(π,-1)` cho `cos(x)`; nếu cho thêm 3 điều kiện biên thì dùng `spline_degree(...,4,...)` hoặc `spline_boundary` | Đề ảnh/PDF không ghi điều kiện biên bậc 4. Với spline C³, còn thiếu 3 điều kiện độc lập; không có một đáp án duy nhất nếu không nêu thêm giả thiết. |
| 20211 câu 1 | Hai bảng năm→nhiệt độ và nhiệt độ→mực nước → nội suy/khớp từng bảng, tính `T` theo STT, lấy đạo hàm theo quy tắc dây chuyền, lập dự báo 2021–2031 | Hai bảng gốc và STT chưa có. Phải chọn mô hình, nêu nguy cơ ngoại suy và tự trình bày quy tắc dây chuyền `dH/dt=(dH/dθ)(dθ/dt)`. |
| 20211 câu 2 | Đọc so sánh ở PDF chương 1 trang 12, 18, 23–26; dùng bảng CSV kiểm tra điều kiện từng công thức | Câu lý thuyết phải viết bằng lời; code không thay thế phần ưu/nhược điểm. |
| 20221 a–b | Có 9 y phụ thuộc `m`; sinh 9 mốc Chebyshev trên `[1,5]`, cộng `1/(1+m)` vào mọi y rồi nội suy Newton bất kỳ/Lagrange | Cần biết m cá nhân. Chebyshev không cách đều nên không dùng Newton tiến/Bessel. Phải trình bày thuật toán chọn mốc. |
| 20221 c | Tạo 9 x đều trên `[-1,1]`, ghép đúng y đã hiệu chỉnh, quét mức `-1.25`, tách từng nhánh, giải riêng | Y có vùng giảm rồi tăng nên không dùng `inverse()` cho toàn bảng. Có thể có nhiều nghiệm; phân nhánh là bắt buộc. |
| 20221 d | `ca.differentiate` với `degree=3`, `order=1` trên cửa sổ 4 điểm ở từng mốc | Trọng số cho biên khác mốc trong. Công thức sai số đạo hàm phải viết riêng. |
| 20222 a–b | Từ số phần tử Y sinh mốc Chebyshev trên `[2,7]`; riêng câu b dùng mốc đều `[2,5]`, chọn 6 mốc quanh `3.ab`, nội suy Bessel/Newton và lấy đạo hàm | File Y và `ab` chưa có. Câu a và b dùng **hai tập mốc khác nhau**; không lấy Y câu a ghép với x câu b mà quên giả thiết của đề. |
| 20222 c–d | Quét nghiệm mức `3.ab` trên dữ liệu câu b, chọn nhánh/chia đôi; ưu nhược điểm các công thức trình bày bằng lời | Dữ liệu chưa có; không suy một nghiệm duy nhất từ bảng nếu chưa có giả thiết. |
| 20231 a–b | File dữ liệu + k → chọn 8 mốc đều gần `5+k·10⁻²`, dựng Bessel, đạo hàm; quét mức `(-1)^k(3+k·10⁻²)` | File đề riêng và k chưa có. Nếu 8 mốc không đều thì Bessel không hợp lệ; phải chuyển Lagrange/Newton tỷ sai phân. |
| 20231 c–d | 5 điểm để suy công thức đạo hàm; Simpson và Newton–Cotes 5 khoảng cho `f(x)ln(2x+k)` trên `[1,6.08]` | `ca.differentiate` và `ca.integrate` tính số khi đã có mẫu hàm trên lưới. Cận sai số đạo hàm/tích phân cần giả thiết và tài liệu chương 2; không lấy phần dư nội suy tại mốc làm cận. |
| 20241 a | `ip.horner([3.5,-2.78,5.63,0,-15.28,1.54],2.57,3)` | `P″=1054.696646`, `P‴=1249.3386`. Bảng có từng lượt chia Horner; khi viết bài cần nêu số dư cấp k nhân `k!`. |
| 20241 b–d | Chọn 9 mốc đều quanh `x̄`; dựng Stirling; dùng `inverse_controlled` với J phù hợp hoặc chia đôi đa thức thuận cho mức `ȳ` | File `20241GKtest2` và giá trị `x̄,ȳ` chưa có. Phải nêu đoạn J và kiểm tra co; nếu không có chứng nhận thì dùng chia đôi hoặc ghi chỉ dừng thực hành. |

## Quy tắc kiểm tra trước khi chép đáp án

- `N` mốc cho đa thức bậc **không quá** `N−1`; spline có `N−1` đoạn. Giữ nguyên đơn vị của x, y và bước h.
- Xem `P(x_i)` có khớp y_i, nhưng nhớ sai lệch tại mốc bằng 0 **không phải** cận sai số so với hàm thật giữa các mốc.
- Hệ số trong `numpy.polynomial.Polynomial.coef` theo **bậc tăng**; giao diện Horner hiện nhập **bậc giảm**. Với dạng `t=(x−x_c)/h`, đạo hàm cấp r theo x có hệ số `h⁻ʳ`.
- Mốc trùng phải xử lý ở dữ liệu gốc. Mô hình log cần kiểm tra dấu và miền; spline bậc p cần `p−1` điều kiện bổ sung. Ngoại suy phải ghi rõ.
- Khi đề hỏi thuật toán hoặc chứng minh, chép đầu vào, đầu ra, công thức chỉ số, vòng lặp, điều kiện dừng và giả thiết từ PDF. Kết quả chạy code chỉ là phần số.
