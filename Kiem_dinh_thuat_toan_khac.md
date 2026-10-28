# Kiểm định các thuật toán còn lại — cập nhật 28/10/2026

Đã đọc lại thân hàm và tuyến gọi trong `common.py`, `interpolation.py`, `chapter1_tools.py`, `calculus.py`, `differential.py` và phần trình bày spline. Việc kiểm tra này gồm cả công thức, chỉ số, điều kiện áp dụng, đổi thang, tiêu chuẩn dừng và cách báo kết quả. Các bộ giải bổ sung nằm trong `linear_algebra.py` và `inverse_model.py`; không thêm thư viện cần cài.

Tài liệu môn học hiện có vẫn chỉ đầy đủ Chương 1. Các thuật toán Chương 2–3 dưới đây được kiểm chứng bằng tính chất đại số, nghiệm chế tạo và nghiệm giải tích; kết quả đó chưa chứng minh rằng ứng dụng bao phủ đủ giáo trình chưa được cung cấp.

## Lỗi tái hiện và sửa đổi

| Nhóm | Đầu vào / hiện tượng trước sửa | Sửa đổi và kết quả kiểm tra |
|---|---|---|
| Biên Robin với p biến thiên | `(p u′)′=-f`, `p=1+x`, `f=-2-4x`, `u′(0)=0`, `u′(1)+u(1)=3`. Nghiệm thật `u=x²`. Sai số lớn nhất ở h=0.1/0.05/0.025 là `0.05128205 / 0.02531646 / 0.01257862`, chỉ giảm khoảng 2 lần | Trước đây dùng p tại nửa bước cho cả đạo hàm ở biên. Đổi sang p **tại biên** trong hạng Robin. Sai số còn `0.0034641768 / 0.0008663364 / 0.0002166024`, giảm `3.99865 / 3.99966` lần, phù hợp bậc hai |
| Trị riêng với trọng số khác thang | p=1, q=0, r=10⁻²⁰⁰ làm tích r_i r_(i+1) hụt số; r=10²⁰⁰ làm tích này tràn. Ngay cả sau sửa căn tích, LAPACK vẫn có thể coi các phần tử rất nhỏ là 0 hoặc không hội tụ | Lập phép đổi đối xứng bằng Decimal, chuẩn hóa toàn ma trận trước LAPACK, đổi lại đơn vị trị riêng. Với h=0.1, ba giá trị `λ*r` ở cả hai thang đều là `9.78869674, 38.19660113, 82.44294954`, đúng phổ sai phân giải tích |
| Bước ODE ẩn vô nghiệm | Euler ẩn cho x′=x, x(0)=10⁻¹², h=1: phương trình bước là `X−10⁻¹²−X=0`, vô nghiệm. Ngưỡng phần dư tối thiểu 10⁻⁹ trước đây vẫn nhận kết quả | Dùng Jacobian giải tích, ngưỡng phần dư theo dung sai/độ chính xác máy và kiểm tra suy biến Jacobian. Ca này báo lỗi; các bước Euler ẩn/hình thang cho x′=−100x vẫn khớp đúng nghiệm của phương trình sai phân |
| Biến trạng thái và lưới thời gian | Nhập một phương trình x′=y nhưng chỉ có x(0): y từng bị gán 0 ngầm. Với t từ 10¹⁶ đến 10¹⁶+4, h=1, một số nút thời gian trùng do float | Chỉ nhận đúng các biến của hệ đã khai báo; từ chối lưới có nút trùng, yêu cầu đổi gốc/đơn vị |
| Hai họ RK3 đặc biệt | Công thức tổng quát từ chối `(c₂,c₃)=(2/3,2/3)` và `(2/3,0)`, dù tồn tại bộ hệ số RK3 hợp lệ | Bổ sung hai quy ước hệ số và in bảng Butcher thực sự. Kiểm tra độc lập `bᵀc=1/2`, `bᵀc²=1/3`, `bᵀAc=1/6`, cùng sai số nghiệm x′=x. Đây là lựa chọn trong họ có tham số tự do, không khẳng định hai c xác định duy nhất toàn bộ bảng |
| Nội suy ngược nhiều mốc | Nhánh có cận co vẫn chia sai phân cho giai thừa bằng float; hơn 171 mốc có thể tràn dù dữ liệu tuyến tính. Nhánh lặp có thể từ chối J hợp lệ vì ước lượng đầu nằm ngoài J | Hai nhánh dùng chung dạng sai phân Decimal, không còn tràn giai thừa. Nếu ánh xạ đã được kiểm tra trên J, có thể chọn tâm J làm điểm đầu. Đã thử đủ 201 mốc y=2x+1, cả tiến/lùi và hai nhánh lặp; mức 1.5 trả x=0.25 |
| Spline bậc 4 và đơn vị x nhỏ | x=[0,10⁻⁸,2·10⁻⁸], y=x, ba đạo hàm trái [1,0,0] xác định duy nhất spline, nhưng xét hạng của hệ chưa đổi thang báo thiếu điều kiện | Dùng tọa độ riêng `(x−x_i)/h_i`, chuẩn hóa phương trình và giải hệ thưa. Trả đường thẳng và đạo hàm 1; không đổi ba điều kiện đầu vào |
| Spline bảng lớn và sai số truyền đạo hàm | Bộ giải cubic thư viện có thể xử lý bảng lớn nhưng phần trình bày lại dựng ma trận đặc và giải lại. Hệ float cũng có thể báo suy biến trên bảng dài dù đã biết đủ đạo hàm ở cùng một đầu | Bảng lớn xuất hệ số từng đoạn, không dựng hệ đặc chỉ để hiển thị. Khi các đạo hàm đều ở cùng đầu, dùng truy hồi Decimal, đổi t=−x cho đầu phải. Đã thử 1001 mốc với cubic tự nhiên, spline bậc 4 cả hai hướng và spline bậc 2 biên tổng quát; kiểm tra thêm đa thức bậc 1–4 với biên đầu phải |
| Đạo hàm spline tại nút | Spline tuyến tính qua (0,0),(1,0),(2,10⁻¹²) có đạo hàm trái 0 và phải 10⁻¹² ở x=1. Ngưỡng tuyệt đối 10⁻¹⁰ từng coi hai phía bằng nhau | Đánh giá hai phía bằng Decimal; dung sai xét theo độ lớn đạo hàm và các hạng của đa thức, không dùng một ngưỡng tuyệt đối che mất bước nhảy nhỏ. Ca này báo không có đạo hàm chung |
| Hồi quy có cột gần phụ thuộc do đơn vị | Cơ sở [1,x] với x=10¹²+[0,1,2,3], y=[1,2,3,4] có hạng đầy đủ nhưng SVD trên ma trận chưa đổi thang có thể báo phụ thuộc | Chuẩn hóa cột/trọng số; dùng QR Decimal khi cần. Dự đoán khớp [1,2,3,4] và trả 1.5 tại 10¹²+0.5. Cơ sở phụ thuộc thật [1,2] vẫn bị từ chối |
| Ma trận chuẩn và chỉ tiêu hồi quy tràn | x=[0,10¹⁰,2·10¹⁰], w=[10³⁰⁰,2·10³⁰⁰,10³⁰⁰] có phần tử ΘᵀWΘ bằng 6·10³²⁰, vượt float dù bài toán khớp và dự đoán hữu hạn | Chia trọng số cho cùng một hằng khi giải, giữ nguyên mục tiêu cực tiểu. Ma trận chuẩn, SSE và các tổng báo cáo dùng Decimal; giá trị quá lớn được xuất chuỗi khoa học. Kết quả y=[1,2,2] vẫn có hệ số 1.25 và 0.5·10⁻¹⁰ |
| Mô hình mũ có tham số riêng lẻ rất lớn | y=2·exp(1000x), x=[0,0.001,0.002,0.003]. Viết A·bˣ làm b=exp(1000) tràn dù mọi giá trị dữ liệu hữu hạn | Dự đoán trên thang log của cả tích; tham số b được giữ dạng `exp(1000)` khi cần. Không báo `inf` như một tham số số thông thường |
| Tối ưu SSE thang y và hình dạng dự đoán | Dự đoán trên mảng 2×2 lỗi kích thước. Đổi y sang đơn vị 10⁻²⁰⁰ có thể làm bộ giải dừng ngay vì mục tiêu tuyệt đối quá nhỏ | Giữ hình dạng mảng đầu vào; chuẩn hóa y khi tối ưu rồi đổi về đơn vị gốc. Đã kiểm tra cùng một nghiệm khớp sau đổi đơn vị. Dữ liệu y=0 được ghi rõ không xác định duy nhất các hệ số số mũ |
| Tích phân có trung gian tràn/hụt số | ∫₀^(10⁻³⁰⁸) 10³⁰⁸ dx = 1 nhưng tổng trước khi nhân h có thể tràn. Cận Simpson với L=10⁻⁷⁰, n=2, M₄=10³⁰⁰ bị tính thành 0 do nhân sai thứ tự | Cộng có đổi thang/Decimal; kết quả tích phân là 1. Cận đúng còn `10⁻⁵⁰/(180·16)≈3.4722222·10⁻⁵⁴`, không bị hiển thị thành 0. Chọn n từ cận bằng tìm kiếm nguyên, tránh tràn khi lấy căn của một tích rất lớn |
| Dừng Runge khi lấy mẫu trùng chu kỳ | Simpson với cos(40πx) trên [0,1], n=10 rồi 20 đều có thể thấy các mẫu bằng 1, báo chênh lệch gần 0 dù tích phân thật bằng 0 | Trước khi tự dừng, đối chiếu thêm lưới không lồng nhau. Ca này cho kết quả gần 0 ở cả tích phân một chiều và tích phân kép. Một kiểm tra bổ sung vẫn không thể loại mọi dạng trùng chu kỳ có thể được tạo ra |
| Đảo cận và đoạn bằng 0 | Các quy tắc từng từ chối a>b hoặc a=b dù định nghĩa tích phân xử lý được | Hỗ trợ đổi dấu khi đảo cận và trả 0 khi độ dài một chiều bằng 0. Giữ điều kiện Simpson cần số khoảng chẵn và Newton–Cotes cần n chia hết cho bậc |
| Gauss/Thomas với đơn vị khác nhau giữa hàng | Hệ diag(10⁻¹⁰⁰,10¹⁰⁰)·u=[2·10⁻¹⁰⁰,3·10¹⁰⁰] từng bị coi suy biến. Tính phần dư của hệ có hạng 10³⁰⁸ cũng có thể tràn trung gian | Đổi thang từng phương trình trước khử, tính phần dư hệ gốc bằng Decimal. Thu được [2,3]; thêm kiểm tra hệ có triệt tiêu các tích lớn và hệ ba đường chéo/ngẫu nhiên với nghiệm được chế tạo trước |
| Điền ô thiếu ở thang y nhỏ | Với bậc ≤1, các giá trị đã biết [0,10⁻¹²,5·10⁻¹²] mâu thuẫn với đường thẳng nhưng ngưỡng dựa trên max(1,|y|) từng bỏ qua | Kiểm tra theo thang y thực, dùng đa thức Newton Decimal cho các điểm đã biết. Mâu thuẫn bị từ chối; không tự đổi bậc hoặc sửa dữ liệu |

## Danh mục kiểm tra theo thuật toán

Các kiểm thử mới nằm trong `tests/test_remaining_algorithms.py`. Các kiểm thử Chương 1, nội suy đa thức và lỗi biên của đợt trước được chạy lại.

| Thuật toán / biến thể | Bằng chứng độc lập chính |
|---|---|
| Chebyshev, Horner, chia Horner, Newton lồng/thêm mốc, Lagrange/Newton/Gauss/Stirling/Bessel | Chạy lại các bộ `test_numerics`, `test_chapter1_audit`, `test_edge_cases`; đa thức biết trước, đạo hàm giải tích, dữ liệu lớn/biên và tổng Lagrange hữu tỉ độc lập |
| Quét mức và nhánh mẫu đơn điệu | Mức đúng tại mốc, đoạn phẳng, khoảng đổi dấu; không gán tính đơn điệu của hàm thật từ mẫu |
| Đảo bảng, lặp tiến/lùi, lặp có cận co, chia đôi | Hàm tăng/giảm, hai hướng, J nhỏ, 201 mốc tuyến tính, nghiệm tại đầu đoạn, nghiệm √2 và khoảng bao nghiệm |
| Spline tuyến tính, bậc hai trái/phải, cubic tự nhiên/kẹp/đạo hàm hai/tuần hoàn/not-a-knot | Tái tạo đa thức phù hợp điều kiện biên, kiểm tra đạo hàm hai đầu, tính liên tục và đạo hàm tại nút |
| Spline bậc 1–4 theo đầu trái và hệ biên tuyến tính tổng quát | Đa thức có hệ số biết trước; nghiệm giải từ hai cách độc lập; đơn vị nhỏ, bảng lớn và điều kiện thiếu/phụ thuộc |
| Hồi quy cơ sở, trọng số, log/lũy thừa/mũ, các dạng nghịch đảo | Nghiệm hệ chuẩn nhỏ tính độc lập, trực giao phần dư có trọng số, quan sát lặp, đổi đơn vị, hạng thật/giả và miền log/chia |
| Tối ưu SSE y gốc cho mũ/lũy thừa | So sánh mục tiêu với khớp log, đổi đơn vị y, dữ liệu âm/0, dự đoán vô hướng/mảng |
| Đạo hàm p+1 điểm | Tái tạo đạo hàm đa thức nhiều cấp, biên, triệt tiêu trên bảng hằng và bậc 100; kiểm thử đợt trước vẫn qua |
| Hình thang, trung điểm, Simpson | Bậc chính xác đại số; bảng và biểu thức cùng kết quả; cận lý thuyết/Runge; đảo cận, cận bằng nhau và thang cực trị |
| Newton–Cotes ghép bậc 1–10 | Tích phân chính xác đa thức tới đúng bậc lý thuyết: m khi m lẻ, m+1 khi m chẵn; điều kiện chia hết và bậc Runge |
| Tích phân kép Simpson/hình thang | Trọng số tensor, đa thức hai biến, đảo cận/miền suy biến và đối chiếu lưới bổ sung |
| Euler hiện, RK2/3/4 | Nghiệm không tự trị x′=x+sin(t), tỷ lệ sai số khi h giảm đôi; hệ 2–3 thành phần |
| RK2 tham số và RK3 tham số, bảng RK nhập tay | Các tham số tương đương, điều kiện bậc, hai họ RK3 đặc biệt; từ chối bảng RK hiện không hợp lệ |
| Euler ẩn/hình thang | Nghiệm sai phân giải tích của x′=−100x; từ chối bước vô nghiệm/suy biến; Jacobian và phần dư |
| Adams 2–5 | Kiểm tra mômen độc lập: AB tích phân đúng t^k tới k=s−1, AM tới k=s; nghiệm chính xác của x′=t² với RK4 khởi tạo |
| Bài toán biên Dirichlet/Robin/Neumann | Nghiệm chế tạo, dấu của nguồn và q, quy đổi đơn vị điều kiện biên, p biến thiên, bậc hội tụ và trường hợp không duy nhất |
| Sturm–Liouville | Phổ sai phân giải tích, một nút trong, nhiều cặp riêng, phần dư A·v−λR·v, trực giao có trọng số và r rất lớn/nhỏ |
| Gauss chọn trụ và Thomas | Hệ cần đổi hàng, suy biến, pivot 0, đơn vị khác thang và nghiệm chế tạo cho nhiều kích thước |

## Cách đọc các bảo đảm

- Bậc chính xác đại số, phần dư hệ rời rạc và sai số so với nghiệm liên tục là các đại lượng khác nhau. ODE/BVP/trị riêng cần xét h và h/2; dung sai của bước ẩn chỉ kiểm tra phép giải bước đó.
- `Adams s` hiện có nghĩa AB dùng s giá trị cũ và AM thêm giá trị mới, với RK4 khởi tạo. Bậc toàn cục còn chịu ảnh hưởng của phần khởi tạo; không tự lấy bậc AM làm bậc của toàn quá trình.
- Mỗi họ RK3 đặc biệt có một quy ước được in bằng bảng Butcher. Nếu đề cho các hệ số tự do khác, dùng **RK tùy chỉnh** để nhập đúng bảng đề yêu cầu.
- Cận co nội suy ngược được tính với nhiều chữ số và có cộng phần làm tròn nghiệm sang float. Chưa dùng số học khoảng để chứng nhận mọi phép làm tròn hệ số; cận vẫn chỉ cho đa thức mô hình.
- Hồi quy phi tuyến vẫn có thể đến cực tiểu cục bộ; QR/Decimal không loại bỏ độ nhạy của bài toán. Trọng số rất nhỏ vẫn là quan sát có trọng số, không được tự xóa.
- Các giới hạn tài nguyên lưới tích phân/ODE/BVP còn giữ. Chưa kiểm tra một triệu khoảng tích phân hoặc kích thước tùy ý lớn. Tài liệu và dữ liệu đề riêng còn thiếu như báo cáo chính đã nêu.

Kết quả kiểm thử cuối cùng và lệnh chạy được cập nhật tại `Bao_cao_kiem_dinh.md`. Những bằng chứng trên không phải chứng minh chương trình đúng cho mọi đầu vào.
