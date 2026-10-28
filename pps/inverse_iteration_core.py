"""Số học hữu tỉ cho lặp nội suy ngược và chứng nhận trên cả đoạn.
Tách từ PPS/pplap.py; không phụ thuộc đường dẫn dự án PPS.
"""
import math
from decimal import Decimal, localcontext
from fractions import Fraction as F

def rat(value):
    try:
        return F(str(value))
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"Số không hợp lệ hoặc không hữu hạn: {value}") from exc


def fmt(value, digits=12):
    if value is None:
        return "—"
    value = rat(value)
    with localcontext() as ctx:
        ctx.prec = max(60, digits + 20)
        return f"{Decimal(value.numerator)/Decimal(value.denominator):.{digits}f}"


def sci(value):
    with localcontext() as ctx:
        ctx.prec = 16
        return f"{Decimal(value.numerator)/Decimal(value.denominator):.6E}"


def table(headers, rows):
    rows = [[str(v) for v in row] for row in rows]
    widths = [max(len(h), *(len(row[j]) for row in rows)) for j, h in enumerate(headers)]
    return "\n".join(" | ".join(v.rjust(w) for v, w in zip(row, widths))
                     for row in [headers, *rows])


def validate(x, y):
    x, y = list(map(rat, x)), list(map(rat, y))
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("Cần ít nhất hai cặp (x_i,y_i), hai cột cùng độ dài.")
    if any(b <= a for a, b in zip(x, x[1:])):
        raise ValueError("Các x_i phải tăng nghiêm ngặt, không trùng nhau.")
    return x, y


def difference_table(y, x=None):
    """Cột r là sai phân bậc r; nếu có x thì là tỷ sai phân."""
    cols = [list(y)]
    for r in range(1, len(y)):
        prev = cols[-1]
        cols.append([(prev[i+1]-prev[i]) / (x[i+r]-x[i] if x else 1)
                     for i in range(len(y)-r)])
    return cols


def mul_linear(poly, root):
    out = [F(0)] * (len(poly)+1)
    for i, c in enumerate(poly):
        out[i] -= root*c
        out[i+1] += c
    return out


def evaluate(poly, u):
    value = F(0)
    for c in reversed(poly):
        value = value*u+c
    return value


def derivative(poly):
    return [i*poly[i] for i in range(1, len(poly))] or [F(0)]


def bernstein_range(poly, a, b):
    """Bao chắc chắn: min b_k <= p(u) <= max b_k, u thuộc [a,b]."""
    degree = len(poly)-1
    # p(a+(b-a)v)=sum c_j v^j; chuyển sang cơ sở Bernstein bậc degree.
    c = [sum((poly[i]*math.comb(i, j)*a**(i-j)*(b-a)**j
              for i in range(j, degree+1)), F(0)) for j in range(degree+1)]
    bernstein = [sum((c[j]*F(math.comb(k, j), math.comb(degree, j))
                      for j in range(k+1)), F(0)) for k in range(degree+1)]
    return min(bernstein), max(bernstein)


def certify(phi, interval, q=None, max_depth=8):
    """Bao bằng phân đoạn đều; không chứng nhận được không có nghĩa là phân kỳ."""
    a, b = interval
    requested = None if q is None else rat(q)
    if requested is not None and not 0 <= requested < 1:
        raise ValueError("Hệ số co q phải thỏa 0 <= q < 1.")
    dp = derivative(phi)
    for depth in range(max_depth+1):
        parts = 2**depth
        ranges, slopes = [], []
        for j in range(parts):
            left, right = a+(b-a)*F(j, parts), a+(b-a)*F(j+1, parts)
            ranges.append(bernstein_range(phi, left, right))
            slopes.append(bernstein_range(dp, left, right))
        low, high = min(v[0] for v in ranges), max(v[1] for v in ranges)
        bound = max(max(abs(l), abs(r)) for l, r in slopes)
        used = bound if requested is None else requested
        if a <= low and high <= b and bound <= used < 1:
            return {"image": (low, high), "bound": bound, "q": used, "parts": parts}
    raise ValueError("Chưa chứng nhận được phi(I) ⊆ I và |phi'| <= q < 1 "
                     "bằng bao Bernstein (tối đa 256 đoạn). Hãy chọn lại mốc, đoạn "
                     "hoặc phép lặp; không kết luận phép lặp phân kỳ.")


def build_model(x, y, target, method):
    if method not in ("tien", "lui", "khong-deu"):
        raise ValueError("Phương pháp phải là tien, lui hoặc khong-deu.")
    h = x[1]-x[0]
    if method != "khong-deu" and any(x[i] != x[0]+i*h for i in range(len(x))):
        raise ValueError("Mốc không cách đều chính xác; dùng --method khong-deu "
                         "hoặc chuẩn hóa dữ liệu đầu vào theo bước lưới đã biết.")
    cols = difference_table(y, x if method == "khong-deu" else None)
    coefficients = [col[-1] if method == "lui" else col[0] for col in cols]
    base = x[-1] if method == "lui" else x[0]
    scale = F(1) if method == "khong-deu" else h
    if coefficients[1] == 0:
        raise ValueError("Sai phân/tỷ sai phân bậc một bằng 0: không dùng được phép tách này.")
    # u=(x-base)/scale. Không đều cũng dùng u=x-x0 (không chuẩn hóa bởi h).
    H, basis = [F(0)]*len(x), [F(1)]
    terms = []
    for r in range(1, len(x)):
        root = x[r-1]-base if method == "khong-deu" else F(-(r-1) if method == "lui" else r-1)
        basis = mul_linear(basis, root)
        if method != "khong-deu":
            basis = [c/r for c in basis]
        terms.append((coefficients[r], basis[:]))
        if r >= 2:
            for i, c in enumerate(basis):
                H[i] += coefficients[r]*c
    phi = [-c/coefficients[1] for c in H]
    phi[0] += (target-coefficients[0])/coefficients[1]
    P = H[:]
    P[0] += coefficients[0]
    P[1] += coefficients[1]
    return dict(base=base, scale=scale, nodes=x, cols=cols, coefficients=coefficients,
                phi=phi, P=P, terms=terms)


def solve(x, y, target, interval, method="tien", epsilon="1e-8", q=None, max_iter=200):
    """interval nhập theo x; không tự chọn nhánh hay giảm số mốc người dùng đã chọn."""
    x, y = validate(x, y)
    target, epsilon = rat(target), rat(epsilon)
    if epsilon <= 0 or isinstance(max_iter, bool) or not isinstance(max_iter, int) or max_iter < 1:
        raise ValueError("epsilon phải dương; max_iter phải là số nguyên dương.")
    a, b = map(rat, interval)
    if not x[0] <= a < b <= x[-1]:
        raise ValueError("Đoạn tìm nghiệm phải nằm trong miền các mốc đã chọn và a < b.")
    if q is not None and not 0 <= rat(q) < 1:
        raise ValueError("Hệ số co q phải thỏa 0 <= q < 1.")
    matches = [xi for xi, yi in zip(x, y) if a <= xi <= b and yi == target]
    if len(matches) > 1:
        raise ValueError("Y trùng nhiều mốc trên đoạn; cần chọn nhánh cụ thể hơn.")
    if matches:
        return dict(exact_node=matches[0], x=x, y=y, target=target, interval=(a,b))
    model = build_model(x, y, target, method)
    base, scale = model['base'], model['scale']
    I = ((a-base)/scale, (b-base)/scale)
    proof = certify(model['phi'], I, q)
    u = (target-model['coefficients'][0])/model['coefficients'][1]
    fallback = not I[0] <= u <= I[1]
    if fallback:
        u = (I[0]+I[1])/2
    history = [(0, u, base+scale*u, None, F(0))]
    # Làm tròn hữu tỉ đến 40 chữ số, đo chính xác sai số từng phép cập nhật.
    grid = 10**40
    for k in range(1, max_iter+1):
        exact_update = evaluate(model['phi'], u)
        new = F((2*exact_update.numerator*grid+exact_update.denominator)//
                (2*exact_update.denominator), grid)
        rho = abs(new-exact_update)
        if not I[0] <= new <= I[1]:
            raise ValueError("Giá trị sau làm tròn ra ngoài đoạn lặp; chưa nhận kết quả.")
        error = scale*(proof['q']*abs(new-u)+rho)/(1-proof['q'])
        history.append((k, new, base+scale*new, error, rho))
        if error <= epsilon:
            break
        u = new
    success = history[-1][3] <= epsilon
    return dict(x=x, y=y, target=target, interval=(a,b), method=method, model=model,
                proof=proof, I=I, history=history, epsilon=epsilon, success=success,
                fallback=fallback, residual=abs(evaluate(model['P'], history[-1][1])-target))


def polynomial_text(poly):
    return " + ".join(f"({c})*u^{i}" for i, c in enumerate(poly) if c) or "0"


def newton_text(model, method):
    terms = [fmt(model['coefficients'][0])]
    for r in range(1, len(model['coefficients'])):
        if method == 'khong-deu':
            factors = ['u' if j == 0 else f"(u-({fmt(model['nodes'][j]-model['base'])}))"
                       for j in range(r)]
            divisor = ''
        else:
            factors = ['u']+[f"(u{'+' if method == 'lui' else '-'}{j})" for j in range(1,r)]
            divisor = f'/{r}!' if r > 1 else ''
        terms.append(f"({fmt(model['coefficients'][r])})*"+'*'.join(factors)+divisor)
    return ' + '.join(terms)


