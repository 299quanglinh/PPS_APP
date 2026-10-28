"""Polynomial, central and inverse interpolation; splines and regression."""
import math
from decimal import Decimal, localcontext
import numpy as np
import sympy as sp
from numpy.polynomial import Polynomial as Poly
from scipy.interpolate import CubicSpline, PPoly
from .common import Report, points, spacing, difference, fmt, polynomial_text, function
from .polynomial_model import NewtonPolynomial, nonnegative_integer, decimal_cell
from .linear_algebra import fit_basis, sum_products, exp_parameter

METHODS = ['Lagrange', 'Newton bất kỳ', 'Newton tiến', 'Newton lùi',
           'Gauss I', 'Gauss II', 'Stirling', 'Bessel']


def _lagrange_exact_decimal_polynomial(x, y):
    """Return the expanded Lagrange polynomial from the displayed decimals.

    Decimal input such as ``1.2`` is intended as the rational number 6/5, not
    the nearby IEEE-754 binary float.  Expanding the basis in rational
    arithmetic prevents the small coefficient drift that is otherwise visible
    when a user compares a hand calculation with the displayed polynomial.
    """
    if len(x) > 20:
        return None
    symbol = sp.Symbol('x')
    try:
        nodes = [(sp.Rational(str(float(xi))), sp.Rational(str(float(yi))))
                 for xi, yi in zip(x, y)]
        coefficients = sp.Poly(sp.interpolate(nodes, symbol), symbol).all_coeffs()
        return Poly([float(value) for value in reversed(coefficients)])
    except (TypeError, ValueError, sp.PolynomialError):
        return None


def interpolate(x, y, method='Lagrange', target=None, derivative=0, count=None):
    """Nội suy toàn bảng hoặc đúng count mốc khi người dùng yêu cầu.

    Bảng lớn dùng dạng tích/Newton lồng thay vì ép khai triển theo x.
    Không tự bỏ mốc hay tự thay bậc đa thức của đề.
    """
    x, y, dropped = points(x, y)
    derivative = nonnegative_integer(derivative, 'Cấp đạo hàm')
    if target is not None and not np.isfinite(target):
        raise ValueError('Điểm tính phải hữu hạn.')
    if method not in METHODS+['Newton bất kỳ lùi']:
        raise ValueError('Phương pháp không hợp lệ.')
    original_count = len(x)
    if count is not None:
        count = nonnegative_integer(count, 'Số mốc')
        if not 1 <= count <= len(x):
            raise ValueError('Số mốc sử dụng phải từ 1 đến số mốc trong bảng.')
        if count < len(x):
            if target is None:
                raise ValueError('Nhập điểm cần tính để chọn nhóm mốc; để trống số mốc nếu muốn dùng toàn bảng.')
            x, y = select_points(x, y, target, count).value
    n = len(x)
    central = method in ('Gauss I', 'Gauss II', 'Stirling', 'Bessel')
    if method == 'Stirling' and n % 2 == 0:
        raise ValueError(f'Stirling cần số mốc lẻ; bảng hiện có {n} mốc. Nhập số mốc dùng (ví dụ 9) để chọn nhóm quanh điểm tính.')
    if method == 'Bessel' and n % 2:
        raise ValueError(f'Bessel cần số mốc chẵn (ít nhất 2); đang dùng {n} mốc. Nhập số mốc dùng (ví dụ 8).')
    if method not in ('Lagrange', 'Newton bất kỳ', 'Newton bất kỳ lùi'):
        spacing(x)
    model = NewtonPolynomial(x, y)
    if n > 32:
        result = _interpolate_compact(x, y, method, target, derivative, model)
    else:
        try:
            with np.errstate(over='raise', invalid='raise', divide='raise'):
                result = _interpolate_expanded(x, y, method, target, derivative, model)
        except (FloatingPointError, OverflowError, ZeroDivisionError):
            result = _interpolate_compact(x, y, method, target, derivative, model)
    result.presentation['selection'] = {'input_count': original_count, 'used_count': n}
    if n != original_count:
        result.lines.insert(0, f'Đã chọn {n}/{original_count} mốc quanh x = {fmt(target)}; đa thức bậc không quá {n-1} chỉ đi qua nhóm mốc này.')
        result.table('Mốc được sử dụng', ['x', 'y'], zip(x, y))
    return result


def _interpolate_compact(x, y, method, target, derivative, model):
    """Bảng lớn giữ hệ số dạng tích; tránh tràn số do khai triển/giai thừa."""
    n = len(x)
    r = Report(method, [
        f'Dùng đủ {n} mốc; đa thức bậc không quá {n-1}.',
        'Giữ dạng tích để tránh tràn số khi khai triển theo x. Đánh giá bằng Horner Newton tương đương trên cùng toàn bộ mốc.',
        f'Số chữ số tính trung gian: {model.precision}; bảng hệ số hiển thị 17 chữ số có nghĩa.',
        'Bậc cao có thể khuếch đại sai số dữ liệu và dao động mạnh giữa các mốc; tăng số chữ số không loại bỏ hiện tượng này.'
    ], value=model, scatter=(x, y))
    tex = r'P(x)=c_0+(z-z_0)\left[c_1+(z-z_1)\left[\cdots\right]\right]'
    r.presentation = {'formulas': [('Dạng tích tương đương', [tex])], 'compact_polynomial': True}
    formulas = {
        'Lagrange': r'P(x)=\sum_{i=0}^{n}D_i\prod_{j\ne i}(x-x_j),\quad D_i=\frac{y_i}{\prod_{j\ne i}(x_i-x_j)}',
        'Newton tiến': r'P(x_0+th)=\sum_{k=0}^{n}\frac{\Delta^k y_0}{k!}\prod_{j=0}^{k-1}(t-j)',
        'Newton lùi': r'P(x_n+th)=\sum_{k=0}^{n}\frac{\Delta^k y_{n-k}}{k!}\prod_{j=0}^{k-1}(t+j)',
        'Stirling': r'P(x_0+th)=y_0+t\frac{\Delta y_{-1}+\Delta y_0}{2}+\frac{t^2}{2!}\Delta^2y_{-1}+\cdots',
        'Bessel': r'P(x_0+th)=\frac{y_0+y_1}{2}+(t-\frac{1}{2})\Delta y_0+\frac{t(t-1)}{2!}\frac{\Delta^2y_{-1}+\Delta^2y_0}{2}+\cdots',
        'Gauss I': r'P(x_0+th)=y_0+t\Delta y_0+\frac{t(t-1)}{2!}\Delta^2y_{-1}+\cdots',
        'Gauss II': r'P(x_0+th)=y_0+t\Delta y_{-1}+\frac{t(t+1)}{2!}\Delta^2y_{-1}+\cdots'
    }
    if method in formulas:
        r.presentation['formulas'].insert(0, ('Công thức '+method+' (n = số mốc − 1)', [formulas[method]]))
    r.lines.append(f'z = (x - ({model.origin})) / ({model.scale}); các c_k ở bảng Newton theo đúng thứ tự z_i.')
    r.table('Hệ số Newton lồng tương đương — dùng toàn bộ mốc', ['k', 'x_i', 'z_i', 'c_k'],
            ([k, format(node, '.17g'), format(z, '.17g'), format(c, '.17g')]
             for k, (node, z, c) in enumerate(zip(model.raw_nodes, model.nodes, model.coefficients))))
    if method == 'Lagrange':
        rows = []
        with localcontext() as ctx:
            ctx.prec = model.precision
            nodes = [Decimal(str(float(v))) for v in x]
            for i, node in enumerate(nodes):
                denominator = Decimal(1)
                for j, other in enumerate(nodes):
                    if i != j: denominator *= node-other
                rows.append([i, x[i], format(denominator, '.17g'), format(Decimal(str(float(y[i])))/denominator, '.17g')])
        r.table('Hệ số Lagrange dạng tích', ['i', 'x_i', 'tích (x_i − x_j)', 'D_i'], rows)
    if method != 'Lagrange':
        normalized = not method.startswith('Newton bất kỳ')
        center = n//2 if method=='Gauss II' else (n-1)//2
        h = spacing(x) if normalized else 1.
        origin = x[center] if method in ('Gauss I', 'Gauss II', 'Stirling', 'Bessel') else x[-1] if method.endswith('lùi') else x[0]
        if normalized:
            r.lines.append(f'Công thức {method}: t = (x - ({fmt(origin)})) / ({fmt(h)}). Bảng sau giữ Δ^k/k! ngoài dạng tích.')
        rows = []
        with localcontext() as ctx:
            ctx.prec = model.precision
            column = [Decimal(str(float(v))) for v in y]
            dx = [Decimal(str(float(v))) for v in x]
            factorial = Decimal(1)
            for k in range(n):
                if k: factorial *= k
                if method == 'Stirling':
                    j = center-k//2
                    value = (column[j-1]+column[j])/2 if k%2 else column[j]
                    source = f'trung bình hàng {j-1} và {j}' if k%2 else f'hàng {j}'
                elif method == 'Bessel':
                    j = center-k//2
                    value = column[j] if k%2 else (column[j]+column[j+1])/2
                    source = f'hàng {j}' if k%2 else f'trung bình hàng {j} và {j+1}'
                else:
                    j = center-k//2 if method=='Gauss I' else center-(k+1)//2 if method=='Gauss II' else n-k-1 if method.endswith('lùi') else 0
                    value = column[j]
                    source = f'hàng {j}'
                rows.append([k, source, format(value, '.17g'), format(value/factorial if normalized else value, '.17g')])
                column = [column[j+1]-column[j] for j in range(len(column)-1)]
                if not normalized:
                    column = [v/(dx[j+k+1]-dx[j]) for j, v in enumerate(column)]
        r.table('Hệ số công thức '+method, ['k', 'vị trí (hàng từ 0)', 'sai phân/tỷ sai phân', 'hệ số sau chia k! (mốc đều)'], rows)
    r.table('Toàn bộ dữ liệu nội suy', ['x', 'y'], zip(x, y))
    if target is not None:
        val = model.deriv(derivative)(target)
        r.lines.append(f'P^({derivative})({fmt(target)}) = {fmt(val)}')
        r.presentation['formulas'].append(('Thay giá trị', [rf'P^{{({derivative})}}({fmt(target)})={fmt(val)}']))
        if target < x[0] or target > x[-1]: r.lines.append('Điểm tính ngoài tập mốc: đang ngoại suy.')
    # Đồ thị là phần phụ; tràn số ở nơi khác không làm mất kết quả tại điểm yêu cầu.
    try:
        grid = np.linspace(x[0], x[-1], 41)
        r.curves = [(grid, model(grid), method)]
    except ValueError:
        r.lines.append('Không vẽ đồ thị: có giá trị ngoài miền float trên lưới vẽ; kết quả tại điểm yêu cầu vẫn được tính riêng.')
    return r


def _interpolate_expanded(x, y, method, target, derivative, model):
    n = len(x)
    if method not in METHODS+['Newton bất kỳ lùi']: raise ValueError('Phương pháp không hợp lệ.')
    central = method in ('Gauss I', 'Gauss II', 'Stirling', 'Bessel')
    if method == 'Bessel' and n % 2:
        raise ValueError('Bessel cần số mốc chẵn, ít nhất 2 mốc.')
    normalized = method not in ('Lagrange', 'Newton bất kỳ', 'Newton bất kỳ lùi')
    h = spacing(x) if normalized else 1.0
    center = (n//2 if method=='Gauss II' else (n-1)//2) if central else (n-1 if method == 'Newton lùi' else 0)
    origin = x[center] if normalized else 0.0
    a = difference(x, y, method.startswith('Newton bất kỳ'))
    q, product = Poly([0.]), Poly([1.])
    steps, products = [], []
    lagrange_data = None
    if method == 'Lagrange':
        for i in range(n):
            basis = Poly.fromroots(np.delete(x, i)) if n>1 else Poly([1.])
            c = y[i]/basis(x[i])
            q += c*basis
            steps.append([i, c, polynomial_text(basis)])
            products.append(basis.coef)
        exact_q = _lagrange_exact_decimal_polynomial(x, y)
        if exact_q is not None:
            q = exact_q
        matrix = [[1.0 if i == j else x[i] - x[j] for j in range(n)]
                  for i in range(n)]
        denominators = np.asarray([np.prod(row) for row in matrix])
        weights = y / denominators
        prefix_products = [Poly([1.])]
        for node in x:
            prefix_products.append(prefix_products[-1] * Poly([-node, 1.]))
        lagrange_data = matrix, weights, prefix_products
    elif method in ('Newton bất kỳ', 'Newton bất kỳ lùi', 'Newton tiến', 'Newton lùi', 'Gauss I', 'Gauss II'):
        if method == 'Gauss I':
            order = [0] + [v for r in range(1, n) for v in (r, -r) if 0<=center+v<n]
        elif method == 'Gauss II':
            order = [0] + [v for r in range(1, n) for v in (-r, r) if 0<=center+v<n]
        elif method == 'Newton bất kỳ lùi': order=list(x[::-1])
        elif method == 'Newton lùi': order = list(-np.arange(n))
        else: order = list(np.arange(n)) if normalized else list(x)
        for k in range(n):
            if k:
                product *= Poly([-order[k-1], 1])
            if method == 'Gauss I': j = center-k//2
            elif method == 'Gauss II': j = center-(k+1)//2
            elif method in ('Newton lùi','Newton bất kỳ lùi'): j = n-k-1
            else: j = 0
            q += (a[j,k]/math.factorial(k) if normalized else a[j,k])*product
            steps.append([k, a[j,k], polynomial_text(product, 't' if normalized else 'x')])
            products.append(product.coef)
    elif method == 'Stirling':
        q = Poly([y[center]])
        products.append(Poly([1.]).coef)
        for i in range(center):
            odd = (a[center-i-1,2*i+1]+a[center-i,2*i+1])/(2*math.factorial(2*i+1))
            even = a[center-i-1,2*i+2]/math.factorial(2*i+2)
            products.extend([(Poly([0,1])*product).coef,
                             (Poly([0,0,1])*product).coef])
            q += Poly([0, odd, even])*product
            steps.append([i, odd, even, polynomial_text(product, 't')])
            product *= Poly([-(i+1)**2,0,1])
    elif method == 'Bessel':
        for i in range(center+1):
            even = (a[center-i,2*i]+a[center-i+1,2*i])/(2*math.factorial(2*i))
            odd = a[center-i,2*i+1]/math.factorial(2*i+1)
            products.extend([product.coef,
                             (Poly([-.5,1])*product).coef])
            q += Poly([even-odd/2,odd])*product
            steps.append([i, odd, even, polynomial_text(product, 't')])
            product *= Poly([i,1])*Poly([-i-1,1])
    p = q(Poly([-origin/h,1/h])) if normalized else q
    if not np.all(np.isfinite(p.coef)):
        raise FloatingPointError('Hệ số khai triển không hữu hạn.')
    result = Report(method, value=model, scatter=(x,y))
    if normalized:
        origin_label = 'x_n' if method == 'Newton lùi' else 'x_0'
        result.lines += [f'{origin_label} = {fmt(origin)}; h = {fmt(h)}; t = (x - ({fmt(origin)}))/{fmt(h)}',
                         f'P({origin_label}+th) = '+polynomial_text(q, 't')]
    if method not in ('Newton tiến', 'Newton lùi'):
        result.lines.append('P(x) = '+polynomial_text(p))
    result.lines.append('Hệ số hiển thị làm tròn; phép tính dùng hệ số đầy đủ.')
    err = np.max(np.abs(q((x-origin)/h) - y)) if normalized else np.max(np.abs(p(x)-y))
    result.lines.append('Sai lệch lớn nhất tại các mốc: '+fmt(err, 12))
    if err > 1e-6 * max(1, np.max(np.abs(y))):
        result.lines.append('CẢNH BÁO: đa thức bậc cao/khoảng mốc gây mất chính xác; hãy giảm số mốc hoặc dùng spline.')
    if target is not None:
        if derivative < 0: raise ValueError('Cấp đạo hàm không âm.')
        val = model.deriv(derivative)(target)
        result.lines.append(f'P^({derivative})({fmt(target)}) = {fmt(val)}')
        if target < x[0] or target > x[-1]: result.lines.append('Điểm tính ngoài tập mốc: đang ngoại suy.')
    if method != 'Lagrange':
        result.table('Bảng tỷ sai phân' if method.startswith('Newton bất kỳ') else 'Bảng sai phân tiến',
                     ['x','y']+[f'bậc {k}' for k in range(1,n)],
                     [[x[j]]+[None if np.isnan(v) else v for v in a[j]] for j in range(n)])
    if method == 'Lagrange':
        matrix, weights, prefix_products = lagrange_data
        width = n + 1
        def pad(coefficients):
            values = coefficients.coef if isinstance(coefficients, Poly) else coefficients
            return list(reversed(values)) + [None] * (width - len(values))
        def pad_prefix_product(coefficients):
            # Bảng tích của ứng dụng gốc ghi hệ số theo bậc tăng dần:
            # cột 0 là hệ số tự do, cột cuối là bậc cao nhất.  Vì vậy
            # w_0 = 1 phải nằm ở cột cuối, tạo tam giác chụm về bên trái.
            values = coefficients.coef if isinstance(coefficients, Poly) else coefficients
            return [None] * (width - len(values)) + list(values)
        rows = [[*row, None, 'Ma trận (x_i - x_j)' if i == 0 else ''] for i, row in enumerate(matrix)]
        rows += [[None] * width + ['']]
        rows += [[*weights, None, 'Hệ số D_i = y_i / ∏(x_i - x_j)']]
        rows += [[None] * width + ['']]
        rows += [[*pad_prefix_product(product), 'Bảng tích' if i == 0 else 'w_(n+1)' if i == n else '']
                 for i, product in enumerate(prefix_products)]
        rows += [[None] * width + ['']]
        rows += [[*pad(basis), 'Bảng thương W_i = w_(n+1)/(x-x_i)' if i == 0 else '']
                 for i, basis in enumerate(products)]
        rows += [[None] * width + ['']]
        rows += [[*pad(p), 'Hệ số đa thức nội suy P(x)']]
        result.table('Các bước lập đa thức Lagrange',
                     [f'Cột {k}' for k in range(width)] + ['Ghi chú'], rows)
    elif products:
        if method.startswith('Newton bất kỳ'):
            table_name = 'Bảng tích: đa thức cơ sở Newton — hệ số theo bậc giảm dần'
        else:
            table_name = 'Bảng tích: chưa chia giai thừa — hệ số theo bậc giảm dần'
        result.table(table_name, ['k']+[f'bậc {k}' for k in range(n-1,-1,-1)],
                     [[i]+list(reversed(list(c)+[0]*(n-len(c)))) for i,c in enumerate(products)])
    grid = np.linspace(x[0],x[-1],400)
    result.curves = [(grid, model(grid), method)]
    from .math_view import interpolation_presentation
    interpolation_presentation(result,method,x,y,a,q,p,origin,h,center,normalized,target,derivative)
    if normalized:
        from .math_view import product_formula
        product_formula(result,method,products,center)
    return result


def inverse(x, y, target, backward=False, tolerance=1e-10, max_steps=1000):
    x,y,_ = points(x,y,2)
    if not np.isfinite(target) or not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('Mức y phải hữu hạn; dung sai phải dương và hữu hạn.')
    h = spacing(x)
    if not (np.all(np.diff(y)>0) or np.all(np.diff(y)<0)):
        raise ValueError('Chọn một đoạn dữ liệu đơn điệu nghiêm ngặt cho nội suy ngược.')
    if not min(y) <= target <= max(y): raise ValueError('Giá trị y cần tìm nằm ngoài đoạn đã chọn.')
    if tolerance <= 0: raise ValueError('Dung sai phải dương.')
    max_steps=nonnegative_integer(max_steps,'Số bước lặp tối đa')
    if not max_steps:raise ValueError('Số bước lặp tối đa phải dương.')
    at_node=np.flatnonzero(y==target)
    if len(at_node):return Report('Nội suy ngược — mức tại mốc',[f'Mức bằng tung độ tại mốc x = {fmt(x[at_node[0]])}.'],value=float(x[at_node[0]]))
    from .inverse_model import InverseMap, evaluate
    model=InverseMap(x,y,target,backward,precision=max(50,-Decimal(str(float(tolerance))).adjusted()+40))
    rows=[]
    with localcontext() as ctx:
        ctx.prec=model.precision
        dtol=Decimal(str(float(tolerance)))
        lo=(Decimal(str(float(x[0])))-model.origin)/model.h
        hi=(Decimal(str(float(x[-1])))-model.origin)/model.h
        t=model.initial()
        if not lo<=t<=hi:t=(lo+hi)/2
        for k in range(1,max_steps+1):
            next_t=model.step(t)
            if not lo<=next_t<=hi:break
            total=evaluate(model.polynomial,t)-model.y0-model.d1*t
            rows.append([k,decimal_cell(t),decimal_cell(total),decimal_cell(next_t),decimal_cell(abs(next_t-t))])
            residual=model.residual(next_t)
            if abs(next_t-t)*model.h<=dtol and residual<=dtol*max(Decimal(1),abs(model.target)):
                exact=model.location(next_t);result=float(exact)
                if abs(Decimal.from_float(result)-exact)>dtol:
                    raise ValueError('Dung sai vị trí nhỏ hơn độ phân giải float; đổi gốc/đơn vị hoặc tăng dung sai.')
                report=Report('Nội suy ngược — lặp Newton '+('lùi' if backward else 'tiến'),
                              [f'x ≈ {fmt(result)}; sai lệch P(x)-y = {fmt(decimal_cell(residual),12)}',
                               'Dừng thực hành theo sai khác bước và phần dư; chưa có cận co/hội tụ được chứng minh.',
                               'Nghiệm là của đa thức nội suy; sai số so với hàm thật cần cận mô hình.'],value=result).table(
                                  'Lịch sử lặp',['k','t trước','tổng bậc cao','t sau','độ lệch'],rows)
                report.presentation={'formulas':[('Nghiệm nội suy ngược',[rf'P(x)={fmt(target)},\qquad x\approx {fmt(result)}'])]}
                return report
            if next_t==t:break
            t=next_t
    raise ValueError('Phép lặp không hội tụ. Chọn ít mốc hơn gần giá trị y cần tìm hoặc đổi hướng tiến/lùi.')


def select_points(x,y,target,count=5,reverse=False):
    x,y,_ = points(x,y,2 if reverse else 1)
    if not np.isfinite(target): raise ValueError('Điểm/mức cần chọn phải hữu hạn.')
    if not reverse: count = nonnegative_integer(count, 'Số mốc')
    r = Report('Chọn mốc nội suy')
    if reverse:
        seen = set()
        for k in range(len(x)-1):
            if min(y[k:k+2]) <= target <= max(y[k:k+2]) and y[k]!=y[k+1]:
                sign = np.sign(y[k+1]-y[k]); left,right=k,k+1
                while left>0 and np.sign(y[left]-y[left-1])==sign: left-=1
                while right<len(x)-1 and np.sign(y[right+1]-y[right])==sign: right+=1
                if (left,right) in seen: continue
                seen.add((left,right))
                r.table(f'Đoạn đơn điệu [{fmt(x[left])}, {fmt(x[right])}]',['x','y'],zip(x[left:right+1],y[left:right+1]))
        if not seen: raise ValueError('Không tìm được đoạn đơn điệu chứa y cần tìm.')
    else:
        if not 1<=count<=len(x): raise ValueError('Số mốc phải nằm trong độ dài bảng.')
        try: spacing(x); equal=True
        except ValueError: equal=False
        if equal:
            # Gần tâm nhóm nhất; khi hòa chọn cửa sổ bên trái.
            start = min(range(len(x)-count+1),key=lambda j: (abs((x[j]+x[j+count-1])/2-target),j))
            chosen=np.arange(start,start+count)
        else:
            chosen=np.sort(np.argsort(abs(x-target),kind='stable')[:count])
        xx,yy=x[chosen],y[chosen]
        r.table('Mốc đã chọn',['x','y'],zip(xx,yy))
        r.lines += ['Dãy x: '+' '.join(fmt(v) for v in xx),'Dãy y: '+' '.join(fmt(v) for v in yy)]
        if not equal:r.lines.append('Mốc không đều: dùng Lagrange hoặc Newton tỷ sai phân.')
        if target<x[0] or target>x[-1]:r.lines.append('Điểm cần tính ngoài bảng: ngoại suy.')
        r.value=(xx,yy)
    return r


def chebyshev(degree,a,b):
    degree=nonnegative_integer(degree, 'Bậc Chebyshev')
    if not np.all(np.isfinite([a,b])) or not a<b: raise ValueError('Cần a < b và hai đầu hữu hạn.')
    x=np.sort((a/2+b/2)+(b/2-a/2)*np.cos((2*np.arange(degree+1)+1)*np.pi/(2*(degree+1))))
    if np.any(np.diff(x)==0):
        raise ValueError('Mốc Chebyshev bị trùng trong độ chính xác float; giảm bậc hoặc đổi gốc/đơn vị trước khi nhập dữ liệu.')
    report=Report('Mốc Chebyshev',value=x).table('Mốc tối ưu',['i','x'],enumerate(x))
    report.presentation={'formulas':[('Công thức mốc Chebyshev',[r'x_i=\frac{a+b}{2}+\frac{b-a}{2}\cos\frac{(2i+1)\pi}{2(n+1)}'])]}
    return report


def horner(coefficients,target,order=0):
    from .math_view import poly_latex,split_equation
    order=nonnegative_integer(order, 'Cấp đạo hàm')
    coefficients=np.asarray(coefficients,dtype=float)
    if coefficients.ndim!=1 or not len(coefficients) or not np.all(np.isfinite(coefficients)) or not np.isfinite(target):
        raise ValueError('Cần hệ số và điểm tính hữu hạn.')
    p=Poly(coefficients[::-1])
    r=Report('Horner', ['Hệ số nhập theo bậc giảm dần.','P(x) = '+polynomial_text(p)])
    r.presentation={'formulas':[('Đa thức',split_equation('P(x)',poly_latex(p)))],'omit_prefixes':('P(x) =','P^(')}
    # Decimal tránh tràn giai thừa trung gian và triệt tiêu trong chia lặp.
    working=[Decimal(str(float(v))) for v in coefficients];factorial=Decimal(1)
    with localcontext() as context:
        context.prec=max(50,2*len(coefficients)+30)
        for k in range(min(order,len(coefficients)-1)+1):
            row=[];remainder=Decimal(0)
            for v in working:
                remainder=remainder*Decimal(str(float(target)))+v
                row.append(remainder)
            derivative_value=float(factorial*remainder)
            if not np.isfinite(derivative_value): raise ValueError('Giá trị Horner/đạo hàm vượt miền float; đổi đơn vị hoặc dùng ít bậc hơn.')
            if k==0:r.value=derivative_value
            r.table(f'Horner lặp lần {k+1} — số dư nhân {k}!',
                    ['hệ số','giá trị tích lũy'],((decimal_cell(a),decimal_cell(b)) for a,b in zip(working,row)))
            r.lines.append(f'P^({k})({fmt(target)}) = {fmt(derivative_value)}')
            label='P' if k==0 else rf'P^{{({k})}}'
            r.presentation['formulas'].append((f'Giá trị — cấp {k}',[rf'{label}({fmt(target)})={fmt(derivative_value)}']))
            working=row[:-1] if len(row)>1 else [Decimal(0)]
            factorial*=k+1
    if order>=len(coefficients):
        r.lines += [f'Mọi đạo hàm từ cấp {len(coefficients)} trở lên đều bằng 0.',f'P^({order})({fmt(target)}) = 0']
    return r


def spline(x,y,kind='Cubic tự nhiên',left=0.,right=0.,target=None):
    from .math_view import poly_latex,split_equation
    x,y,dropped=points(x,y,2); h=np.diff(x); slopes=np.diff(y)/h
    if not np.all(np.isfinite([left,right])):
        raise ValueError('Giá trị điều kiện biên phải hữu hạn.')
    if target is not None and not x[0]<=target<=x[-1]:
        raise ValueError('Điểm tính nằm ngoài miền spline; không tự ngoại suy.')
    valid=('Tuyến tính','Bậc hai — biên trái','Bậc hai — biên phải',
           'Cubic tự nhiên','Cubic kẹp','Cubic đạo hàm bậc hai',
           'Cubic tuần hoàn','Cubic not-a-knot')
    if kind not in valid: raise ValueError('Loại spline không được hỗ trợ.')
    if kind=='Tuyến tính': pp=PPoly(np.array([slopes,y[:-1]]),x,extrapolate=False)
    elif kind.startswith('Bậc hai'):
        d=np.zeros(len(x))
        if kind.endswith('phải'):
            d[-1]=right
            for i in range(len(h)-1,-1,-1): d[i]=2*slopes[i]-d[i+1]
        else:
            d[0]=left
            for i in range(len(h)): d[i+1]=2*slopes[i]-d[i]
        pp=PPoly(np.array([(slopes-d[:-1])/h,d[:-1],y[:-1]]),x,extrapolate=False)
    else:
        if kind=='Cubic tuần hoàn' and not np.isclose(y[0],y[-1],rtol=1e-15,atol=1e-15):
            raise ValueError('Spline tuần hoàn cần y đầu bằng y cuối.')
        if kind=='Cubic not-a-knot' and len(x)<4:
            raise ValueError('Spline not-a-knot cần ít nhất 4 mốc.')
        bc={'Cubic tự nhiên':'natural','Cubic kẹp':((1,left),(1,right)),
            'Cubic đạo hàm bậc hai':((2,left),(2,right)),
            'Cubic tuần hoàn':'periodic','Cubic not-a-knot':'not-a-knot'}[kind]
        pp=CubicSpline(x,y,bc_type=bc,extrapolate=False)
    if not np.all(np.isfinite(pp.c)):raise ValueError('Hệ số spline không hữu hạn; đổi đơn vị hoặc kiểm tra điều kiện biên.')
    r=Report('Spline — '+kind, ['Mỗi đoạn dùng biến (x - mốc trái).'],scatter=(x,y),value=pp)
    r.presentation={'formulas':[],'omit_prefixes':('[','S(')}
    if len(x)>33:
        r.lines.append(f'Dùng đủ {len(x)} mốc và {len(x)-1} đoạn; trình bày hệ số từng đoạn để tránh dựng hệ đặc chỉ phục vụ hiển thị.')
        degree=pp.c.shape[0]-1
        r.table('Hệ số spline theo u = x − mốc trái',['đoạn','x trái','x phải']+[f'u^{k}' for k in range(degree+1)],
                ([i,x[i],x[i+1],*pp.c[::-1,i]] for i in range(len(h))))
        r.presentation['formulas']=[('Công thức mỗi đoạn',[rf'S_i(x)=\sum_{{k=0}}^{{{degree}}}c_{{ik}}(x-x_i)^k'])]
        if target is not None:
            value=float(pp(target))
            if not np.isfinite(value):raise ValueError('Giá trị spline vượt miền float.')
            r.lines.append(f'S({fmt(target)}) = {fmt(value)}')
            r.presentation['formulas'].append(('Thay giá trị',[rf'S({fmt(target)})={fmt(value)}']))
        grid=np.linspace(x[0],x[-1],400);values=pp(grid)
        if not np.all(np.isfinite(values)):raise ValueError('Giá trị spline trên lưới vẽ vượt miền float; đổi đơn vị dữ liệu.')
        r.curves=[(grid,values,kind)]
        return r
    if kind in ('Cubic tự nhiên','Cubic kẹp','Cubic đạo hàm bậc hai'):
        from .math_view import cubic_spline_presentation
        cubic_spline_presentation(r,x,y,pp,kind,left,right,target)
        grid=np.linspace(x[0],x[-1],400); r.curves=[(grid,pp(grid),kind)]
        return r
    for i in range(len(h)):
        r.lines.append(f'[{fmt(x[i])}, {fmt(x[i+1])}]: '+polynomial_text(Poly(pp.c[:,i][::-1]),f'(x-({fmt(x[i])}))'))
        r.presentation['formulas'].append((f'Đoạn [{fmt(x[i])}, {fmt(x[i+1])}]',split_equation(rf'S_{{{i}}}(x)',poly_latex(Poly(pp.c[:,i][::-1]),rf'(x-({fmt(x[i])}))'))))
    if target is not None:
        r.lines.append(f'S({fmt(target)}) = {fmt(pp(target))}')
        r.presentation['formulas'].append(('Thay giá trị',[rf'S({fmt(target)})={fmt(pp(target))}']))
        if not x[0]<=target<=x[-1]: r.lines.append('Điểm tính ngoài bảng: ngoại suy bằng đoạn đầu/cuối.')
    grid=np.linspace(x[0],x[-1],400); r.curves=[(grid,pp(grid),kind)]
    return r


def least_squares(x,y,kind='Tuyến tính theo cơ sở',basis=('1','x'),weights=None):
    from .math_view import symbolic_latex,split_equation
    # Repeated abscissas are separate observations in regression, not duplicate
    # interpolation nodes. Dropping them changes the least-squares objective.
    x,y=np.asarray(x,dtype=float),np.asarray(y,dtype=float)
    if x.ndim!=1 or y.ndim!=1 or len(x)!=len(y) or len(x)<2:
        raise ValueError('Dãy x và y phải cùng độ dài, có ít nhất 2 điểm.')
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError('Dữ liệu phải hữu hạn.')
    indices=np.argsort(x,kind='stable');x,y=x[indices],y[indices]
    if kind not in ('Tuyến tính theo cơ sở','Lũy thừa','Hàm mũ'):
        raise ValueError('Mô hình bình phương tối thiểu không hợp lệ.')
    if weights is None: weights=np.ones(len(x))
    else:
        weights=np.asarray(weights,dtype=float)
        if weights.shape!=(len(x),) or not np.all(np.isfinite(weights)) or np.any(weights<=0):
            raise ValueError('Trọng số phải dương, hữu hạn và cùng độ dài với dữ liệu.')
        weights=weights[indices]
    sign=1.
    if kind!='Tuyến tính theo cơ sở':
        if np.any(y==0) or not (np.all(y>0) or np.all(y<0)):
            raise ValueError('Mô hình log cần mọi y khác 0 và cùng dấu; không tự dịch dữ liệu.')
        sign=float(np.sign(y[0])); transformed=np.log(abs(y))
    else: transformed=y
    if kind=='Lũy thừa':
        if np.any(x<=0): raise ValueError('Mô hình lũy thừa cần x > 0.')
        funcs=[lambda v: np.ones_like(v),lambda v: np.log(v)]
        labels=['1','ln(x)']
    else:
        labels=list(basis)
        if kind=='Hàm mũ': labels=['1']+[v for v in labels if v.strip()!='1']
        if not labels: raise ValueError('Cần ít nhất một hàm cơ sở.')
        funcs=[function(v) for v in labels]
    theta=np.column_stack([f(x) for f in funcs])
    if len(x)<theta.shape[1]: raise ValueError('Số điểm ít hơn số hàm cơ sở.')
    coeffs,high_precision=fit_basis(theta,transformed,weights)
    def pred(v):
        val=sum(c*f(np.asarray(v)) for c,f in zip(coeffs,funcs))
        result=val if kind=='Tuyến tính theo cơ sở' else sign*np.exp(val)
        if not np.all(np.isfinite(result)): raise ValueError('Dự đoán tràn số hoặc ngoài miền xác định.')
        return result
    text=' + '.join(f'({fmt(c)})*({s})' for c,s in zip(coeffs,labels))
    if kind!='Tuyến tính theo cơ sở': text=f'({fmt(sign)})*exp({text})'
    residual=np.asarray(pred(x)-y)
    if not np.all(np.isfinite(residual)):raise ValueError('Phần dư hồi quy vượt miền float; đổi đơn vị y.')
    sse=sum_products(residual,residual);mse=sse/len(x)
    transformed_residual=theta@coeffs-transformed
    if not np.all(np.isfinite(transformed_residual)):raise ValueError('Phần dư thang biến đổi vượt miền float; đổi tỉ lệ cơ sở.')
    transformed_sse=sum_products(transformed_residual,transformed_residual)
    transformed_objective=sum_products(weights,transformed_residual,transformed_residual)
    with localcontext() as context:
        context.prec=80
        centered=[Decimal(str(float(v)))-Decimal(str(float(y[0]))) for v in y]
        mean=sum(centered,Decimal(0))/len(y)
        total=sum(((v-mean)**2 for v in centered),Decimal(0))
    r=Report('Bình phương tối thiểu — '+kind,['f(x) = '+text,
        'SSE thang gốc = '+fmt(decimal_cell(sse))+'; MSE = '+fmt(decimal_cell(mse))+'; RMSE = '+fmt(decimal_cell(mse.sqrt())),
        'R² = '+(fmt(decimal_cell(1-sse/total)) if total>0 else 'không xác định (y hằng)'),
        'Tổng bình phương có trọng số = '+fmt(decimal_cell(sum_products(weights,residual,residual))),
        'SSE thang biến đổi = '+fmt(decimal_cell(transformed_sse))+'; mục tiêu có trọng số = '+fmt(decimal_cell(transformed_objective)) if kind!='Tuyến tính theo cơ sở' else 'Khớp trực tiếp trên thang y.',
        'Mô hình log cực tiểu hóa sai lệch log, không phải SSE trên thang y.' if kind!='Tuyến tính theo cơ sở' else '',
        'Giải bằng SVD sau đổi thang cột/trọng số; dùng QR Decimal cho cơ sở gần phụ thuộc.'+
        (' Đã dùng QR Decimal; kết quả vẫn nhạy với sai số dữ liệu.' if high_precision else '')],value=pred,scatter=(x,y))
    if kind=='Lũy thừa':
        r.lines.append(f'Tham số mô hình y = A*x^b: A = {fmt(exp_parameter(coeffs[0],sign))}; b = {fmt(coeffs[1])}.')
    elif kind=='Hàm mũ':
        r.lines.append(f'Tham số mô hình y = A*exp(tổng b_j*phi_j): A = {fmt(exp_parameter(coeffs[0],sign))}; b_j = '+', '.join(fmt(v) for v in coeffs[1:])+'.')
    r.table('Hệ số',['cơ sở','hệ số'],zip(labels,coeffs))
    r.presentation={'formulas':[('Hàm xấp xỉ',split_equation('f(x)',symbolic_latex(text))),
                                ('Bình phương tối thiểu',[r'\Theta^T W\Theta\,a=\Theta^T W z',r'\mathrm{MSE}=\frac{1}{N}\sum_{i=0}^{N-1}(f(x_i)-y_i)^2'])],
                    'omit_prefixes':('f(x) =',)}
    r.table('Ma trận thiết kế Θ',labels,theta)
    gram=[[decimal_cell(sum_products(weights,theta[:,i],theta[:,j])) for j in range(len(coeffs))] for i in range(len(coeffs))]
    r.table('Ma trận ΘᵀWΘ',[str(i) for i in range(len(coeffs))],gram)
    r.table('Vế phải ΘᵀWz',['cơ sở','giá trị'],zip(labels,[decimal_cell(sum_products(weights,column,transformed)) for column in theta.T]))
    r.table('Đối chiếu',['x','y','trọng số','dự đoán','sai lệch'],zip(x,y,weights,pred(x),residual))
    grid=np.linspace(x[0],x[-1],400)
    try:r.curves=[(grid,pred(grid),kind)]
    except ValueError:r.lines.append('Không vẽ toàn dải: mô hình không hữu hạn ở một số điểm giữa các quan sát. Các kết quả khớp tại dữ liệu vẫn được giữ.')
    return r
