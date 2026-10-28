"""Các thao tác Chương 1 dùng trực tiếp khi giải bài thi ngoại tuyến.

Mảng hệ số đa thức trong mô-đun này luôn theo bậc tăng.
"""
import math
import warnings
from decimal import Decimal, localcontext

import numpy as np
from numpy.polynomial import Polynomial
from scipy.interpolate import PPoly
from scipy.optimize import least_squares as nonlinear_least_squares
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve, MatrixRankWarning

from .common import Report, points, spacing, difference, fmt, function
from .polynomial_model import NewtonPolynomial, decimal_cell, nonnegative_integer, working_precision
from .linear_algebra import fit_basis, sum_products, exp_parameter


def _float_bound(value):
    """Làm tròn lên khi chuyển cận Decimal sang float; tránh báo cận 0 do underflow."""
    result=float(value)
    if not np.isfinite(result):raise ValueError('Cận sai số vượt miền float; đổi đơn vị dữ liệu.')
    if Decimal.from_float(result)<value:result=float(np.nextafter(result,np.inf))
    if not np.isfinite(result):raise ValueError('Cận sai số vượt miền float; đổi đơn vị dữ liệu.')
    return result


def chebyshev_value(degree, z):
    """A01: T_n(z) bằng truy hồi, dùng được cả ngoài [-1, 1]."""
    if not isinstance(degree,(int,np.integer)) or degree<0:
        raise ValueError('Bậc Chebyshev phải là số nguyên không âm.')
    if not np.isfinite(z): raise ValueError('Điểm tính phải hữu hạn.')
    previous, current=1.0, float(z)
    if degree==0: return previous
    for _ in range(1,degree):
        previous,current=current,2*z*current-previous
        if not np.isfinite(current):raise ValueError('T_n(z) tràn số; ngoài [-1,1] đa thức Chebyshev có thể tăng rất nhanh.')
    return current


def horner_division(coefficients, root):
    """A03: chia P cho x-root; hệ số đầu vào và thương theo bậc tăng."""
    coefficients=np.asarray(coefficients,float)
    if coefficients.ndim!=1 or len(coefficients)<2 or not np.all(np.isfinite(coefficients)) or not np.isfinite(root):
        raise ValueError('Cần đa thức bậc ít nhất một và dữ liệu hữu hạn.')
    quotient=np.empty(len(coefficients)-1)
    quotient[-1]=coefficients[-1]
    for j in range(len(quotient)-1,0,-1):
        quotient[j-1]=coefficients[j]+root*quotient[j]
    remainder=coefficients[0]+root*quotient[0]
    if not np.all(np.isfinite(quotient)) or not np.isfinite(remainder):
        raise ValueError('Chia Horner tràn số; đổi đơn vị đa thức.')
    return quotient,float(remainder)


def newton_horner(nodes, coefficients, target):
    """A13: đánh giá dạng Newton từ hệ số tỷ sai phân và đúng thứ tự mốc."""
    nodes=np.asarray(nodes,float); coefficients=np.asarray(coefficients,float)
    if nodes.ndim!=1 or coefficients.shape!=nodes.shape or not len(nodes) or not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(coefficients)) or not np.isfinite(target):
        raise ValueError('Mốc, hệ số Newton và điểm tính phải hữu hạn; hai mảng cùng độ dài.')
    if len(np.unique(nodes))!=len(nodes):raise ValueError('Mốc Newton phải phân biệt.')
    with localcontext() as ctx:
        ctx.prec=max(50,2*len(nodes)+30)
        value=Decimal(str(float(coefficients[-1])))
        for k in range(len(nodes)-2,-1,-1):
            value=Decimal(str(float(coefficients[k])))+(Decimal(str(float(target)))-Decimal(str(float(nodes[k]))))*value
        result=float(value)
    if not np.isfinite(result):raise ValueError('Giá trị Newton vượt miền float; chọn nhóm mốc nhỏ hơn hoặc đổi đơn vị.')
    return result


def newton_append(nodes, coefficients, new_x, new_y):
    """A15: thêm mốc vào dạng Newton, giữ thứ tự mốc gốc."""
    nodes=np.asarray(nodes,float); coefficients=np.asarray(coefficients,float)
    if nodes.ndim!=1 or coefficients.shape!=nodes.shape or not len(nodes) or not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(coefficients)) or not np.isfinite(new_x) or not np.isfinite(new_y):
        raise ValueError('Mốc và hệ số Newton phải cùng độ dài, hữu hạn.')
    if len(np.unique(nodes))!=len(nodes) or np.any(nodes==new_x):
        raise ValueError('Mốc Newton phải phân biệt.')
    with localcontext() as ctx:
        ctx.prec=max(50,2*len(nodes)+30)
        z=Decimal(str(float(new_x)))
        value=Decimal(str(float(coefficients[-1])))
        for k in range(len(nodes)-2,-1,-1):
            value=Decimal(str(float(coefficients[k])))+(z-Decimal(str(float(nodes[k]))))*value
        product=Decimal(1)
        for node in nodes:product*=z-Decimal(str(float(node)))
        new_coefficient=float((Decimal(str(float(new_y)))-value)/product)
    if not np.isfinite(new_coefficient):
        raise ValueError('Hệ số Newton bổ sung vượt miền float; đổi đơn vị mốc hoặc giữ dạng tích với ít mốc hơn.')
    return np.r_[nodes,new_x],np.r_[coefficients,new_coefficient]


def interpolation_error_bound(x, target, derivative_bound):
    """B01: M_N |prod(target-x_i)| / N!, với M_N được cung cấp bởi đề."""
    x, _, _ = points(x, np.zeros(len(x)))
    if not np.isfinite(target) or not np.isfinite(derivative_bound) or derivative_bound < 0:
        raise ValueError('Điểm tính và cận đạo hàm phải hữu hạn; cận không âm.')
    with localcontext() as ctx:
        ctx.prec=max(50,2*len(x)+30)
        value=Decimal(str(float(derivative_bound)))
        z=Decimal(str(float(target)))
        for k,node in enumerate(x,1):
            value*=abs(z-Decimal(str(float(node))))/k
        bound=_float_bound(value)
    lines = [f'Cận tại x = {fmt(target)}: {fmt(bound, 12)}.',
             'Cần chứng minh M_N chặn đạo hàm cấp N trên đoạn chứa mọi mốc và điểm tính.']
    if target < x[0] or target > x[-1]:
        lines.append('Ngoại suy: cận đạo hàm phải đúng cả trên đoạn mở rộng.')
    return Report('Cận sai số nội suy tại một điểm', lines, value=bound)


def interpolation_global_bound(x, derivative_bound, pattern):
    """Cận toàn đoạn cho lưới đều hoặc mốc Chebyshev khi đã biết M_N."""
    x,_,_=points(x,np.zeros(len(x)),2)
    if not np.isfinite(derivative_bound) or derivative_bound<0:
        raise ValueError('Cận đạo hàm phải hữu hạn, không âm.')
    n=len(x)
    if pattern=='cách đều':
        spacing(x)
        with localcontext() as ctx:
            ctx.prec=max(50,2*n+30)
            h=(Decimal(str(float(x[-1])))-Decimal(str(float(x[0]))))/(n-1)
            bound=_float_bound(Decimal(str(float(derivative_bound)))*h**n/(4*n))
    elif pattern=='Chebyshev':
        # Hai đầu khoảng [a,b] không trùng mốc Chebyshev loại nghiệm;
        # phải truyền thêm khoảng, không thể suy ra chính xác từ min/max mốc.
        raise ValueError('Cận Chebyshev toàn đoạn cần a,b; dùng chebyshev_global_bound.')
    else:raise ValueError('Chọn cách đều hoặc dùng chebyshev_global_bound.')
    return Report('Cận sai số nội suy toàn dải mốc',
                  ['Chỉ hợp lệ nếu M_N chặn |f^(N)| trên toàn dải mốc.'],value=bound)


def chebyshev_global_bound(count, a, b, derivative_bound):
    """Cận toàn [a,b] cho N nghiệm của đa thức Chebyshev T_N."""
    if not isinstance(count,(int,np.integer)) or count<1 or not all(np.isfinite(v) for v in (a,b,derivative_bound)) or a>=b or derivative_bound<0:
        raise ValueError('Cần N≥1; a<b; M_N hữu hạn không âm.')
    with localcontext() as ctx:
        ctx.prec=max(50,2*count+30)
        value=2*Decimal(str(float(derivative_bound)))
        width=(Decimal(str(float(b)))-Decimal(str(float(a))))/4
        for k in range(1,count+1):value*=width/k
        bound=_float_bound(value)
    return Report('Cận sai số Chebyshev trên toàn đoạn',
                  ['Cần dùng đúng N nghiệm Chebyshev trên [a,b] và M_N hợp lệ.'],value=bound)


def fill_missing(x, y, degree, tolerance=1e-9):
    """B02: điền các ô None/nan khi biết dữ liệu là đa thức bậc không quá degree."""
    x = np.asarray(x, float)
    if x.ndim != 1 or len(x) != len(y) or not np.all(np.isfinite(x)):
        raise ValueError('Dãy x hữu hạn và cùng độ dài với y.')
    if len(np.unique(x)) != len(x):
        raise ValueError('Mốc x trùng.')
    if not isinstance(degree, (int, np.integer)) or degree < 0 or degree >= len(x):
        raise ValueError('Bậc đa thức phải là số nguyên từ 0 đến N−1.')
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('Dung sai kiểm tra phải dương, hữu hạn.')
    yy = np.asarray([np.nan if v is None else v for v in y], float)
    known = np.isfinite(yy)
    if np.any(np.isinf(yy)) or np.sum(known) < degree + 1:
        raise ValueError('Không đủ d+1 tung độ đã biết hữu hạn để xác định đa thức.')
    selected = np.flatnonzero(known)[:degree + 1]
    interpolant = NewtonPolynomial(x[selected], yy[selected])
    predicted = np.asarray(interpolant(x), float)
    scale = float(np.max(np.abs(yy[known])))
    if np.max(np.abs(predicted[known] - yy[known])) > tolerance * scale:
        raise ValueError('Dữ liệu đã biết mâu thuẫn với giả thiết bậc đa thức.')
    result = yy.copy()
    result[~known] = predicted[~known]
    return Report('Điền dữ liệu theo giả thiết bậc đa thức',
                  [f'Đã điền {np.sum(~known)} ô; giả thiết bậc ≤ {degree}.',
                   'Kết quả phụ thuộc vào giả thiết bậc được cho.'], value=result).table(
                       'Bảng hoàn chỉnh', ['x', 'y'], zip(x, result))


def scan_level(x, y, target):
    """A21: nghiệm tại mốc và các khoảng đổi dấu nghiêm ngặt trên bảng."""
    x, y, _ = points(x, y, 2)
    if not np.isfinite(target): raise ValueError('Mức y phải hữu hạn.')
    exact = x[y == target].tolist()
    intervals = [(float(x[i]), float(x[i+1])) for i in range(len(x)-1)
                 if (y[i]-target)*(y[i+1]-target) < 0]
    flats = [(float(x[i]), float(x[i+1])) for i in range(len(x)-1)
             if y[i] == target and y[i+1] == target]
    return Report('Quét mức trên bảng',
                  [f'Nghiệm tại mốc: {exact}; khoảng đổi dấu: {intervals}; đoạn bằng mức: {flats}.',
                   'Chỉ quét bảng; có thể bỏ sót nghiệm tiếp xúc giữa hai mốc.',
                   'Tính duy nhất của nghiệm thật cần giả thiết bổ sung.'],
                  value=(exact, intervals, flats))


def monotone_branches(x, y):
    """A22: trả các nhánh chỉ số đơn điệu và các đoạn phẳng của bảng."""
    x,y,_=points(x,y,2)
    branches=[]; flats=[]; start=0; direction=0
    for i in range(len(x)-1):
        step=int(np.sign(y[i+1]-y[i]))
        if step==0:
            if start<i:branches.append((start,i))
            flats.append((i,i+1)); start=i+1; direction=0
        elif direction==0:
            direction=step
        elif step!=direction:
            branches.append((start,i)); start=i; direction=step
    if start<len(x)-1:branches.append((start,len(x)-1))
    return Report('Nhánh đơn điệu trên bảng',
                  [f'Nhánh chỉ số: {branches}; đoạn phẳng: {flats}.',
                   'Các nhánh chỉ phản ánh dãy mẫu; cần giả thiết riêng cho tính đơn điệu của hàm thật.'],
                  value=(branches,flats))


def inverse_by_swap(x, y, target, count=None):
    """A23: nội suy x như hàm của y trên một nhánh mẫu đơn điệu."""
    x, y, _ = points(x, y, 2)
    if not (np.all(np.diff(y) > 0) or np.all(np.diff(y) < 0)):
        raise ValueError('Cần chọn riêng một nhánh có y đơn điệu nghiêm ngặt.')
    if not np.isfinite(target) or not min(y) <= target <= max(y):
        raise ValueError('Mức y nằm ngoài nhánh đã chọn.')
    count = len(x) if count is None else count
    if not isinstance(count, (int, np.integer)) or not 2 <= count <= len(x):
        raise ValueError('Số mốc đảo phải từ 2 đến kích thước nhánh.')
    chosen = np.sort(np.argsort(abs(y-target), kind='stable')[:count])
    estimate = float(NewtonPolynomial(y[chosen], x[chosen])(target))
    forward = float(NewtonPolynomial(x[chosen], y[chosen])(estimate))
    lines = [f'x xấp xỉ = {fmt(estimate, 12)}; phần dư đa thức thuận = {fmt(forward-target, 12)}.',
             'Đây là nghiệm của đa thức nội suy bảng đảo; nói chung khác nghiệm đa thức thuận.']
    if not x[chosen[0]] <= estimate <= x[chosen[-1]]:
        lines.append('CẢNH BÁO: nghiệm mô hình đảo nằm ngoài nhóm mốc đã chọn.')
    return Report('Nội suy hàm ngược', lines, value=estimate).table(
        'Mốc đảo được chọn', ['x', 'y'], zip(x[chosen], y[chosen]))


def bisect_model(model, target, a, b, tolerance=1e-10, max_steps=1000):
    """A25: chia đôi P(x)=target; trả nghiệm và khoảng bao nghiệm mô hình."""
    if not all(np.isfinite(v) for v in (target, a, b, tolerance)) or a >= b or tolerance <= 0:
        raise ValueError('Cần a < b và dung sai dương, mọi tham số hữu hạn.')
    if not isinstance(max_steps,(int,np.integer)) or max_steps<1:
        raise ValueError('Số bước chia đôi tối đa phải là số nguyên dương.')
    def residual(z):
        value = float(model(z)) - target
        if not np.isfinite(value): raise ValueError('Mô hình không hữu hạn trong khoảng chia đôi.')
        return value
    fa, fb = residual(a), residual(b)
    if fa == 0: return Report('Chia đôi nghiệm mô hình',[f'Nghiệm đúng tại đầu trái: x = {fmt(a)}.'],value=(a, (a, a)))
    if fb == 0: return Report('Chia đôi nghiệm mô hình',[f'Nghiệm đúng tại đầu phải: x = {fmt(b)}.'],value=(b, (b, b)))
    if np.sign(fa) == np.sign(fb):
        raise ValueError('Hai đầu không trái dấu; chưa đủ điều kiện chia đôi.')
    for step in range(max_steps):
        c = a/2 + b/2
        fc = residual(c)
        if fc == 0: a = b = c; break
        if c==a or c==b:
            if b/2-a/2<=tolerance:break
            raise ValueError('Dung sai nhỏ hơn độ phân giải float của khoảng này; đổi gốc/đơn vị hoặc tăng dung sai.')
        if np.sign(fa) != np.sign(fc): b, fb = c, fc
        else: a, fa = c, fc
        if b/2-a/2 <= tolerance: break
    else:
        raise ValueError('Chưa đạt dung sai sau số bước tối đa.')
    root = a/2+b/2
    return Report('Chia đôi nghiệm mô hình',
                  [f'x ≈ {fmt(root, 12)}; khoảng [{fmt(a, 12)}, {fmt(b, 12)}].',
                   f'Cận vị trí nghiệm mô hình ≤ {fmt(b/2-a/2, 12)}.',
                   'Sai số so với nghiệm hàm thật cần thêm cận sai số mô hình và độ dốc.'],
                  value=(root, (a, b)))


def inverse_controlled(x, y, target, interval, backward=False,
                       tolerance=1e-10, residual_tolerance=1e-10, max_steps=1000):
    """A24: chỉ báo cận vị trí khi đã chặn được ánh xạ co trên J.

    Chứng nhận dùng cận hệ số đa thức trên [-1,1], thường bảo thủ. Cận áp dụng
    cho nghiệm của đa thức nội suy tính bằng số thực dấu phẩy động.
    """
    x,y,_=points(x,y,2)
    h=spacing(x)
    if len(interval)!=2 or not all(np.isfinite(v) for v in (target,*interval,tolerance,residual_tolerance)):
        raise ValueError('Mức, khoảng và dung sai phải hữu hạn.')
    left,right=map(float,interval)
    if not x[0]<=left<right<=x[-1] or tolerance<=0 or residual_tolerance<=0:
        raise ValueError('J phải nằm trong dải mốc; dung sai dương.')
    if not isinstance(max_steps,(int,np.integer)) or max_steps<1:
        raise ValueError('Số bước tối đa phải dương.')
    at_node=np.flatnonzero((y==target)&(x>=left)&(x<=right))
    if len(at_node):return Report('Nội suy ngược — mức tại mốc',[f'Mức bằng tung độ tại mốc x = {fmt(x[at_node[0]])}.'],value=float(x[at_node[0]]))
    from .inverse_model import InverseMap, compose_interval
    precision=max(50,-Decimal(str(float(tolerance))).adjusted()+40,-Decimal(str(float(residual_tolerance))).adjusted()+40)
    model=InverseMap(x,y,target,backward,precision)
    rows=[]
    with localcontext() as ctx:
        ctx.prec=model.precision
        lo=(Decimal(str(left))-model.origin)/model.h;hi=(Decimal(str(right))-model.origin)/model.h
        center=(lo+hi)/2;radius=(hi-lo)/2
        change=compose_interval(model.phi,center,radius)
        deviation=sum((abs(v) for v in change[1:]),Decimal(0))
        center_value=change[0]
        q=sum((k*abs(v) for k,v in enumerate(change)),Decimal(0))/radius
        if q>=1 or center_value-deviation<lo or center_value+deviation>hi:
            raise ValueError('Chưa chứng minh được phi(J) nằm trong J và |phi′|<1; chọn J khác hoặc chia đôi.')
        t=model.initial()
        if not lo<=t<=hi:t=center
        dtol=Decimal(str(float(tolerance)));rtol=Decimal(str(float(residual_tolerance)))
        for k in range(1,max_steps+1):
            next_t=model.step(t)
            if not lo<=next_t<=hi:raise ValueError('Bước lặp ra khỏi J.')
            residual=model.residual(next_t)
            bound=model.h*q/(1-q)*abs(next_t-t)
            exact=model.location(next_t);value=float(exact)
            rounding=abs(Decimal.from_float(value)-exact)
            total_bound=bound+rounding
            rows.append((k,decimal_cell(t),decimal_cell(next_t),decimal_cell(residual),decimal_cell(total_bound)))
            if total_bound<=dtol and residual<=rtol:
                return Report('Nội suy ngược — lặp có cận co',
                              [f'x ≈ {fmt(value,12)}; cận nghiệm mô hình ≤ {fmt(decimal_cell(total_bound),12)}.',
                               f'Đã kiểm tra trên J=[{fmt(left)},{fmt(right)}]: q ≤ {fmt(decimal_cell(q),12)} < 1.',
                               f'Phần dư |P(x)-y*| = {fmt(decimal_cell(residual),12)}.',
                               'Cận gồm sai số lặp và chuyển nghiệm sang float; chưa chứng nhận làm tròn hệ số bằng số học khoảng.',
                               'Cận này chỉ xét đa thức nội suy; sai số hàm thật cần cận mô hình và độ dốc.'],
                              value=value).table('Lịch sử lặp',
                                  ['k','t trước','t sau','phần dư','cận vị trí x'],rows)
            if bound<=dtol and residual<=rtol and rounding>dtol:
                raise ValueError('Dung sai vị trí nhỏ hơn độ phân giải float; đổi gốc/đơn vị hoặc tăng dung sai.')
            if next_t==t:break
            t=next_t
    raise ValueError('Hết số bước lặp mà chưa đạt cả hai dung sai.')


def spline_degree(x, y, degree, left_derivatives=()):
    """A29: spline bậc 1–4, biết p−1 đạo hàm ở đầu trái."""
    x, y, _ = points(x, y, 2)
    if not isinstance(degree, (int, np.integer)) or not 1 <= degree <= 4:
        raise ValueError('Bậc spline phải là 1, 2, 3 hoặc 4.')
    v = np.asarray(left_derivatives, float)
    if v.shape != (degree-1,) or not np.all(np.isfinite(v)):
        raise ValueError('Cần đúng p−1 đạo hàm đầu trái hữu hạn.')
    columns = []
    with localcontext() as ctx:
        ctx.prec=max(80,working_precision(x))
        derivatives=[Decimal(str(float(value))) for value in v]
        for i in range(len(x)-1):
            h=Decimal(str(float(x[i+1])))-Decimal(str(float(x[i])))
            coefficients=[Decimal(str(float(y[i])))]+[derivatives[k-1]/math.factorial(k) for k in range(1,degree)]
            known=sum((value*h**k for k,value in enumerate(coefficients)),Decimal(0))
            coefficients.append((Decimal(str(float(y[i+1])))-known)/h**degree)
            columns.append([float(value) for value in reversed(coefficients)])
            derivatives=[sum((coefficients[k]*(math.factorial(k)//math.factorial(k-order))*h**(k-order)
                              for k in range(order,degree+1)),Decimal(0)) for order in range(1,degree)]
    if not np.all(np.isfinite(columns)):
        raise ValueError('Hệ số spline tràn số; đổi đơn vị mốc/đạo hàm biên.')
    pp = PPoly(np.array(columns).T, x, extrapolate=False)
    return Report(f'Spline bậc {degree} — đạo hàm đầu trái',
                  ['Mỗi cột hệ số dùng u = x - mốc trái; thứ tự trong bảng là bậc tăng.'],
                  value=pp).table('Hệ số từng đoạn',
                                  ['đoạn']+[f'u^{k}' for k in range(degree+1)],
                                  [[i]+list(row[::-1]) for i, row in enumerate(columns)])


def thomas(lower, diagonal, upper, rhs, pivot_tolerance=1e-12):
    """A28: truy đuổi hệ ba đường chéo, trả nghiệm và phần dư hệ gốc."""
    lower=np.asarray(lower,float);diagonal=np.asarray(diagonal,float)
    upper=np.asarray(upper,float);rhs=np.asarray(rhs,float)
    n=len(diagonal)
    if n<1 or any(v.shape!=(n,) or not np.all(np.isfinite(v)) for v in (lower,diagonal,upper,rhs)):
        raise ValueError('Bốn dãy ba đường chéo cùng độ dài và hữu hạn.')
    if lower[0]!=0 or upper[-1]!=0 or not np.isfinite(pivot_tolerance) or pivot_tolerance<=0:
        raise ValueError('Cần lower[0]=upper[-1]=0 và dung sai pivot dương.')
    row_scale=np.maximum.reduce([abs(lower),abs(diagonal),abs(upper)])
    if np.any(row_scale==0):raise ValueError('Có phương trình toàn hệ số 0; hệ không xác định duy nhất.')
    original_lower=lower.copy();original_upper=upper.copy()
    lower=lower/row_scale;upper=upper/row_scale
    b=diagonal/row_scale;d=rhs/row_scale;rows=[]
    if not np.all(np.isfinite(d)):raise ValueError('Vế phải sau đổi thang vượt miền float.')
    scale=max(float(np.max(abs(b))),float(np.max(abs(lower))),float(np.max(abs(upper))))
    for i in range(1,n):
        if abs(b[i-1])<=pivot_tolerance*scale:raise ValueError('Pivot Thomas gần 0; cần bộ giải có chọn trụ.')
        factor=lower[i]/b[i-1]
        b[i]-=factor*upper[i-1];d[i]-=factor*d[i-1]
        rows.append((i,factor,b[i],d[i]))
    if abs(b[-1])<=pivot_tolerance*scale:raise ValueError('Pivot Thomas gần 0; cần bộ giải có chọn trụ.')
    solution=np.empty(n);solution[-1]=d[-1]/b[-1]
    for i in range(n-2,-1,-1):
        solution[i]=(d[i]-upper[i]*solution[i+1])/b[i]
    if not np.all(np.isfinite(solution)):
        raise ValueError('Giải Thomas tràn số; đổi tỉ lệ hệ hoặc dùng bộ giải có chọn trụ.')
    residual=[]
    for i in range(n):
        coefficients=[diagonal[i]];values=[solution[i]]
        if i:coefficients.append(original_lower[i]);values.append(solution[i-1])
        if i<n-1:coefficients.append(original_upper[i]);values.append(solution[i+1])
        residual.append(sum_products(coefficients,values)-Decimal(str(float(rhs[i]))))
    return Report('Thomas — hệ ba đường chéo',
                  [f'Phần dư lớn nhất = {fmt(decimal_cell(max(map(abs,residual))),12)}.',
                   'Đã chia mỗi phương trình cho hệ số có độ lớn lớn nhất của hàng trước khi khử.'],
                  value=solution).table('Khử tiến',
                      ['i','hệ số khử','đường chéo sau khử','vế phải sau khử'],rows)


def gaussian_elimination(matrix, rhs, pivot_tolerance=1e-12):
    """A30: khử Gauss có chọn pivot theo cột, rồi thế ngược."""
    a=np.asarray(matrix,float);b=np.asarray(rhs,float)
    if a.ndim!=2 or a.shape[0]<1 or a.shape[0]!=a.shape[1] or b.shape!=(a.shape[0],) or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError('Cần ma trận vuông và vế phải cùng cỡ, hữu hạn.')
    if not np.isfinite(pivot_tolerance) or pivot_tolerance<=0:raise ValueError('Dung sai pivot phải dương.')
    original=a.copy();original_rhs=b.copy();n=len(b);rows=[]
    row_scale=np.max(abs(a),axis=1)
    if np.any(row_scale==0):raise ValueError('Có phương trình toàn hệ số 0; hệ không xác định duy nhất.')
    a=a/row_scale[:,None];b=b/row_scale
    if not np.all(np.isfinite(b)):raise ValueError('Vế phải sau đổi thang vượt miền float.')
    scale=float(np.max(abs(a)))
    for k in range(n):
        pivot=k+int(np.argmax(abs(a[k:,k])))
        if abs(a[pivot,k])<=pivot_tolerance*scale:
            raise ValueError('Hệ chưa xác định duy nhất trong dung sai; cần xét hạng.')
        if pivot!=k:
            a[[k,pivot]]=a[[pivot,k]];b[[k,pivot]]=b[[pivot,k]]
        rows.append((k,pivot,float(a[k,k])))
        for i in range(k+1,n):
            factor=a[i,k]/a[k,k]
            a[i,k:]-=factor*a[k,k:];b[i]-=factor*b[k]
    solution=np.zeros(n)
    for i in range(n-1,-1,-1):
        solution[i]=(b[i]-a[i,i+1:]@solution[i+1:])/a[i,i]
    if not np.all(np.isfinite(solution)):
        raise ValueError('Khử Gauss tràn số; đổi tỉ lệ hệ trước khi giải.')
    residual=[sum_products(row,solution)-Decimal(str(float(v))) for row,v in zip(original,original_rhs)]
    return Report('Khử Gauss có chọn trụ',
                  [f'Phần dư lớn nhất = {fmt(decimal_cell(max(map(abs,residual))),12)}.',
                   'Đã đổi thang từng hàng trước khi chọn trụ và khử.'],
                  value=solution).table('Chọn pivot',
                      ['bước','hàng được chọn trong hệ đang khử','giá trị pivot'],rows)


def spline_boundary(x, y, degree, conditions):
    """B03: p−1 điều kiện biên đơn hoặc tổ hợp tuyến tính đạo hàm.

    Điều kiện đơn: ('trái'/'phải', cấp, giá trị).
    Tổ hợp: ([(phía, cấp, hệ số), ...], giá trị).
    """
    x, y, _ = points(x, y, 2)
    if not isinstance(degree, (int, np.integer)) or not 1 <= degree <= 4:
        raise ValueError('Bậc spline phải từ 1 đến 4.')
    if len(conditions) != degree-1:
        raise ValueError('Cần đúng p−1 điều kiện biên độc lập.')
    simple=[]
    for condition in conditions:
        if len(condition)!=3 or not isinstance(condition[0],str):break
        side,order,value=condition
        if side not in ('trái','phải') or not isinstance(order,(int,np.integer)) or not 1<=order<degree or not np.isfinite(value):break
        simple.append((side,int(order),float(value)))
    if len(simple)==degree-1 and len({side for side,_,_ in simple})<=1 and {order for _,order,_ in simple}==set(range(1,degree)):
        derivatives=[next(value for _,order,value in simple if order==k) for k in range(1,degree)]
        if not simple or simple[0][0]=='trái':
            result=spline_degree(x,y,degree,derivatives)
            result.title=f'Spline bậc {degree} — biên hai đầu'
            result.lines.append('Các đạo hàm cùng ở đầu trái: dùng truy hồi Decimal tương đương, tránh hệ float quá nhạy trên bảng dài.')
            return result
        # Đổi t=-x để dùng cùng truy hồi ổn định chữ số với dữ liệu đầu phải.
        reverse=spline_degree(-x[::-1],y[::-1],degree,[(-1)**k*value for k,value in enumerate(derivatives,1)]).value
        columns=[]
        with localcontext() as ctx:
            ctx.prec=max(80,working_precision(x))
            for i in range(len(x)-1):
                h=Decimal(str(float(x[i+1])))-Decimal(str(float(x[i])))
                old=[Decimal(str(float(v))) for v in reverse.c[::-1,len(x)-2-i]]
                new=[float((-1)**k*sum((old[j]*math.comb(j,k)*h**(j-k) for j in range(k,degree+1)),Decimal(0))) for k in range(degree+1)]
                columns.append(new)
        if not np.all(np.isfinite(columns)):raise ValueError('Hệ số spline theo x vượt miền float; đổi đơn vị.')
        pp=PPoly(np.array(columns)[:,::-1].T,x,extrapolate=False)
        return Report(f'Spline bậc {degree} — biên hai đầu',
                      ['Các đạo hàm cùng ở đầu phải: đổi t=-x và dùng truy hồi Decimal.',
                       'Hệ số theo u=x−mốc trái, bậc tăng.'],value=pp).table(
                          'Hệ số từng đoạn',['đoạn']+[f'u^{k}' for k in range(degree+1)],
                          ([i,*row] for i,row in enumerate(columns)))
    n = len(x)-1; p = degree; size = n*(p+1); steps=np.diff(x)
    if not np.all(np.isfinite(steps)):raise ValueError('Độ dài đoạn spline vượt miền float; đổi đơn vị x.')
    row_indices=[];column_indices=[];entries=[];rhs=[]
    def derivative_row(segment, fraction, order):
        # Ẩn theo v=(x-x_i)/h_i. Chuẩn hóa từng phương trình trước khi
        # đưa vào bộ giải thưa để đơn vị h không làm sai việc xét hạng.
        with localcontext() as ctx:
            ctx.prec=80
            h=Decimal(str(float(steps[segment])))
            return {segment*(p+1)+power:Decimal(math.factorial(power)//math.factorial(power-order))/h**order
                    for power in range(order,p+1) if fraction or power==order}
    def append_row(terms,value):
        scale=max((abs(v) for v in terms.values()),default=Decimal(0))
        if not scale:raise ValueError('Điều kiện spline không độc lập hoặc không tương thích.')
        row=len(rhs)
        for column,value_a in terms.items():
            if value_a:
                row_indices.append(row);column_indices.append(column);entries.append(float(value_a/scale))
        rhs.append(float(Decimal(str(float(value)))/scale))
    for i in range(n):
        append_row(derivative_row(i,0,0),y[i])
        append_row(derivative_row(i,1,0),y[i+1])
    for i in range(n-1):
        for order in range(1,p):
            terms=derivative_row(i,1,order)
            for column,value in derivative_row(i+1,0,order).items():terms[column]=-value
            append_row(terms,0.)
    for condition in conditions:
        if len(condition)==3 and isinstance(condition[0],str):
            side,order,value=condition
            terms=[(side,order,1.)]
        elif len(condition)==2:
            terms,value=condition
        else:raise ValueError('Điều kiện biên phải là đạo hàm đơn hoặc tổ hợp tuyến tính.')
        if not np.isfinite(value) or not terms:
            raise ValueError('Điều kiện biên cần hệ số và vế phải hữu hạn.')
        boundary_row={}
        for side,order,weight in terms:
            if side not in ('trái','phải') or not isinstance(order,(int,np.integer)) or not 0<=order<=p or not np.isfinite(weight):
                raise ValueError('Mỗi hạng biên là (trái/phải, cấp 0..p, hệ số hữu hạn).')
            segment,fraction=(0,0) if side=='trái' else (n-1,1)
            for column,value_a in derivative_row(segment,fraction,order).items():
                boundary_row[column]=boundary_row.get(column,Decimal(0))+Decimal(str(float(weight)))*value_a
        append_row(boundary_row,value)
    matrix=coo_matrix((entries,(row_indices,column_indices)),shape=(size,size)).tocsc()
    rhs=np.asarray(rhs)
    if not np.all(np.isfinite(rhs)):raise ValueError('Vế phải spline vượt miền float sau chuẩn hóa; đổi đơn vị.')
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('error',MatrixRankWarning)
            normalized=spsolve(matrix,rhs)
    except (MatrixRankWarning,RuntimeError) as exc:
        raise ValueError('Các điều kiện không xác định spline duy nhất trong độ chính xác máy.') from exc
    if not np.all(np.isfinite(normalized)):raise ValueError('Hệ spline suy biến hoặc tràn số.')
    residual=float(np.max(abs(matrix@normalized-rhs)))
    if residual>1e-9*max(1.,float(np.max(abs(rhs)))):
        raise ValueError('Hệ spline không đạt phần dư yêu cầu; kiểm tra điều kiện biên và thang dữ liệu.')
    normalized=normalized.reshape(n,p+1)
    with localcontext() as ctx:
        ctx.prec=80
        coefficients=np.array([[float(Decimal(str(float(v)))/Decimal(str(float(h)))**k)
                                for k,v in enumerate(row)] for h,row in zip(steps,normalized)])
    if not np.all(np.isfinite(coefficients)):raise ValueError('Hệ số spline theo x vượt miền float; đổi đơn vị x.')
    pp=PPoly(coefficients[:,::-1].T,x,extrapolate=False)
    return Report(f'Spline bậc {p} — biên hai đầu',
                  [f'Phần dư hệ phương trình sau chuẩn hóa: {fmt(residual,12)}.',
                   'Hệ số từng mảnh theo bậc tăng của u = x - mốc trái.'],
                  value=pp).table('Hệ số từng đoạn',
                                  ['đoạn']+[f'u^{k}' for k in range(p+1)],
                                  [[i]+list(row) for i,row in enumerate(coefficients)])


def spline_at(pp, target, order=0):
    """B04: kiểm tra đạo hàm hai phía ở nút trước khi báo giá trị chung."""
    x=pp.x
    if not np.isfinite(target) or target<x[0] or target>x[-1]:
        raise ValueError('Điểm tính nằm ngoài miền spline; không tự ngoại suy.')
    if not isinstance(order,(int,np.integer)) or order<0:
        raise ValueError('Cấp đạo hàm phải là số nguyên không âm.')
    i=int(np.searchsorted(x,target,side='right')-1)
    i=min(i,len(x)-2)
    degree=pp.c.shape[0]-1
    def piece(segment,k):
        if k>degree:return Decimal(0),Decimal(0)
        offset=Decimal(str(float(target)))-Decimal(str(float(x[segment])))
        result=Decimal(0);magnitude=Decimal(0)
        for power in range(degree,k-1,-1):
            coefficient=Decimal(str(float(pp.c[degree-power,segment])))*(math.factorial(power)//math.factorial(power-k))
            result=result*offset+coefficient;magnitude=magnitude*abs(offset)+abs(coefficient)
        return result,magnitude
    with localcontext() as ctx:
        ctx.prec=80
        value,_=piece(i,order)
        if 0<i<len(x)-1 and target==x[i]:
            for k in range(min(order,degree)+1):
                left,lmag=piece(i-1,k);right,rmag=piece(i,k)
                tolerance=Decimal('1e-8')*max(abs(left),abs(right))+Decimal(str(128*np.finfo(float).eps))*max(lmag,rmag)
                if abs(left-right)>tolerance:
                    raise ValueError(f'Đạo hàm cấp {k} không chung tại nút: trái={fmt(decimal_cell(left))}, phải={fmt(decimal_cell(right))}.')
    value=float(value)
    if not np.isfinite(value):raise ValueError('Giá trị spline/đạo hàm không hữu hạn.')
    return value


def transformed_fit(x, y, model):
    """A32: các dạng biến đổi còn lại của bảng mô hình trang 49.

    model: '1/(a+b*x)', 'x/(a+b*x)', 'a+b*ln(x)', 'A*b^x'.
    """
    x=np.asarray(x,float); y=np.asarray(y,float)
    if x.ndim!=1 or y.shape!=x.shape or len(x)<2 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError('Cần ít nhất hai cặp dữ liệu hữu hạn.')
    if model=='1/(a+b*x)':
        if np.any(y==0):raise ValueError('Phép nghịch đảo cần mọi y khác 0.')
        u=x; v=1/y; names=('a','b')
    elif model=='x/(a+b*x)':
        if np.any(x==0) or np.any(y==0):raise ValueError('Mô hình cần x và y khác 0.')
        u=1/x; v=1/y; names=('b','a')
    elif model=='a+b*ln(x)':
        if np.any(x<=0):raise ValueError('Log cần x > 0.')
        u=np.log(x); v=y; names=('a','b')
    elif model=='A*b^x':
        if np.any(y==0) or not (np.all(y>0) or np.all(y<0)):
            raise ValueError('Log cần mọi y cùng dấu và khác 0.')
        u=x; v=np.log(abs(y)); names=('ln|A|','ln(b)')
    else:
        raise ValueError('Tên mô hình biến đổi không hợp lệ.')
    theta=np.column_stack((np.ones(len(x)),u))
    if not np.all(np.isfinite(theta)) or not np.all(np.isfinite(v)):
        raise ValueError('Biến đổi log/nghịch đảo tràn số; đổi đơn vị dữ liệu.')
    coefficients,high_precision=fit_basis(theta,v)
    c0,c1=map(float,coefficients)
    if model=='1/(a+b*x)':
        params={'a':c0,'b':c1}
        def predict(z):return 1/(c0+c1*np.asarray(z,float))
    elif model=='x/(a+b*x)':
        params={'a':c1,'b':c0}
        def predict(z):
            z=np.asarray(z,float)
            return z/(c1+c0*z)
    elif model=='a+b*ln(x)':
        params={'a':c0,'b':c1}
        def predict(z):return c0+c1*np.log(np.asarray(z,float))
    else:
        params={'A':exp_parameter(c0,float(np.sign(y[0]))),'b':exp_parameter(c1)}
        def predict(z):return np.sign(y[0])*np.exp(c0+c1*np.asarray(z,float))
    def checked_predict(z):
        if not np.all(np.isfinite(z)):raise ValueError('Điểm dự đoán phải hữu hạn.')
        with np.errstate(all='ignore'):result=np.asarray(predict(z),float)
        if not np.all(np.isfinite(result)):
            raise ValueError('Dự đoán ngoài miền xác định hoặc tràn số.')
        return result
    prediction=checked_predict(x)
    residual=prediction-y;transformed_residual=theta@coefficients-v
    if not np.all(np.isfinite(residual)) or not np.all(np.isfinite(transformed_residual)):
        raise ValueError('Phần dư hồi quy vượt miền float; đổi đơn vị dữ liệu.')
    sse=sum_products(residual,residual);transformed_sse=sum_products(transformed_residual,transformed_residual)
    report=Report('Khớp mô hình biến đổi — '+model,
                  [f'SSE thang y = {fmt(decimal_cell(sse),12)}; MSE = {fmt(decimal_cell(sse/len(x)),12)}; RMSE = {fmt(decimal_cell((sse/len(x)).sqrt()),12)}.',
                   f'SSE thang biến đổi = {fmt(decimal_cell(transformed_sse),12)}.',
                   'Tham số cực tiểu hóa sai lệch ở thang biến đổi; không mặc nhiên cực tiểu SSE thang y.'],
                  value=checked_predict,scatter=(x,y))
    report.table('Tham số',['tên','giá trị'],params.items())
    report.table('Bảng biến đổi và dự đoán',['x','y','u','v','dự đoán'],zip(x,y,u,v,prediction))
    return report


def original_scale_fit(x, y, model='Hàm mũ', basis=('x',), initial=None):
    """Tối thiểu hóa SSE ở thang y cho A*exp(sum b_j*psi_j(x)).

    Với model='Lũy thừa', psi_1(x)=ln(x), tức y=A*x^b.
    Bộ giải là lặp phi tuyến; báo trạng thái và phần dư, không bảo đảm cực tiểu toàn cục.
    """
    x=np.asarray(x,float);y=np.asarray(y,float)
    if x.ndim!=1 or y.shape!=x.shape or len(x)<2 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError('Cần ít nhất hai cặp dữ liệu hữu hạn.')
    if model=='Lũy thừa':
        if np.any(x<=0):raise ValueError('Lũy thừa với số mũ thực cần x>0.')
        labels=['ln(x)']; funcs=[lambda z:np.log(z)]
    elif model=='Hàm mũ':
        labels=list(basis)
        if not labels:raise ValueError('Cần ít nhất một cơ sở trong số mũ.')
        funcs=[function(label) for label in labels]
    else:raise ValueError('Chọn Hàm mũ hoặc Lũy thừa.')
    theta=np.column_stack([f(x) for f in funcs])
    design=np.column_stack([np.ones(len(x)),theta])
    fit_basis(design,np.zeros(len(x)))
    yscale=float(np.max(abs(y))) or 1.
    scaled_y=y/yscale
    if initial is None:
        if np.all(y>0) or np.all(y<0):
            coefficients,_=fit_basis(design,np.log(abs(y))-np.log(yscale))
            initial=np.r_[np.sign(y[0])*np.exp(coefficients[0]),coefficients[1:]]
        else:initial=np.r_[float(np.mean(scaled_y)) if np.mean(scaled_y)!=0 else 1.,np.zeros(theta.shape[1])]
    else:
        initial=np.asarray(initial,float).copy()
        if initial.ndim==1 and len(initial):initial[0]/=yscale
    initial=np.asarray(initial,float)
    if initial.shape!=(theta.shape[1]+1,) or not np.all(np.isfinite(initial)):
        raise ValueError('Giá trị đầu phải gồm A và các b_j hữu hạn.')
    def prediction(parameters, values):
        with np.errstate(over='ignore',invalid='ignore'):
            result=np.zeros(len(values)) if parameters[0]==0 else np.sign(parameters[0])*np.exp(np.log(abs(parameters[0]))+values@parameters[1:])
        return result
    def residual(parameters):
        predicted=prediction(parameters,theta)
        if not np.all(np.isfinite(predicted)):
            raise ValueError('Tham số lặp làm mô hình tràn số.')
        return predicted-scaled_y
    try:solution=nonlinear_least_squares(residual,initial,max_nfev=10000,
                                          xtol=1e-12,ftol=1e-12,gtol=1e-12)
    except ValueError as exc:raise ValueError('Tối ưu phi tuyến không tính được tại giá trị đầu; chọn lại tham số đầu.') from exc
    if not solution.success or not np.all(np.isfinite(solution.x)):
        raise ValueError('Tối ưu phi tuyến chưa hội tụ.')
    parameters=solution.x
    def predict(z):
        z=np.asarray(z,float)
        if not np.all(np.isfinite(z)):raise ValueError('Điểm dự đoán phải hữu hạn.')
        if model=='Lũy thừa' and np.any(z<=0):raise ValueError('Dự đoán lũy thừa cần x>0.')
        features=np.column_stack([f(z.ravel()) for f in funcs])
        with np.errstate(over='ignore',invalid='ignore'):
            result=np.zeros(len(features)) if parameters[0]==0 else np.sign(parameters[0])*np.exp(np.log(abs(parameters[0]))+np.log(yscale)+features@parameters[1:])
        if not np.all(np.isfinite(result)):raise ValueError('Dự đoán tràn số.')
        return float(result[0]) if z.ndim==0 else result.reshape(z.shape)
    physical_residual=predict(x)-y
    if not np.all(np.isfinite(physical_residual)):raise ValueError('Phần dư thang y vượt miền float; đổi đơn vị.')
    sse=sum_products(physical_residual,physical_residual)
    report=Report('Bình phương tối thiểu phi tuyến — thang y',
                  [f'SSE thang y = {fmt(decimal_cell(sse),12)}; MSE = {fmt(decimal_cell(sse/len(x)),12)}.',
                   f'Bộ giải kết thúc: {solution.message}',
                   'Kết quả là cực tiểu cục bộ theo tham số đầu; cần xem phần dư và thử giá trị đầu khác nếu nhạy.'],
                  value=predict,scatter=(x,y))
    amplitude=0. if parameters[0]==0 else exp_parameter(np.log(abs(parameters[0]))+np.log(yscale),float(np.sign(parameters[0])))
    report.table('Tham số',['tên','giá trị'],zip(['A']+[f'b_{j+1}' for j in range(len(labels))],[amplitude,*parameters[1:]]))
    report.table('Đối chiếu',['x','y','dự đoán','phần dư'],zip(x,y,predict(x),physical_residual))
    if np.all(y==0):report.lines.append('Dữ liệu y=0: A=0 phù hợp nhưng các tham số số mũ không xác định duy nhất.')
    return report
