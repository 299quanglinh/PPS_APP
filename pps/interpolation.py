"""Polynomial, central and inverse interpolation; splines and regression."""
import math
import numpy as np
from numpy.polynomial import Polynomial as Poly
from scipy.interpolate import CubicSpline, PPoly
from .common import Report, points, spacing, difference, fmt, polynomial_text, function

METHODS = ['Lagrange', 'Newton bất kỳ', 'Newton tiến', 'Newton lùi',
           'Gauss I', 'Gauss II', 'Stirling', 'Bessel']


def interpolate(x, y, method='Lagrange', target=None, derivative=0):
    x, y, dropped = points(x, y)
    n = len(x)
    if n > 100:
        raise ValueError('Chọn tối đa 100 mốc cho đa thức; dùng spline cho bảng lớn.')
    if method not in METHODS+['Newton bất kỳ lùi']: raise ValueError('Phương pháp không hợp lệ.')
    central = method in ('Gauss I', 'Gauss II', 'Stirling', 'Bessel')
    if central and method != 'Bessel' and n % 2 == 0:
        raise ValueError('Chọn số mốc lẻ, đối xứng quanh x_0 cho Gauss/Stirling.')
    if method == 'Bessel' and n % 2:
        raise ValueError('Bessel cần số mốc chẵn, ít nhất 2 mốc.')
    normalized = method not in ('Lagrange', 'Newton bất kỳ', 'Newton bất kỳ lùi')
    h = spacing(x) if normalized else 1.0
    center = (n-1)//2 if central else (n-1 if method == 'Newton lùi' else 0)
    origin = x[center] if normalized else 0.0
    a = difference(x, y, method.startswith('Newton bất kỳ'))
    q, product = Poly([0.]), Poly([1.])
    steps, products = [], []
    if method == 'Lagrange':
        for i in range(n):
            basis = Poly.fromroots(np.delete(x, i))
            c = y[i]/basis(x[i])
            q += c*basis
            steps.append([i, c, polynomial_text(basis)])
            products.append(basis.coef)
    elif method in ('Newton bất kỳ', 'Newton bất kỳ lùi', 'Newton tiến', 'Newton lùi', 'Gauss I', 'Gauss II'):
        if method == 'Gauss I':
            order = [0] + [v for r in range(1, center+1) for v in (r, -r)]
        elif method == 'Gauss II':
            order = [0] + [v for r in range(1, center+1) for v in (-r, r)]
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
    result = Report(method, value=p, scatter=(x,y))
    if dropped: result.lines.append(f'Đã bỏ {dropped} mốc x trùng, giữ tung độ đầu tiên.')
    if normalized:
        result.lines += [f'x_0 = {fmt(origin)}; h = {fmt(h)}; t = (x - ({fmt(origin)}))/{fmt(h)}',
                         'P(x_0+th) = '+polynomial_text(q, 't')]
    result.lines += ['P(x) = '+polynomial_text(p),
                     'Hệ số hiển thị làm tròn; phép tính dùng hệ số đầy đủ.']
    err = np.max(np.abs(q((x-origin)/h) - y)) if normalized else np.max(np.abs(p(x)-y))
    result.lines.append('Sai lệch lớn nhất tại các mốc: '+fmt(err, 12))
    if err > 1e-6 * max(1, np.max(np.abs(y))):
        result.lines.append('CẢNH BÁO: đa thức bậc cao/khoảng mốc gây mất chính xác; hãy giảm số mốc hoặc dùng spline.')
    if target is not None:
        if derivative < 0: raise ValueError('Cấp đạo hàm không âm.')
        val = q.deriv(derivative)((target-origin)/h)/h**derivative if normalized else p.deriv(derivative)(target)
        result.lines.append(f'P^({derivative})({fmt(target)}) = {fmt(val)}')
        if target < x[0] or target > x[-1]: result.lines.append('Điểm tính ngoài tập mốc: đang ngoại suy.')
    result.table('Bảng tỷ sai phân' if method.startswith('Newton bất kỳ') else 'Bảng sai phân tiến',
                 ['x','y']+[f'bậc {k}' for k in range(1,n)],
                 [[x[j]]+[None if np.isnan(v) else v for v in a[j]] for j in range(n)])
    if products:
        if method == 'Lagrange':
            table_name = 'Bảng tích: tử số cơ sở Lagrange — hệ số theo bậc giảm dần'
        elif method.startswith('Newton bất kỳ'):
            table_name = 'Bảng tích: đa thức cơ sở Newton — hệ số theo bậc giảm dần'
        else:
            table_name = 'Bảng tích: chưa chia giai thừa — hệ số theo bậc giảm dần'
        result.table(table_name, ['k']+[f'bậc {k}' for k in range(n-1,-1,-1)],
                     [[i]+list(reversed(list(c)+[0]*(n-len(c)))) for i,c in enumerate(products)])
    grid = np.linspace(x[0],x[-1],400)
    result.curves = [(grid, q((grid-origin)/h) if normalized else p(grid), method)]
    from .math_view import interpolation_presentation
    interpolation_presentation(result,method,x,y,a,q,p,origin,h,center,normalized,target,derivative)
    if normalized:
        from .math_view import product_formula
        product_formula(result,method,products,center)
    return result


def inverse(x, y, target, backward=False, tolerance=1e-10):
    x,y,_ = points(x,y,2)
    h = spacing(x)
    if not (np.all(np.diff(y)>0) or np.all(np.diff(y)<0)):
        raise ValueError('Chọn một đoạn dữ liệu đơn điệu nghiêm ngặt cho nội suy ngược.')
    if not min(y) <= target <= max(y): raise ValueError('Giá trị y cần tìm nằm ngoài đoạn đã chọn.')
    if tolerance <= 0: raise ValueError('Dung sai phải dương.')
    a = difference(x,y)
    j = len(x)-1 if backward else 0
    d1 = a[-2,1] if backward else a[0,1]
    t = (target-y[j])/d1
    rows = []
    for k in range(1,1001):
        prod, total = t, 0.0
        for r in range(2,len(x)):
            prod *= (t+r-1 if backward else t-r+1)/r
            total += a[len(x)-r-1 if backward else 0,r]*prod
        next_t = (target-y[j]-total)/d1
        rows.append([k,t,total,next_t,abs(next_t-t)])
        if not np.isfinite(next_t) or abs(next_t)>1e8: break
        if abs(next_t-t)*h < tolerance:
            result = x[j]+h*next_t
            residual = interpolate(x,y,'Newton bất kỳ').value(result)-target
            if abs(residual)>max(tolerance*10,1e-8)*max(1,abs(target)): break
            report=Report('Nội suy ngược — lặp Newton '+('lùi' if backward else 'tiến'),
                          [f'x ≈ {fmt(result)}; sai lệch P(x)-y = {fmt(residual,12)}'],value=result).table(
                              'Lịch sử lặp',['k','t trước','tổng bậc cao','t sau','độ lệch'],rows)
            report.presentation={'formulas':[('Nghiệm nội suy ngược',[rf'P(x)={fmt(target)},\qquad x\approx {fmt(result)}'])]}
            return report
        t = next_t
    raise ValueError('Phép lặp không hội tụ. Chọn ít mốc hơn gần giá trị y cần tìm hoặc đổi hướng tiến/lùi.')


def select_points(x,y,target,count=5,reverse=False):
    x,y,_ = points(x,y,2)
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
        spacing(x)
        if not 1<=count<=len(x): raise ValueError('Số mốc phải nằm trong độ dài bảng.')
        if not x[0]<=target<=x[-1]: raise ValueError('Điểm chọn mốc phải nằm trong bảng.')
        # Contiguous window minimizing distance, with leftmost tie break.
        start = min(range(len(x)-count+1),key=lambda j: np.sum(abs(x[j:j+count]-target)))
        xx,yy=x[start:start+count],y[start:start+count]
        r.table('Mốc đã chọn',['x','y'],zip(xx,yy))
        r.lines += ['Dãy x: '+' '.join(fmt(v) for v in xx),'Dãy y: '+' '.join(fmt(v) for v in yy)]
        r.value=(xx,yy)
    return r


def chebyshev(degree,a,b):
    if not 0<=degree<=1000 or not a<b: raise ValueError('Bậc từ 0 đến 1000; a < b.')
    x=np.sort((a+b)/2+(b-a)/2*np.cos((2*np.arange(degree+1)+1)*np.pi/(2*(degree+1))))
    report=Report('Mốc Chebyshev',value=x).table('Mốc tối ưu',['i','x'],enumerate(x))
    report.presentation={'formulas':[('Công thức mốc Chebyshev',[r'x_i=\frac{a+b}{2}+\frac{b-a}{2}\cos\frac{(2i+1)\pi}{2(n+1)}'])]}
    return report


def horner(coefficients,target,order=0):
    from .math_view import poly_latex,split_equation
    if order<0 or order>100: raise ValueError('Cấp đạo hàm từ 0 đến 100.')
    p=Poly(np.asarray(coefficients)[::-1])
    r=Report('Horner', ['Hệ số nhập theo bậc giảm dần.','P(x) = '+polynomial_text(p)],value=p(target))
    r.presentation={'formulas':[('Đa thức',split_equation('P(x)',poly_latex(p)))],'omit_prefixes':('P(x) =','P^(')}
    for k in range(order+1):
        c=p.deriv(k).coef[::-1]; row=[]; value=0.
        for v in c: value=value*target+v; row.append(value)
        r.table(f'Horner cho đạo hàm cấp {k}', ['hệ số','giá trị tích lũy'],zip(c,row))
        r.lines.append(f'P^({k})({fmt(target)}) = {fmt(value)}')
        label='P' if k==0 else rf'P^{{({k})}}'
        r.presentation['formulas'].append((f'Giá trị — cấp {k}',[rf'{label}({fmt(target)})={fmt(value)}']))
    return r


def spline(x,y,kind='Cubic tự nhiên',left=0.,right=0.,target=None):
    from .math_view import poly_latex,split_equation
    x,y,dropped=points(x,y,2); h=np.diff(x); slopes=np.diff(y)/h
    if kind=='Tuyến tính': pp=PPoly(np.array([slopes,y[:-1]]),x)
    elif kind.startswith('Bậc hai'):
        d=np.zeros(len(x))
        if kind.endswith('phải'):
            d[-1]=right
            for i in range(len(h)-1,-1,-1): d[i]=2*slopes[i]-d[i+1]
        else:
            d[0]=left
            for i in range(len(h)): d[i+1]=2*slopes[i]-d[i]
        pp=PPoly(np.array([(slopes-d[:-1])/h,d[:-1],y[:-1]]),x)
    else:
        bc='natural' if kind=='Cubic tự nhiên' else ((1,left),(1,right)) if kind=='Cubic kẹp' else ((2,left),(2,right))
        pp=CubicSpline(x,y,bc_type=bc)
    r=Report('Spline — '+kind, ['Mỗi đoạn dùng biến (x - mốc trái).'],scatter=(x,y),value=pp)
    r.presentation={'formulas':[],'omit_prefixes':('[','S(')}
    for i in range(len(h)):
        r.lines.append(f'[{fmt(x[i])}, {fmt(x[i+1])}]: '+polynomial_text(Poly(pp.c[:,i][::-1]),f'(x-({fmt(x[i])}))'))
        r.presentation['formulas'].append((f'Đoạn [{fmt(x[i])}, {fmt(x[i+1])}]',split_equation(rf'S_{{{i}}}(x)',poly_latex(Poly(pp.c[:,i][::-1]),rf'(x-({fmt(x[i])}))'))))
    if target is not None:
        r.lines.append(f'S({fmt(target)}) = {fmt(pp(target))}')
        r.presentation['formulas'].append(('Thay giá trị',[rf'S({fmt(target)})={fmt(pp(target))}']))
        if not x[0]<=target<=x[-1]: r.lines.append('Điểm tính ngoài bảng: ngoại suy bằng đoạn đầu/cuối.')
    grid=np.linspace(x[0],x[-1],400); r.curves=[(grid,pp(grid),kind)]
    return r


def least_squares(x,y,kind='Tuyến tính theo cơ sở',basis=('1','x')):
    from .math_view import symbolic_latex,split_equation
    x,y,_=points(x,y,2); yshift=xshift=0.; sign=1.
    if kind!='Tuyến tính theo cơ sở':
        if min(y)<=0<=max(y):
            up=abs(min(y))+1; down=-(abs(max(y))+1)
            yshift=down if abs(down)<abs(up) else up
        sign=-1. if max(y+yshift)<0 else 1.
        transformed=np.log(abs(y+yshift))
    else: transformed=y
    if kind=='Lũy thừa':
        xshift=abs(min(x))+1 if min(x)<=0 else 0
        funcs=[lambda v: np.ones_like(v),lambda v: np.log(v+xshift)]
        labels=['1',f'ln(x+{fmt(xshift)})']
    else:
        labels=list(basis)
        if kind=='Hàm mũ': labels=['1']+[v for v in labels if v.strip()!='1']
        funcs=[function(v) for v in labels]
    theta=np.column_stack([f(x) for f in funcs])
    if len(x)<theta.shape[1]: raise ValueError('Số điểm ít hơn số hàm cơ sở.')
    coeffs,_,rank,_=np.linalg.lstsq(theta,transformed,rcond=None)
    if rank<theta.shape[1]: raise ValueError('Các hàm cơ sở phụ thuộc tuyến tính trên tập mốc.')
    def pred(v):
        val=sum(c*f(np.asarray(v)) for c,f in zip(coeffs,funcs))
        return val if kind=='Tuyến tính theo cơ sở' else sign*np.exp(val)-yshift
    text=' + '.join(f'({fmt(c)})*({s})' for c,s in zip(coeffs,labels))
    if kind!='Tuyến tính theo cơ sở': text=f'({fmt(sign)})*exp({text}) - ({fmt(yshift)})'
    r=Report('Bình phương tối thiểu — '+kind,['f(x) = '+text,
        'MSE trên dữ liệu gốc = '+fmt(np.mean((pred(x)-y)**2)),
        f'Dịch x = {fmt(xshift)}; dịch y = {fmt(yshift)}.',
        'Giải bằng SVD để tránh bình phương số điều kiện của ma trận.'],value=pred,scatter=(x,y))
    r.table('Hệ số',['cơ sở','hệ số'],zip(labels,coeffs))
    r.presentation={'formulas':[('Hàm xấp xỉ',split_equation('f(x)',symbolic_latex(text))),
                                ('Bình phương tối thiểu',[r'\Theta^T\Theta\,a=\Theta^T y',r'\mathrm{MSE}=\frac{1}{N}\sum_{i=0}^{N-1}(f(x_i)-y_i)^2'])],
                    'omit_prefixes':('f(x) =',)}
    r.table('Ma trận thiết kế Θ',labels,theta)
    r.table('Ma trận ΘᵀΘ',[str(i) for i in range(len(coeffs))],theta.T@theta)
    r.table('Đối chiếu',['x','y','dự đoán','sai lệch'],zip(x,y,pred(x),pred(x)-y))
    grid=np.linspace(x[0],x[-1],400); r.curves=[(grid,pred(grid),kind)]
    return r
