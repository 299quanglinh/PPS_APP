"""Finite differences, composite quadrature and rectangular double integrals."""
import math
import numpy as np
from numpy.polynomial import Polynomial as Poly
from scipy.integrate import newton_cotes
from .common import Report, points, spacing, fmt, function


def differentiate(x,y,target,degree=2,order=1):
    x,y,_=points(x,y,2)
    if not 1<=degree<len(x) or not 1<=order<=degree:
        raise ValueError('Cần 1 ≤ cấp đạo hàm ≤ bậc p < số mốc.')
    center=int(np.argmin(abs(x-target)))
    start=max(0,min(center-degree//2,len(x)-degree-1))
    xx,yy=x[start:start+degree+1],y[start:start+degree+1]
    scale=max(abs(xx-target))
    nodes=(xx-target)/scale
    weights=[]
    for i in range(len(xx)):
        poly=Poly.fromroots(np.delete(nodes,i))
        weights.append(poly.deriv(order)(0)/poly(nodes[i])/scale**order)
    result=float(np.dot(weights,yy))
    r=Report('Đạo hàm bằng công thức p+1 điểm',[
        f'p = {degree}; dùng {degree+1} mốc. Cấp đạo hàm = {order}.',
        f'f^({order})({fmt(target)}) ≈ {fmt(result)}',
        'Kết quả là tổng hệ số trọng số nhân tung độ; chưa có cận sai số nếu không biết đạo hàm cao hơn.'],value=result)
    r.table('Trọng số đạo hàm',['x','y','trọng số','đóng góp'],zip(xx,yy,weights,np.array(weights)*yy))
    r.presentation={'formulas':[('Công thức đạo hàm từ đa thức nội suy',[
        rf'f^{{({order})}}(\bar{{x}})\approx\sum_{{i=0}}^{{{degree}}}w_i y_i',
        rf'f^{{({order})}}({fmt(target)})\approx {fmt(result)}'])],'omit_prefixes':('f^(',)}
    if target<xx[0] or target>xx[-1]: r.lines.append('Điểm tính ngoài đoạn mốc đã chọn: ngoại suy.')
    return r


def quadrature_weights(n,method,order=4):
    if n<1: raise ValueError('Số khoảng chia phải dương.')
    if method=='Trung điểm': return np.ones(n),2
    if method=='Hình thang':
        w=np.ones(n+1); w[[0,-1]]=.5
        return w,2
    if method=='Simpson':
        if n%2: raise ValueError('Simpson cần số khoảng chia chẵn.')
        w=np.ones(n+1); w[1:-1:2]=4; w[2:-1:2]=2
        return w/3,4
    if method!='Newton-Cotes': raise ValueError('Phương pháp tích phân không hợp lệ.')
    if not 1<=order<=10 or n%order:
        raise ValueError('Bậc Newton-Cotes từ 1 đến 10; số khoảng chia phải chia hết cho bậc.')
    local,_=newton_cotes(order,equal=1)
    w=np.zeros(n+1)
    for i in range(0,n,order): w[i:i+order+1]+=local
    return w,order+2 if order%2==0 else order+1


def integrate(method='Simpson',expr='sin(x)',a=0.,b=1.,n=20,order=4,
              tolerance=None,x=None,y=None,bound=None):
    if not np.isfinite(a+b) or not a<b: raise ValueError('Cần a < b, hai cận hữu hạn.')
    if tolerance is not None and tolerance<=0: raise ValueError('Dung sai phải dương.')
    rows=[]; estimate=None; selection_note=''
    if bound is not None and (not np.isfinite(bound) or bound<0): raise ValueError('Cận đạo hàm M phải hữu hạn, không âm.')
    def calc(count):
        if not 1<=count<=1000000: raise ValueError('Giới hạn 1.000.000 khoảng chia.')
        weights,p=quadrature_weights(count,method,order)
        h=(b-a)/count
        nodes=a+(np.arange(count)+.5)*h if method=='Trung điểm' else np.linspace(a,b,count+1)
        vals=function(expr)(nodes)
        return h*np.dot(weights,vals),nodes,vals,weights,p
    if x is not None:
        if tolerance is not None: raise ValueError('Tự chia lưới cần biểu thức hàm, không dùng với bảng hữu hạn.')
        xx,yy,_=points(x,y,1 if method=='Trung điểm' else 2)
        mask=(xx>=a-1e-12)&(xx<=b+1e-12)
        xx,yy=xx[mask],yy[mask]
        n=len(xx) if method=='Trung điểm' else len(xx)-1
        weights,p=quadrature_weights(n,method,order); h=(b-a)/n
        expected=a+(np.arange(n)+.5)*h if method=='Trung điểm' else np.linspace(a,b,n+1)
        if not np.allclose(xx,expected,rtol=1e-8,atol=1e-12):
            raise ValueError('Bảng x phải khớp lưới [a,b]; với trung điểm, x là tâm của từng khoảng.')
        value=h*np.dot(weights,yy)
        if method!='Trung điểm' and n%2==0:
            try:
                coarse_weights,_=quadrature_weights(n//2,method,order)
                estimate=abs(value-2*h*np.dot(coarse_weights,yy[::2]))/(2**p-1)
            except ValueError: pass
    else:
        if tolerance is not None and bound is not None:
            _,p=quadrature_weights(order if method=='Newton-Cotes' else 2,method,order)
            if method=='Hình thang': constant=(b-a)**3*bound/12
            elif method=='Trung điểm': constant=(b-a)**3*bound/24
            elif method=='Simpson': constant=(b-a)**5*bound/180
            else:
                _,B=newton_cotes(order,equal=1)
                constant=(b-a)**(p+1)*abs(B)*bound/order
            needed=max(1,math.ceil((constant/tolerance)**(1/p)))
            block=2 if method=='Simpson' else order if method=='Newton-Cotes' else 1
            n=math.ceil(needed/block)*block
            selection_note=f'Chọn n = {n} từ cận đạo hàm M để cận lý thuyết không vượt dung sai.'
        value,xx,yy,weights,p=calc(n)
        if tolerance is not None and bound is None:
            for _ in range(18):
                new,xx,yy,weights,p=calc(n*2)
                estimate=abs(new-value)/(2**p-1)
                n*=2; rows.append([n,new,estimate]); value=new
                if estimate<=tolerance: break
            else: raise ValueError('Chưa đạt dung sai sau 18 lần chia đôi bước lưới.')
        else:
            refined,*_=calc(n*2)
            estimate=abs(refined-value)/(2**p-1)
    r=Report('Tích phân — '+method,[f'I ≈ {fmt(value)}; số khoảng chia n = {n}.'],value=float(value))
    if selection_note:r.lines.append(selection_note)
    if estimate is not None:
        r.lines += [f'Ước lượng Runge: {fmt(estimate,12)}.',
                    'Runge là ước lượng tiệm cận, không phải cận sai số được chứng minh.']
        if x is None and (tolerance is None or bound is not None): r.lines.append('Ước lượng Runge ở trên áp dụng cho kết quả lưới 2n; kết quả chính vẫn dùng n đã chọn.')
    if bound is not None:
        if bound<0: raise ValueError('Cận đạo hàm M phải không âm.')
        h=(b-a)/n
        if method=='Hình thang': err=(b-a)*h*h*bound/12
        elif method=='Trung điểm': err=(b-a)*h*h*bound/24
        elif method=='Simpson': err=(b-a)*h**4*bound/180
        else:
            _,B=newton_cotes(order,equal=1)
            err=n/order*abs(B)*h**(p+1)*bound
        r.lines.append(f'Cận lý thuyết ≤ {fmt(err,12)}, nếu M là cận đúng của |f^({p})| trên [a,b].')
    r.table('Các giá trị và trọng số',['x','f(x)','trọng số'],zip(xx,yy,weights))
    if rows: r.table('Chia đôi lưới',['n','tích phân','ước lượng Runge'],rows)
    r.curves=[(xx,yy,'f(x)')]
    formula={'Hình thang':r'I_h=\frac{h}{2}\left[y_0+2\sum_{i=1}^{n-1}y_i+y_n\right]',
             'Trung điểm':r'I_h=h\sum_{i=0}^{n-1}f\left(a+\left(i+\frac{1}{2}\right)h\right)',
             'Simpson':r'I_h=\frac{h}{3}\left[y_0+4\sum_{i=1,3,\ldots}^{n-1}y_i+2\sum_{i=2,4,\ldots}^{n-2}y_i+y_n\right]',
             'Newton-Cotes':r'I_h=h\sum_{i=0}^{n}w_i y_i'}[method]
    r.presentation={'formulas':[('Công thức tích phân',[formula]),('Kết quả',[
        rf'I=\int_{{{fmt(a)}}}^{{{fmt(b)}}}f(x)\,dx\approx {fmt(value)}',
        rf'h=\frac{{b-a}}{{n}},\qquad n={n}'])]}
    return r


def double_integral(expr,a,b,c,d,n=10,method='Simpson',tolerance=1e-7):
    if not a<b or not c<d or tolerance<=0: raise ValueError('Cần a<b, c<d và dung sai dương.')
    if method not in ('Simpson','Hình thang'): raise ValueError('Chọn Simpson hoặc hình thang.')
    f=function(expr,('x','y')); rows=[]
    def calc(k):
        if k>2000: raise ValueError('Lưới tích phân kép vượt 2000 × 2000; tăng dung sai.')
        w,p=quadrature_weights(k,method)
        xx=np.linspace(a,b,k+1); yy=np.linspace(c,d,k+1)
        vals=f(xx[:,None],yy[None,:])
        return float(w@vals@w*(b-a)*(d-c)/k**2),p
    old,p=calc(n)
    for _ in range(12):
        n*=2; value,p=calc(n); err=abs(value-old)/(2**p-1)
        rows.append([n,value,err])
        if err<=tolerance:
            report=Report('Tích phân kép trên miền chữ nhật',
                          [f'I ≈ {fmt(value)}; Runge ≈ {fmt(err,12)}.',
                           'Tự tăng số khoảng chia ở cả hai chiều. Runge không phải cận bảo đảm.'],value=value).table(
                               'Lịch sử chia lưới',['n mỗi chiều','I','Runge'],rows)
            report.presentation={'formulas':[('Tích phân kép',[
                rf'I=\int_{{{fmt(a)}}}^{{{fmt(b)}}}\int_{{{fmt(c)}}}^{{{fmt(d)}}}f(x,y)\,dy\,dx\approx {fmt(value)}'])]}
            return report
        old=value
    raise ValueError('Tích phân kép chưa hội tụ với dung sai đã nhập.')
