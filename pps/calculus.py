"""Finite differences, composite quadrature and rectangular double integrals."""
import math
from decimal import Decimal, localcontext
import numpy as np
from scipy.integrate import newton_cotes
from .common import Report, points, spacing, fmt, function
from .polynomial_model import nonnegative_integer, working_precision, decimal_cell


def _weighted_integral(weights, values, measure):
    """Cộng có đổi thang; không tràn tổng trước khi nhân bước rất nhỏ."""
    weights=np.asarray(weights,float);values=np.asarray(values,float)
    scale=float(np.max(abs(values)))
    if scale==0:return 0.
    nonzero=abs(values[values!=0])
    spread=max(0,math.log10(scale)-math.log10(float(np.min(nonzero))))
    with localcontext() as ctx:
        ctx.prec=max(80,int(spread)+50)
        if values.size<=512:
            total=sum((Decimal(str(float(w)))*Decimal(str(float(v)))
                       for w,v in zip(weights.flat,values.flat)),Decimal(0))
        else:
            normalized=math.fsum((weights*(values/scale)).flat)
            if spread>100 and abs(normalized)<1e-12:
                total=sum((Decimal(str(float(w)))*Decimal(str(float(v)))
                           for w,v in zip(weights.flat,values.flat)),Decimal(0))
            else:total=Decimal(str(normalized))*Decimal(str(scale))
        result=float(total*measure)
    if not np.isfinite(result):raise ValueError('Giá trị tích phân vượt miền float; đổi đơn vị dữ liệu.')
    return result


def _bound_constant(method, order, length, bound):
    """Trả C,p trong cận C/n^p, giữ các tích rất nhỏ/lớn bằng Decimal."""
    with localcontext() as ctx:
        ctx.prec=80
        length=Decimal(str(float(length)));bound=Decimal(str(float(bound)))
        if method=='Hình thang':return length**3*bound/12,2
        if method=='Trung điểm':return length**3*bound/24,2
        if method=='Simpson':return length**5*bound/180,4
        _,p=quadrature_weights(order,method,order)
        _,constant=newton_cotes(order,equal=1)
        return length**(p+1)*Decimal(str(abs(constant)))*bound/order,p


def differentiate(x,y,target,degree=2,order=1):
    x,y,_=points(x,y,2)
    degree=nonnegative_integer(degree,'Bậc đa thức')
    order=nonnegative_integer(order,'Cấp đạo hàm')
    if not np.isfinite(target): raise ValueError('Điểm tính đạo hàm phải hữu hạn.')
    if not 1<=degree<len(x) or not 1<=order<=degree:
        raise ValueError('Cần 1 ≤ cấp đạo hàm ≤ bậc p < số mốc.')
    center=int(np.argmin(abs(x-target)))
    start=max(0,min(center-degree//2,len(x)-degree-1))
    xx,yy=x[start:start+degree+1],y[start:start+degree+1]
    weights=[];contributions=[]
    with localcontext() as ctx:
        ctx.prec=working_precision(xx)
        nodes=[Decimal(str(float(v))) for v in xx]
        z=Decimal(str(float(target)))
        factorial=Decimal(math.factorial(order))
        for i,node in enumerate(nodes):
            # Chỉ sinh hệ số tới cấp cần lấy của L_i(target+u).
            # Không khai triển toàn bộ đa thức theo x bằng float.
            coefficients=[Decimal(1)]+[Decimal(0)]*order
            denominator=Decimal(1)
            for j,other in enumerate(nodes):
                if i==j:continue
                factor=z-other
                for k in range(order,0,-1):
                    coefficients[k]=factor*coefficients[k]+coefficients[k-1]
                coefficients[0]*=factor
                denominator*=node-other
            weight=coefficients[order]*factorial/denominator
            weights.append(weight)
            contributions.append(weight*Decimal(str(float(yy[i]))))
        result=float(sum(contributions,Decimal(0)))
    if not np.isfinite(result):
        raise ValueError('Giá trị đạo hàm vượt miền float; đổi đơn vị dữ liệu.')
    r=Report('Đạo hàm bằng công thức p+1 điểm',[
        f'p = {degree}; dùng {degree+1} mốc. Cấp đạo hàm = {order}.',
        f'f^({order})({fmt(target)}) ≈ {fmt(result)}',
        'Kết quả là tổng hệ số trọng số nhân tung độ; chưa có cận sai số nếu không biết đạo hàm cao hơn.'],value=result)
    r.table('Trọng số đạo hàm',['x','y','trọng số','đóng góp'],
            ([xv,yv,decimal_cell(w),decimal_cell(c)] for xv,yv,w,c in zip(xx,yy,weights,contributions)))
    r.presentation={'formulas':[('Công thức đạo hàm từ đa thức nội suy',[
        rf'f^{{({order})}}(\bar{{x}})\approx\sum_{{i=0}}^{{{degree}}}w_i y_i',
        rf'f^{{({order})}}({fmt(target)})\approx {fmt(result)}'])],'omit_prefixes':('f^(',)}
    if target<xx[0] or target>xx[-1]: r.lines.append('Điểm tính ngoài đoạn mốc đã chọn: ngoại suy.')
    return r


def quadrature_weights(n,method,order=4):
    n=nonnegative_integer(n,'Số khoảng chia')
    order=nonnegative_integer(order,'Bậc Newton–Cotes')
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
    if not np.all(np.isfinite([a,b,b-a])): raise ValueError('Hai cận và độ dài đoạn phải hữu hạn.')
    if tolerance is not None and (not np.isfinite(tolerance) or tolerance<=0): raise ValueError('Dung sai phải dương và hữu hạn.')
    rows=[]; checks=[]; estimate=None; selection_note=''
    if bound is not None and (not np.isfinite(bound) or bound<0): raise ValueError('Cận đạo hàm M phải hữu hạn, không âm.')
    if (x is None)!=(y is None):raise ValueError('Phải nhập đồng thời hai dãy x và y khi dùng bảng.')
    order=nonnegative_integer(order,'Bậc Newton–Cotes')
    quadrature_weights(order if method=='Newton-Cotes' else 2,method,order)
    original_a,original_b=a,b;direction=1 if b>=a else -1
    if a>b:a,b=b,a
    fun=function(expr) if x is None else None
    if a==b:
        if x is not None:points(x,y)
        return Report('Tích phân — '+method,['Hai cận bằng nhau: I = 0.'],value=0.)
    constant,p_bound=_bound_constant(method,order,b-a,bound) if bound is not None else (None,None)
    def calc(count):
        if not 1<=count<=1000000: raise ValueError('Giới hạn 1.000.000 khoảng chia.')
        weights,p=quadrature_weights(count,method,order)
        h=(b-a)/count
        nodes=a+(np.arange(count)+.5)*h if method=='Trung điểm' else np.linspace(a,b,count+1)
        if h==0 or np.any(np.diff(nodes)==0):raise ValueError('Lưới tích phân không phân biệt được các mốc bằng float; đổi gốc/đơn vị hoặc giảm n.')
        vals=fun(nodes)
        result=_weighted_integral(weights,vals,Decimal(str(b-a))/count)
        return result,nodes,vals,weights,p
    if x is not None:
        if tolerance is not None: raise ValueError('Tự chia lưới cần biểu thức hàm, không dùng với bảng hữu hạn.')
        xx,yy,_=points(x,y,1 if method=='Trung điểm' else 2)
        mask=(xx>=a-1e-12)&(xx<=b+1e-12)
        xx,yy=xx[mask],yy[mask]
        n=len(xx) if method=='Trung điểm' else len(xx)-1
        weights,p=quadrature_weights(n,method,order); h=(b-a)/n
        expected=a+(np.arange(n)+.5)*h if method=='Trung điểm' else np.linspace(a,b,n+1)
        grid_tolerance=min(abs(h)/8,max(abs(h)*1e-8,4*max(abs(np.spacing(a)),abs(np.spacing(b)))))
        if not np.allclose(xx,expected,rtol=0,atol=grid_tolerance):
            raise ValueError('Bảng x phải khớp lưới [a,b]; với trung điểm, x là tâm của từng khoảng.')
        value=_weighted_integral(weights,yy,Decimal(str(b-a))/n)
        if method!='Trung điểm' and n%2==0:
            try:
                coarse_weights,_=quadrature_weights(n//2,method,order)
                coarse=_weighted_integral(coarse_weights,yy[::2],Decimal(str(b-a))/(n//2))
                estimate=abs(value-coarse)/(2**p-1)
            except ValueError: pass
    else:
        if tolerance is not None and bound is not None:
            block=2 if method=='Simpson' else order if method=='Newton-Cotes' else 1
            limit=1000000//block
            requested=Decimal(str(float(tolerance)))
            if constant/Decimal(limit*block)**p_bound>requested:
                raise ValueError('Cận M yêu cầu hơn 1.000.000 khoảng; đổi phương pháp hoặc tăng dung sai.')
            lo,hi=1,limit
            while lo<hi:
                mid=(lo+hi)//2
                if constant/Decimal(mid*block)**p_bound<=requested:hi=mid
                else:lo=mid+1
            n=lo*block
            selection_note=f'Chọn n = {n} từ cận đạo hàm M để cận lý thuyết không vượt dung sai.'
        value,xx,yy,weights,p=calc(n)
        if tolerance is not None and bound is None:
            for _ in range(18):
                new,xx,yy,weights,p=calc(n*2)
                estimate=abs(new-value)/(2**p-1)
                n*=2; rows.append([n,new,estimate]); value=new
                if estimate<=tolerance:
                    block=2 if method=='Simpson' else order if method=='Newton-Cotes' else 1
                    probe_n=n+block if n+block<=1000000 else n-block
                    probe,px,py,pw,_=calc(probe_n)
                    mismatch=abs(probe-value)
                    checks.append([n,probe_n,probe,mismatch])
                    if mismatch<=tolerance:break
                    n=probe_n;value=probe;xx,yy,weights=px,py,pw
            else: raise ValueError('Chưa đạt dung sai sau 18 lần chia đôi bước lưới.')
        else:
            if n*2<=1000000:
                refined,*_=calc(n*2)
                estimate=abs(refined-value)/(2**p-1)
            else:
                selection_note += ' Không tính Runge vì lưới 2n vượt giới hạn tài nguyên 1.000.000 khoảng.'
    value*=direction
    if direction<0:
        rows=[[count,-v,e] for count,v,e in rows]
        selection_note+=' Đảo cận: đổi dấu tích phân; bảng trọng số dùng đoạn tăng dần.'
    r=Report('Tích phân — '+method,[f'I ≈ {fmt(value)}; số khoảng chia n = {n}.'],value=float(value))
    if selection_note:r.lines.append(selection_note)
    if estimate is not None:
        r.lines += [f'Ước lượng Runge: {fmt(estimate,12)}.',
                    'Runge là ước lượng tiệm cận, không phải cận sai số được chứng minh.']
        if x is None and (tolerance is None or bound is not None): r.lines.append('Ước lượng Runge ở trên áp dụng cho kết quả lưới 2n; kết quả chính vẫn dùng n đã chọn.')
    if bound is not None:
        err=constant/Decimal(n)**p_bound
        r.lines.append(f'Cận lý thuyết ≤ {fmt(decimal_cell(err),12)}, nếu M là cận đúng của |f^({p_bound})| trên [a,b].')
    r.table('Các giá trị và trọng số',['x','f(x)','trọng số'],zip(xx,yy,weights))
    if rows: r.table('Chia đôi lưới',['n','tích phân','ước lượng Runge'],rows)
    if checks:
        r.table('Kiểm tra thêm trên lưới không lồng nhau',['n','n kiểm tra','I kiểm tra','độ chênh'],
                ([n0,n1,direction*v,e] for n0,n1,v,e in checks))
        r.lines.append('Đã đối chiếu thêm lưới không lồng nhau để giảm bỏ sót dao động; phép kiểm tra này vẫn không phải cận sai số bảo đảm.')
    r.curves=[(xx,yy,'f(x)')]
    formula={'Hình thang':r'I_h=\frac{h}{2}\left[y_0+2\sum_{i=1}^{n-1}y_i+y_n\right]',
             'Trung điểm':r'I_h=h\sum_{i=0}^{n-1}f\left(a+\left(i+\frac{1}{2}\right)h\right)',
             'Simpson':r'I_h=\frac{h}{3}\left[y_0+4\sum_{i=1,3,\ldots}^{n-1}y_i+2\sum_{i=2,4,\ldots}^{n-2}y_i+y_n\right]',
             'Newton-Cotes':r'I_h=h\sum_{i=0}^{n}w_i y_i'}[method]
    r.presentation={'formulas':[('Công thức tích phân',[formula]),('Kết quả',[
        rf'I=\int_{{{fmt(original_a)}}}^{{{fmt(original_b)}}}f(x)\,dx\approx {fmt(value)}',
        rf'|h|=\frac{{|b-a|}}{{n}},\qquad n={n}'])]}
    return r


def double_integral(expr,a,b,c,d,n=10,method='Simpson',tolerance=1e-7):
    n=nonnegative_integer(n,'Số khoảng chia mỗi chiều')
    if not np.all(np.isfinite([a,b,c,d,b-a,d-c,tolerance])) or tolerance<=0:
        raise ValueError('Các cận/độ dài hữu hạn và dung sai phải dương hữu hạn.')
    if method not in ('Simpson','Hình thang'): raise ValueError('Chọn Simpson hoặc hình thang.')
    f=function(expr,('x','y')); rows=[];checks=[]
    original=(a,b,c,d)
    if a==b or c==d:return Report('Tích phân kép trên miền chữ nhật',['Miền có một chiều bằng 0: I = 0.'],value=0.)
    direction=(1 if a<b else -1)*(1 if c<d else -1)
    if a>b:a,b=b,a
    if c>d:c,d=d,c
    def calc(k):
        if k>2000: raise ValueError('Lưới tích phân kép vượt 2000 × 2000; tăng dung sai.')
        w,p=quadrature_weights(k,method)
        xx=np.linspace(a,b,k+1); yy=np.linspace(c,d,k+1)
        vals=f(xx[:,None],yy[None,:])
        if np.any(np.diff(xx)==0) or np.any(np.diff(yy)==0):raise ValueError('Các nút tích phân kép bị trùng do độ phân giải float; đổi gốc/đơn vị.')
        measure=Decimal(str(b-a))*Decimal(str(d-c))/k**2
        result=direction*_weighted_integral(w[:,None]*w[None,:],vals,measure)
        return result,p
    old,p=calc(n)
    for _ in range(12):
        n*=2; value,p=calc(n); err=abs(value-old)/(2**p-1)
        rows.append([n,value,err])
        if err<=tolerance:
            block=2 if method=='Simpson' else 1
            probe_n=n+block if n+block<=2000 else n-block
            probe,_=calc(probe_n);mismatch=abs(probe-value)
            checks.append([n,probe_n,probe,mismatch])
            if mismatch>tolerance:
                n=probe_n;old=probe
                continue
            report=Report('Tích phân kép trên miền chữ nhật',
                          [f'I ≈ {fmt(value)}; Runge ≈ {fmt(err,12)}.',
                           'Tự tăng số khoảng chia ở cả hai chiều. Runge không phải cận bảo đảm.'],value=value).table(
                               'Lịch sử chia lưới',['n mỗi chiều','I','Runge'],rows)
            report.table('Kiểm tra thêm trên lưới không lồng nhau',['n','n kiểm tra','I kiểm tra','độ chênh'],checks)
            report.presentation={'formulas':[('Tích phân kép',[
                rf'I=\int_{{{fmt(original[0])}}}^{{{fmt(original[1])}}}\int_{{{fmt(original[2])}}}^{{{fmt(original[3])}}}f(x,y)\,dy\,dx\approx {fmt(value)}'])]}
            return report
        old=value
    raise ValueError('Tích phân kép chưa hội tụ với dung sai đã nhập.')
