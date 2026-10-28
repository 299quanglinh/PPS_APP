"""Euler, Runge-Kutta, Adams, finite-difference BVP and Sturm-Liouville."""
import numpy as np
import sympy as sp
from decimal import Decimal, localcontext
from scipy.optimize import root
from scipy.linalg import eigh_tridiagonal, solve_banded
from .common import Report, expression, function, fmt
from .polynomial_model import nonnegative_integer

RK = {
    'Euler hiện': ([[0]], [1], [0]),
    'RK2': ([[0,0],[1,0]], [.5,.5], [0,1]),
    'RK3': ([[0,0,0],[.5,0,0],[-1,2,0]], [1/6,2/3,1/6], [0,.5,1]),
    'RK4': ([[0,0,0,0],[.5,0,0,0],[0,.5,0,0],[0,0,1,0]], [1/6,1/3,1/3,1/6], [0,.5,.5,1])
}
AB={2:np.array([3,-1])/2,3:np.array([23,-16,5])/12,
    4:np.array([55,-59,37,-9])/24,5:np.array([1901,-2774,2616,-1274,251])/720}
AM={2:np.array([5,8,-1])/12,3:np.array([9,19,-5,1])/24,
    4:np.array([251,646,-264,106,-19])/720,
    5:np.array([475,1427,-798,482,-173,27])/1440}


def custom_rk2(parameter,value):
    if parameter in ('alpha2','beta11'): c=value
    elif parameter=='r2':
        if value==0: raise ValueError('r2 phải khác 0.')
        c=1/(2*value)
    elif parameter=='r1':
        if value==1: raise ValueError('r1 phải khác 1.')
        c=1/(2*(1-value))
    else: raise ValueError('Tham số RK2: alpha2, beta11, r1 hoặc r2.')
    if c==0: raise ValueError('alpha2 phải khác 0.')
    return tableau([[0,0],[c,0]],[1-1/(2*c),1/(2*c)],[0,c])


def custom_rk3(c2,c3):
    if not np.all(np.isfinite([c2,c3])) or c2==0:
        raise ValueError('alpha2 phải khác 0; các tham số RK3 phải hữu hạn.')
    if c2==c3 or c3==0:
        if not np.isclose(c2,2/3,rtol=0,atol=2e-15):
            raise ValueError('Khi alpha3=0 hoặc alpha3=alpha2, RK3 cần alpha2=2/3.')
        # Hai họ suy biến của công thức tổng quát vẫn có bộ hệ số bậc ba.
        if c3==0:return tableau([[0,0,0],[2/3,0,0],[-1,1,0]],[0,.75,.25],[0,2/3,0])
        return tableau([[0,0,0],[2/3,0,0],[0,2/3,0]],[.25,.375,.375],[0,2/3,2/3])
    b2=(3*c3-2)/(6*c2*(c3-c2)); b3=(2-3*c2)/(6*c3*(c3-c2))
    if b3==0: raise ValueError('Bộ tham số này làm r3 = 0; đổi alpha2.')
    a32=1/(6*b3*c2)
    return tableau([[0,0,0],[c2,0,0],[c3-a32,a32,0]],[1-b2-b3,b2,b3],[0,c2,c3])


def grid(a,b,h,limit=20000,allow_reverse=False):
    if not np.all(np.isfinite([a,b,h])) or h==0 or a==b or ((b>a)!=(h>0)) or (not allow_reverse and a>b):
        raise ValueError('Bước h phải cùng chiều a → b; bài toán biên cần a < b và h > 0. Các số phải hữu hạn.')
    ratio=(b-a)/h
    if not np.isfinite(ratio) or ratio>limit+1e-6:
        raise ValueError(f'Lưới vượt giới hạn tài nguyên {limit} khoảng; tăng |h|.')
    n=round(ratio)
    if not 1<=n<=limit: raise ValueError(f'Số khoảng chia phải từ 1 đến {limit}.')
    if not np.isclose(ratio,n,rtol=1e-9,atol=1e-9):
        raise ValueError('Chọn h sao cho (b-a)/h là số nguyên để lưới kết thúc đúng tại b.')
    nodes=np.linspace(a,b,n+1)
    if np.any(np.diff(nodes)==0):
        raise ValueError('Bước thời gian nhỏ hơn độ phân giải float tại gốc này; đổi gốc/đơn vị thời gian.')
    return nodes


def tableau(a,b,c):
    a,b,c=np.asarray(a,float),np.asarray(b,float),np.asarray(c,float)
    if b.ndim!=1: raise ValueError('Trọng số b của bảng RK phải là dãy một chiều.')
    s=len(b)
    if not 1<=s<=10 or a.shape!=(s,s) or c.shape!=(s,): raise ValueError('Bảng RK cần A s×s, b và c có s phần tử (s≤10).')
    if not all(np.all(np.isfinite(v)) for v in (a,b,c)): raise ValueError('Hệ số RK phải hữu hạn.')
    if np.any(np.triu(a)!=0): raise ValueError('RK hiện cần A tam giác dưới nghiêm ngặt.')
    if not np.allclose(a.sum(axis=1),c) or not np.isclose(sum(b),1):
        raise ValueError('Bảng RK phải thỏa c = tổng hàng A và tổng b = 1.')
    return a,b,c


def rkstep(f,t,y,h,tab):
    a,b,c=tab; stages=[]
    for i in range(len(b)):
        dy=sum((a[i][j]*stages[j] for j in range(i)),np.zeros_like(y))
        stages.append(f(t+c[i]*h,y+h*dy))
    return y+h*sum((b[j]*stages[j] for j in range(len(b))),np.zeros_like(y))


def ode(expressions,y0,a,b,h,method='RK4',tolerance=1e-10,custom=None):
    y0=np.asarray(y0,float)
    if y0.ndim!=1 or not 1<=len(y0)<=3 or len(expressions)!=len(y0): raise ValueError('Nhập 1–3 phương trình và cùng số giá trị ban đầu.')
    if not np.isfinite(tolerance) or tolerance<=0 or not np.all(np.isfinite(y0)): raise ValueError('Dung sai dương hữu hạn, giá trị đầu hữu hạn.')
    if a==b and np.all(np.isfinite([a,b,h])) and h!=0:ts=np.array([a],float)
    else:ts=grid(a,b,h,allow_reverse=True)
    if len(ts)>1:h=(b-a)/(len(ts)-1)
    # Chỉ cho phép các biến trạng thái thực sự có trong hệ đã nhập.
    names=['x','y','z'][:len(y0)]
    variables=('t',*names)
    funcs=[function(s,variables) for s in expressions]
    def f(t,state):
        return np.array([float(fun(t,*state)) for fun in funcs])
    vals=[y0.copy()]; stats=[]
    if method=='RK tùy chỉnh':
        if custom is None: raise ValueError('Cần nhập bảng hệ số RK.')
        tab=tableau(*custom)
    elif method in RK: tab=RK[method]
    elif method not in ('Euler ẩn','Hình thang') and method not in [f'Adams {s}' for s in AB]:
        raise ValueError('Phương pháp ODE không hợp lệ.')
    implicit_method=method in ('Euler ẩn','Hình thang') or method.startswith('Adams')
    if implicit_method:
        parsed=[expression(s,variables) for s in expressions]
        symbols=parsed[0][1]
        raw_jac=sp.lambdify(symbols,sp.Matrix([item[0] for item in parsed]).jacobian(symbols[1:]),modules='numpy')
        def jac(t,state):
            with np.errstate(all='ignore'):value=np.asarray(raw_jac(t,*state))
            if np.iscomplexobj(value) or not np.all(np.isfinite(value)):
                raise ValueError('Jacobian của bước ẩn không hữu hạn; kiểm tra miền hoặc giảm bước.')
            return value.astype(float)
    def implicit(residual,guess,jacobian):
        sol=root(residual,guess,jac=jacobian,tol=tolerance)
        err=float(np.max(abs(residual(sol.x))))
        if not np.all(np.isfinite(sol.x)) or not np.isfinite(err) or err>max(tolerance,32*np.finfo(float).eps)*(1+np.max(abs(sol.x))):
            raise ValueError('Bước ẩn không hội tụ; hãy giảm h hoặc kiểm tra hàm.')
        singular=np.linalg.svd(jacobian(sol.x),compute_uv=False)
        if singular[-1]<=64*np.finfo(float).eps*max(1.,singular[0]):
            raise ValueError('Jacobian bước ẩn suy biến hoặc gần suy biến; không xác nhận nghiệm, hãy giảm |h|.')
        return sol.x,err
    for i,t in enumerate(ts[:-1]):
        y=vals[-1]; tn=ts[i+1]; err=0.
        if method in RK or method=='RK tùy chỉnh': new=rkstep(f,t,y,h,tab)
        elif method in ('Euler ẩn','Hình thang'):
            initial=f(t,y); theta=1. if method=='Euler ẩn' else .5
            new,err=implicit(lambda v:v-y-h*((1-theta)*initial+theta*f(tn,v)),y+h*initial,
                             lambda v:np.eye(len(y))-h*theta*jac(tn,v))
        else:
            s=int(method.split()[-1])
            if i<s-1: new=rkstep(f,t,y,h,RK['RK4'])
            else:
                history=[f(ts[i-j],vals[i-j]) for j in range(s)]
                predicted=y+h*sum((c*v for c,v in zip(AB[s],history)),np.zeros_like(y))
                tail=sum((c*v for c,v in zip(AM[s][1:],history)),np.zeros_like(y))
                new,err=implicit(lambda v:v-y-h*(AM[s][0]*f(tn,v)+tail),predicted,
                                 lambda v:np.eye(len(y))-h*AM[s][0]*jac(tn,v))
        if not np.all(np.isfinite(new)): raise ValueError('Nghiệm tràn số; giảm bước h hoặc kiểm tra phương trình.')
        vals.append(new.copy()); stats.append([i+1,tn,err])
    vals=np.array(vals)
    names=['x','y','z'][:len(y0)]
    r=Report('Bài toán Cauchy — '+method,
        ['Biến độc lập t; các thành phần nghiệm x(t), y(t), z(t).',
         'Các bước ẩn giải hệ phương trình và kiểm tra phần dư.',
         'Lưới cố định: dung sai bước ẩn không phải cận sai số toàn cục hoặc bảo đảm ổn định; cần xét thêm h và h/2.',
         'Adams dùng RK4 khởi tạo, AB dự báo và AM hiệu chỉnh đến hội tụ.' if method.startswith('Adams') else ''],value=(ts,vals))
    if method.startswith('Adams'):
        s=int(method.split()[-1])
        r.lines.append(f'Adams {s}: AB dùng {s} giá trị cũ, AM dùng thêm giá trị mới; bậc toàn cục còn phụ thuộc RK4 khởi tạo.')
    r.table('Nghiệm trên lưới',['t']+names,np.column_stack([ts,vals]))
    r.table('Phần dư bước ẩn',['bước','t','phần dư'],stats)
    if method in RK or method=='RK tùy chỉnh':
        aa,bb,cc=map(np.asarray,tab)
        r.table('Bảng hệ số RK — mỗi hàng gồm c_i và A_ij',
                ['c_i']+[f'A_{j+1}' for j in range(len(bb))],np.column_stack([cc,aa]))
        r.table('Trọng số RK',['i','b_i'],enumerate(bb,1))
    r.curves=[(ts,vals[:,j],names[j]+'(t)') for j in range(len(y0))]
    from .math_view import symbolic_latex
    r.presentation={'formulas':[('Phương trình và điều kiện đầu',[
        rf'\frac{{d{names[j]}}}{{dt}}='+symbolic_latex(expressions[j],('t','x','y','z'))+rf',\qquad {names[j]}({fmt(a)})={fmt(y0[j])}' for j in range(len(y0))])]}
    return r


def bvp(p,q,f,a,b,h,left=(0,1,0),right=(0,1,0)):
    xx=grid(a,b,h,10000); n=len(xx); h=(b-a)/(n-1); hh=h*h
    p,q,f=function(p),function(q),function(f)
    low=np.zeros(n); diag=np.zeros(n); upper=np.zeros(n); rhs=np.zeros(n)
    for i in range(1,n-1):
        low[i]=p(xx[i]-h/2); upper[i]=p(xx[i]+h/2)
        diag[i]=-(low[i]+upper[i]+hh*q(xx[i])); rhs[i]=-hh*f(xx[i])
    # Preserve upstream half-cell Robin discretization.
    for index,bc,direction in ((0,left,1),(-1,right,-1)):
        if len(bc)!=3 or not np.all(np.isfinite(bc)): raise ValueError('Biên gồm alpha beta gamma hữu hạn.')
        alpha,beta,gamma=bc
        if alpha==0:
            if beta==0: raise ValueError('Alpha và beta không thể đồng thời bằng 0.')
            diag[index]=1; rhs[index]=gamma/beta
        else:
            ph=float(p(xx[index]+direction*h/2))
            boundary_p=float(p(xx[index]))
            diag[index]=-ph-hh*q(xx[index])/2+direction*h*beta/alpha*boundary_p
            rhs[index]=-hh*f(xx[index])/2+direction*h*gamma/alpha*boundary_p
            if index==0: upper[index]=ph
            else: low[index]=ph
    if left[1]==right[1]==0 and np.all(q(xx)==0):
        raise ValueError('Neumann thuần với q=0 không xác định nghiệm duy nhất; cần thêm điều kiện chuẩn hóa hoặc một biên giá trị.')
    band=np.zeros((3,n)); band[0,1:]=upper[:-1]; band[1]=diag; band[2,:-1]=low[1:]
    try: yy=solve_banded((1,1),band,rhs)
    except np.linalg.LinAlgError as exc: raise ValueError('Hệ biên suy biến; kiểm tra điều kiện biên (Neumann thuần có thể không duy nhất).') from exc
    residual=diag*yy+np.r_[0,low[1:]*yy[:-1]]+np.r_[upper[:-1]*yy[1:],0]-rhs
    if not np.all(np.isfinite(yy)) or not np.all(np.isfinite(residual)):
        raise ValueError('Nghiệm hệ biên tràn số; kiểm tra hệ số và điều kiện biên.')
    r=Report('Bài toán biên — sai phân hữu hạn',[
        "[p(x)u']' - q(x)u = -f(x)","Biên: alpha*u' + beta*u = gamma.",
        'Sai phân tại nút trong: A*u[i-1] + B*u[i] + C*u[i+1] = F.',
         'Biên Robin dùng công thức nửa ô: p tại biên nhân đạo hàm biên; p tại nửa bước nhân thông lượng trong ô.',
         'Phần dư nhỏ xác nhận hệ sai phân; sai số so với nghiệm liên tục cần xét khi giảm h.',
        f'Phần dư hệ tuyến tính lớn nhất: {fmt(max(abs(residual)),12)}',
        f'Min u = {fmt(min(yy))}; max u = {fmt(max(yy))}'],value=(xx,yy))
    r.table('Bảng hệ số và nghiệm',['x','A','B','C','F','u'],zip(xx,low,diag,upper,rhs,yy))
    r.curves=[(xx,yy,'u(x)')]
    r.presentation={'formulas':[('Bài toán biên',[
        r'\frac{d}{dx}\left(p(x)\frac{du}{dx}\right)-q(x)u=-f(x)',
        rf'({fmt(left[0])})u^\prime(a)+({fmt(left[1])})u(a)={fmt(left[2])}',
        rf'({fmt(right[0])})u^\prime(b)+({fmt(right[1])})u(b)={fmt(right[2])}']),
        ('Hệ sai phân',[r'A_i u_{i-1}+B_i u_i+C_i u_{i+1}=F_i'])],
        'omit_prefixes':('[p(x)',)}
    return r


def eigen(p,q,r,a,b,h,count=5):
    xx=grid(a,b,h,10000); inside=xx[1:-1]
    h=(b-a)/(len(xx)-1)
    count=nonnegative_integer(count,'Số trị riêng')
    if not 1<=count<=len(inside): raise ValueError('Số trị riêng phải từ 1 đến số nút trong.')
    p,q,r=function(p),function(q),function(r)
    weights=r(inside); faces=p(xx[:-1]/2+xx[1:]/2)
    if np.any(weights<=0) or np.any(faces<=0): raise ValueError('Bài toán Sturm–Liouville này yêu cầu p > 0 và r > 0 trên lưới.')
    q_values=q(inside)
    # Lập ma trận đối xứng và đổi thang trước LAPACK. Nếu chỉ đổi công
    # thức sqrt(r_i*r_j), LAPACK vẫn có thể coi phần tử rất nhỏ là 0.
    with localcontext() as ctx:
        ctx.prec=80
        dh=Decimal(str(float(h)))
        dw=[Decimal(str(float(v))) for v in weights]
        dp=[Decimal(str(float(v))) for v in faces]
        raw_diag=[(dp[i]+dp[i+1])/(dh*dh*dw[i])+Decimal(str(float(v)))/dw[i] for i,v in enumerate(q_values)]
        raw_off=[-dp[i+1]/(dh*dh*(dw[i]*dw[i+1]).sqrt()) for i in range(len(dw)-1)]
        matrix_scale=max((abs(v) for v in raw_diag+raw_off),default=Decimal(1)) or Decimal(1)
        scaled_diag=np.array([float(v/matrix_scale) for v in raw_diag])
        scaled_off=np.array([float(v/matrix_scale) for v in raw_off])
        eig,vec=eigh_tridiagonal(scaled_diag,scaled_off,select='i',select_range=(0,count-1))
        eig=np.array([float(Decimal(str(float(v)))*matrix_scale) for v in eig])
    if not np.all(np.isfinite(eig)):raise ValueError('Trị riêng vượt miền float; đổi đơn vị bài toán.')
    vec=vec/np.sqrt(weights[:,None]); vec=np.pad(vec,((1,1),(0,0)))
    vec/=np.max(abs(vec),axis=0)
    result=Report('Trị riêng Sturm–Liouville',[
        "-(p(x)u')' + q(x)u = lambda*r(x)*u; u(a)=u(b)=0.",
        'Dùng dấu của phép tính trong mã gốc; sửa nhãn phương trình vốn không nhất quán.',
        'Hàm riêng được chuẩn hóa max |u| = 1; dấu hàm riêng tùy ý.',
        'Đây là trị riêng của hệ sai phân; phải xét giảm h để đánh giá sai số so với bài toán liên tục.'],value=(eig,xx,vec))
    result.table('Trị riêng',['k','lambda'],enumerate(eig,1))
    result.table('Hàm riêng',['x']+[f'u{k+1}' for k in range(count)],np.column_stack([xx,vec]))
    result.curves=[(xx,vec[:,j],f'u{j+1}, lambda={fmt(eig[j])}') for j in range(count)]
    result.presentation={'formulas':[('Bài toán trị riêng',[
        r'-\frac{d}{dx}\left(p(x)\frac{du}{dx}\right)+q(x)u=\lambda r(x)u',r'u(a)=u(b)=0'])],
        'omit_prefixes':('-(p(x)',)}
    return result
