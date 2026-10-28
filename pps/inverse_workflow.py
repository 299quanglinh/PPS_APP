"""Quét mức, chọn mốc theo hướng, chứng nhận và trình bày phép lặp."""
from fractions import Fraction as F
import numpy as np
from .common import Report, points, spacing
from .polynomial_model import nonnegative_integer
from . import inverse_iteration_core as core


def number(v):
    return core.fmt(v, 12)


def latex(v):
    """Phân số ngắn hoặc thập phân gọn; không in tử/mẫu hàng chục chữ số."""
    from decimal import Decimal, localcontext, ROUND_DOWN
    v=core.rat(v)
    if v.denominator==1:return str(v.numerator)
    d=v.denominator
    for factor in (2,5):
        while d%factor==0:d//=factor
    if d==1:
        value=core.fmt(v,12).rstrip('0').rstrip('.')
        if core.rat(value)==v:return value
    if abs(v.numerator)<10000 and v.denominator<10000:
        return rf'\frac{{{v.numerator}}}{{{v.denominator}}}'
    with localcontext() as ctx:
        ctx.prec=max(60,len(str(abs(v.numerator)))+20)
        value=Decimal(v.numerator)/Decimal(v.denominator)
        shown=value.quantize(Decimal('1e-10'),rounding=ROUND_DOWN)
    return str(shown)+r'\ldots'


def scientific_latex(value):
    mantissa,exponent=core.sci(value).split('E')
    mantissa=mantissa.rstrip('0').rstrip('.')
    return mantissa+r'\times10^{'+str(int(exponent))+'}'


def poly_latex(coefficients):
    out=[]
    for i in range(len(coefficients)-1,-1,-1):
        c=coefficients[i]
        if not c:continue
        power='' if i==0 else 't' if i==1 else rf't^{{{i}}}'
        coefficient='' if abs(c)==1 and i else latex(abs(c))
        out.append(('- ' if c<0 else '+ ' if out else '')+coefficient+power)
    return ' '.join(out) or '0'


def newton_latex(model, backward):
    coefficients=model['coefficients']
    terms=[latex(coefficients[0])]
    for r in range(1,len(coefficients)):
        c=coefficients[r]
        if not c:continue
        factors='t'+''.join(rf'(t{"+" if backward else "-"}{j})' for j in range(1,r))
        factor=latex(abs(c)) if r==1 else rf'\frac{{{latex(abs(c))}}}{{{r}!}}'
        terms.append(('- ' if c<0 else '+ ')+factor+r'\,'+factors)
    return ' '.join(terms)


def iteration_workflow(x,y,target,backward=False,degree=2,tolerance=1e-8,max_steps=200,
                       alpha=None,beta=None):
    degree=nonnegative_integer(degree,'Bậc nội suy')
    max_steps=nonnegative_integer(max_steps,'Số bước tối đa')
    if degree<1 or max_steps<1:raise ValueError('Bậc nội suy và số bước tối đa phải dương.')
    x,y,_=points(x,y,2); spacing(x)
    if not np.isfinite(target) or not np.isfinite(tolerance) or tolerance<=0:
        raise ValueError('Mức đích hữu hạn; sai số lặp phải dương, hữu hạn.')
    if (alpha is None)!=(beta is None):raise ValueError('Nhập cả hai đầu đoạn I hoặc để trống cả hai.')
    if alpha is not None and (not np.isfinite(alpha) or not np.isfinite(beta) or alpha>=beta):
        raise ValueError('Đoạn I theo t phải hữu hạn và alpha < beta.')
    # Lưới đầu vào đã qua kiểm tra spacing; triệt sai số biểu diễn float của mốc.
    start,end=core.rat(x[0]),core.rat(x[-1]); h=(end-start)/(len(x)-1)
    xx=[start+i*h for i in range(len(x))]; yy=list(map(core.rat,y)); target=core.rat(target)
    nodes=[i for i,v in enumerate(yy) if v==target]
    intervals=[i for i in range(len(xx)-1) if (yy[i]-target)*(yy[i+1]-target)<0]
    report=Report('Phương pháp lặp nội suy ngược — Newton '+('lùi' if backward else 'tiến'),scatter=(x,y))
    blocks=[]; results=[]
    origin_symbol=r'x_{k+1}' if backward else 'x_k'
    value_symbol=r'y_{k+1}' if backward else 'y_k'
    first_difference=r'\nabla y_{k+1}' if backward else r'\Delta y_k'
    map_symbol=r'\psi' if backward else r'\varphi'
    def text(s):blocks.append(('text',s));report.lines.append(s)
    def head(s):blocks.append(('heading',s))
    def formula(s):
        s=s.replace('x_g',origin_symbol).replace('y_g',value_symbol).replace('c_1',first_difference)
        s=s.replace(r'F\prime',map_symbol+r'\prime').replace('F(',map_symbol+'(')
        s=s.replace(r'\frac',r'\dfrac')
        from .math_view import split_equation
        if '=' in s and len(s)>95:
            lhs,rhs=s.split('=',1)
            for line in split_equation(lhs,rhs,85):blocks.append(('formula',line))
        else:blocks.append(('formula',s))
    def table(name,headers,rows):
        report.table(name,headers,rows);blocks.append(('table',len(report.tables)-1))
    head('Đầu vào')
    text(f'Có {len(x)} mốc; bậc n={degree}, cần đúng {degree+1} mốc mỗi nhánh. '
         f'Sai số lặp cho phép {tolerance:g}; tối đa {max_steps} bước.')
    formula(rf'\bar y={latex(target)},\qquad h={latex(h)}')
    head('Đầu ra')
    text('Các nghiệm trùng mốc và từng khoảng đổi dấu, nghiệm mô hình, chặn sai số, số bước hoặc lý do chưa đạt.')
    head('1. Tìm khoảng chứa nghiệm')
    formula(r'g_i=y_i-\bar y,\qquad g_i=0\quad\mathrm{hoac}\quad g_i g_{i+1}<0')
    table('Quét toàn bộ bảng',['loại','trái','phải'],
          [('Trùng mốc',number(xx[i]),number(xx[i])) for i in nodes]+
          [('Đổi dấu',number(xx[i]),number(xx[i+1])) for i in intervals])
    text('Khoảng đổi dấu có ít nhất một nghiệm nếu f liên tục. Tính duy nhất của hàm gốc cần giả thiết đơn điệu trên khoảng; bảng rời rạc chưa chứng minh điều đó. Quét đổi dấu có thể bỏ sót nghiệm tiếp xúc.')
    for i in nodes:results.append(dict(status='TRUNG_MOC',value=float(xx[i]),steps=0,bound=0.))
    if not nodes and not intervals:text('Không phát hiện mốc trùng hoặc khoảng đổi dấu với mức đã cho.')
    for k in intervals:
        head(f'Nhánh ({number(xx[k])}; {number(xx[k+1])})')
        begin=k+1-degree if backward else k
        stop=k+2 if backward else k+degree+1
        result=dict(interval=(float(xx[k]),float(xx[k+1])),status='CHUA_DAT')
        results.append(result)
        if begin<0 or stop>len(xx):
            result['status']='THIEU_MOC';result['reason']=f'Không đủ {degree+1} mốc về phía '+('trái' if backward else 'phải')+'. Chọn hướng hoặc bậc khác.'
            text(result['reason']);continue
        sx,sy=xx[begin:stop],yy[begin:stop]
        method='lui' if backward else 'tien'
        head('2. Chọn mốc và xây dựng công thức lặp')
        cols=core.difference_table(sy)
        if backward:
            table('Bộ mốc và bảng sai phân lùi',['x_i']+[f'∇^{j} y_i' for j in range(degree+1)],
                  [[number(v)]+[number(c[i-r]) if i>=r else '' for r,c in enumerate(cols)] for i,v in enumerate(sx)])
            text('Newton lùi lấy gốc ở đầu phải khoảng: x_(k+1). Hệ số lấy trên hàng cuối của bảng sai phân lùi.')
        else:
            table('Bộ mốc và bảng sai phân tiến',['x_i']+[f'Δ^{j} y_i' for j in range(degree+1)],
                  [[number(v)]+[number(c[i]) if i<len(c) else '' for c in cols] for i,v in enumerate(sx)])
            text('Newton tiến lấy gốc ở đầu trái khoảng: x_k. Hệ số lấy trên hàng đầu của bảng sai phân tiến.')
        try:model=core.build_model(sx,sy,target,method)
        except ValueError as exc:
            result.update(status='KHONG_HOP_LE',reason=str(exc));text(str(exc));continue
        origin=model['base']
        formula(rf'x_g={latex(origin)},\qquad y_g={latex(model["coefficients"][0])}')
        formula(r't=\frac{x-x_g}{h},\qquad x=x_g+ht')
        formula(rf't=\frac{{x-({latex(origin)})}}{{{latex(h)}}},\qquad x=({latex(origin)})+({latex(h)})t')
        # Newton form as a sum; P below gives all actual numeric coefficients.
        basis=r'\prod_{j=0}^{r-1}(t+j)' if backward else r'\prod_{j=0}^{r-1}(t-j)'
        diff=r'\nabla^r y_g' if backward else r'\Delta^r y_g'
        formula(rf'P_n(x_g+ht)=y_g+c_1t+\sum_{{r=2}}^n\frac{{{basis}}}{{r!}}{diff}')
        formula(rf'P_n(x_g+ht)={newton_latex(model,backward)}')
        text('Thu gọn đa thức theo biến t:')
        formula(rf'P_n(x_g+ht)={poly_latex(model["P"])}')
        formula(r'P_n(x_g+ht)=y_g+c_1t+H(t)')
        formula(rf'H(t)=\sum_{{r=2}}^n\frac{{{basis}}}{{r!}}{diff}')
        formula(r'F(t)=\frac{\bar y-y_g-H(t)}{c_1}')
        formula(rf'c_1={latex(model["coefficients"][1])}')
        formula(rf'F(t)={poly_latex(model["phi"])}')
        t0=(target-model['coefficients'][0])/model['coefficients'][1]
        formula(rf't_0=\frac{{\bar y-y_g}}{{c_1}}={latex(t0)}')
        head('3. Kiểm tra hội tụ')
        # Thử cả khoảng chuẩn trước; nếu chưa đủ, tìm đoạn nhỏ chứa t0.
        lo,hi=(F(-1),F(0)) if backward else (F(0),F(1))
        candidates=[(core.rat(alpha),core.rat(beta))] if alpha is not None else [(lo,hi)]
        if alpha is None:
            t1=core.evaluate(model['phi'],t0)
            for margin in [F(1,10000),F(1,1000),F(1,100),F(1,20)]:
                candidates.append((max(lo,min(t0,t1)-margin),min(hi,max(t0,t1)+margin)))
            for radius in [F(1,1000),F(1,100),F(1,50),F(1,20),F(1,10),F(1,5),F(2,5)]:
                candidates.append((max(lo,t0-radius),min(hi,t0+radius)))
        solved=None;reason='Chưa chứng minh được ánh xạ co trên các đoạn đã thử.'
        for I in candidates:
            if not lo<=I[0]<I[1]<=hi or not I[0]<=t0<=I[1]:
                reason='Đoạn I phải chứa t0 và nằm trong đoạn chuẩn '+('[-1,0].' if backward else '[0,1].');continue
            try:
                proof=core.certify(model['phi'],I,max_depth=4)
                physical=(origin+h*I[0],origin+h*I[1])
                solved=core.solve(sx,sy,target,physical,method,tolerance,proof['q'],max_steps)
                break
            except ValueError as exc:reason=str(exc)
        if solved is None:
            result.update(status='CHUA_CHUNG_NHAN',reason=reason);text(reason)
            text('Chọn đoạn I khác, bậc hoặc hướng khác; chưa nhận kết quả bảo đảm cho nhánh này.');continue
        I=solved['I'];proof=solved['proof']
        formula(rf'I=[{latex(I[0])},{latex(I[1])}],\qquad t_0\in I')
        formula(rf'F(I)\subseteq[{latex(proof["image"][0])},{latex(proof["image"][1])}]\subseteq I')
        formula(rf'\max_{{t\in I}}|F\prime(t)|\leq q={latex(proof["q"])}<1')
        text(f'q ≤ {number(proof["q"])}. Chặn trên toàn đoạn được kiểm tra bằng số học hữu tỉ; các số thập phân hiển thị được làm tròn. Đoạn I là theo t, khác đoạn x.')
        if degree==2:
            sign='+' if backward else '-'
            slope=-model['coefficients'][2]/(2*model['coefficients'][1])
            formula(rf'F\prime(t)=({latex(slope)})(2t{sign}1)')
        text('Phép lặp hội tụ đến nghiệm duy nhất của đa thức trong đoạn I đã chứng nhận.')
        head('4. Lặp và kiểm tra dừng')
        formula(r't_{j+1}=F(t_j),\qquad x_{j+1}=x_g+ht_{j+1}')
        formula(r'B_j=|h|\frac{q}{1-q}|t_j-t_{j-1}|,\qquad B_j\leq\varepsilon')
        text('Tính nội bộ bằng phân số; cập nhật làm tròn đến 40 chữ số. Chặn báo cáo cộng thêm |h|ρ_j/(1−q); định nghĩa ρ_j theo hàm lặp tương ứng ở mục Đánh giá sai số.')
        table('Bảng lặp',['j','t_j','x_j','chặn sai số x'],
              [(j,number(t),number(xv),'—' if E is None else core.sci(E)) for j,t,xv,E,rho in solved['history']])
        j,t,value,E,rho=solved['history'][-1]
        result.update(status='DAT' if solved['success'] else 'CHUA_DAT',value=float(value),steps=j,bound=float(E),q=float(proof['q']))
        head('5. Đánh giá sai số')
        text('Gọi x̂ là nghiệm chính xác của đa thức nội suy trên đoạn đã chứng nhận; x_j là giá trị sau j bước.')
        formula(r'|x_j-\widehat{x}|\leq B_j=|h|\frac{q}{1-q}|t_j-t_{j-1}|')
        formula(r'|x_j-\widehat{x}|\leq A_j=|h|\frac{q^j}{1-q}|t_1-t_0|')
        text('B_j là chặn hậu nghiệm, dùng để dừng; A_j là chặn tiên nghiệm cho phép lặp chính xác. Khi tính cập nhật có làm tròn, chương trình dùng chặn E_j sau đây thay B_j:')
        formula(r'\rho_j=|t_j-F(t_{j-1})|,\qquad E_j=|h|\frac{q|t_j-t_{j-1}|+\rho_j}{1-q}')
        previous_t=solved['history'][-2][1]
        delta=abs(t-previous_t)
        basic=h*proof['q']/(1-proof['q'])*delta
        table('Đánh giá tại bước cuối',['đại lượng','giá trị'],[
            ('Bước j',j),('h',number(h)),('q',number(proof['q'])),
            ('|t_j − t_(j−1)|',core.sci(delta)),
            ('Chặn hậu nghiệm B_j',core.sci(basic)),
            ('Sai số cập nhật ρ_j',core.sci(rho)),
            ('Chặn dùng kiểm tra E_j',core.sci(E)),
            ('Sai số yêu cầu ε',f'{tolerance:g}'),
            ('Phần dư |P_n(x_j) − ȳ|',core.sci(solved['residual']))])
        relation=r'\leq' if solved['success'] else '>'
        formula(rf'E_{{{j}}}={scientific_latex(E)}{relation}\varepsilon={scientific_latex(core.rat(tolerance))}')
        if j>1:
            previous_E=solved['history'][-2][3]
            formula(rf'E_{{{j-1}}}={scientific_latex(previous_E)}>\varepsilon')
        text('Đã đạt điều kiện dừng theo sai số vị trí.' if solved['success'] else
             'Chưa đạt sai số yêu cầu sau số bước tối đa; giá trị cuối chỉ là kết quả tạm.')
        text('Phần dư là sai lệch theo y; không đồng nhất với sai số vị trí theo x.')
        text('Nếu x† là nghiệm hàm gốc, |f−P_n| ≤ E_n và |f′| ≥ m > 0 trên đoạn chứa cả hai nghiệm, thì:')
        formula(r'|x_j-x^\dagger|\leq |x_j-\widehat{x}|+|\widehat{x}-x^\dagger|\leq E_j+\frac{E_n}{m}')
        text('Chưa có cận E_n và m từ dữ liệu này, nên chỉ chứng nhận sai số giải đa thức. E_n phải kể cả sai số bảng nếu muốn đánh giá sai số tổng.')
        head('6. Kết luận nhánh')
        text(f'{result["status"]}: x ≈ {number(value)}; {j} bước; chặn sai số lặp nội bộ {core.sci(E)}.')
        displayed=core.rat(number(value))
        display_rounding=abs(displayed-value)
        shown_bound=E+display_rounding
        text(f'Khi ghi x với 12 chữ số thập phân như trên, cộng sai số làm tròn {core.sci(display_rounding)}; '
             f'chặn cho số đã ghi là {core.sci(shown_bound)}.')
        result['displayed_bound']=float(shown_bound)
        result['residual']=float(solved['residual'])
    head('Tổng hợp kết quả')
    table('Kết quả từng nhánh',['trạng thái','x','bước','chặn sai số'],
          [(r['status'],number(r['value']) if 'value' in r else '—',r.get('steps','—'),f'{r["bound"]:.6g}' if 'bound' in r else r.get('reason','')) for r in results])
    text('Sai số trên dành cho nghiệm đa thức nội suy đã chọn, chưa bao gồm sai số nội suy và dữ liệu. Không tự giảm bậc hay đổi hướng khi thiếu mốc hoặc chưa chứng nhận hội tụ.')
    report.value=results;report.presentation={'blocks':blocks,'inverse_results':results,'math_size':16}
    return report
