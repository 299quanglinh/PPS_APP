"""Offline LaTeX/MathText presentation, independent of numeric storage order."""
import base64
import html
import io
import re
from functools import lru_cache
import tkinter as tk
from tkinter import ttk
import numpy as np
from PIL import Image, ImageTk, ImageDraw, ImageFont
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties
from matplotlib import rc_context
from .common import fmt


@lru_cache(maxsize=512)
def math_png(latex, size=14):
    out=io.BytesIO()
    with rc_context({'savefig.transparent':True}):
        math_to_image('$'+latex+'$',out,prop=FontProperties(size=size,math_fontfamily='stix'),
                      dpi=130,format='png',color='#162d39')
    return out.getvalue()


def math_image(latex,size=14):
    return Image.open(io.BytesIO(math_png(latex,size))).convert('RGBA')


def poly_latex(poly,var='x'):
    terms=[]
    for power in range(len(poly.coef)-1,-1,-1):
        c=poly.coef[power]; value=fmt(abs(c))
        if value=='0':continue
        body=value if power==0 else ('' if value=='1' else value+r'\,')+var+(f'^{{{power}}}' if power>1 else '')
        terms.append(('-' if c<0 else '+',body))
    if not terms:return '0'
    return ('-' if terms[0][0]=='-' else '')+terms[0][1]+''.join(f' {s} {body}' for s,body in terms[1:])


def split_equation(lhs,rhs,width=100):
    """Break only at additive separators emitted by poly_latex, not inside powers."""
    terms=re.split(r' (?=[+-] )',rhs); lines=[]; current=lhs+' = '
    for term in terms:
        if len(current)+len(term)>width and not current.endswith(' = '):
            lines.append(current);current=r'\qquad '+term
        else:current+=term+' '
    if current.strip():lines.append(current.strip())
    return lines


def cubic_spline_presentation(report,x,y,pp,kind,left,right,target):
    """A worked solution following the alpha, theta, v notation in the course notes."""
    from numpy.polynomial import Polynomial
    from .common import polynomial_text
    h=np.diff(x);dy=np.diff(y);d=dy/h;n=len(h)
    matrix=np.zeros((n+1,n+1));rhs=np.zeros(n+1)
    for i in range(1,n):
        matrix[i,i-1:i+2]=[h[i-1],2*(h[i-1]+h[i]),h[i]]
        rhs[i]=6*(d[i]-d[i-1])
    clamped=kind=='Cubic kẹp'
    if clamped:
        matrix[0,:2]=[2*h[0],h[0]];rhs[0]=6*(d[0]-left)
        matrix[-1,-2:]=[h[-1],2*h[-1]];rhs[-1]=6*(right-d[-1])
    else:
        matrix[0,0]=matrix[-1,-1]=1
        rhs[0]=0 if kind=='Cubic tự nhiên' else left
        rhs[-1]=0 if kind=='Cubic tự nhiên' else right
    alpha=np.linalg.solve(matrix,rhs)
    theta=d-h*(alpha[1:]-alpha[:-1])/6
    v=y[:-1]-alpha[:-1]*h*h/6
    blocks=[]
    def heading(s):blocks.append(('heading',s))
    def body(s):blocks.append(('text',s))
    def formula(s):blocks.append(('formula',s))
    def table(name,headers,rows):
        report.table(name,headers,rows);blocks.append(('table',len(report.tables)-1))
    def equation(row,b):
        terms=[]
        for j,c in enumerate(row):
            if c:
                term=('' if fmt(abs(c))=='1' else fmt(abs(c)))+rf'\alpha_{{{j}}}'
                terms.append(('-' if c<0 else '+',term))
        return (('-' if terms[0][0]=='-' else '')+terms[0][1]+''.join(' '+s+' '+t for s,t in terms[1:])+'='+fmt(b))
    heading('1. Sắp các cặp dữ liệu theo x tăng dần')
    body('Giữ nguyên từng cặp (x, y); các hệ số chỉ làm tròn khi hiển thị.')
    body(f'Có {n+1} mốc x_0 đến x_{n}, tạo thành {n} đoạn và {n-1} nút trong; n = {n}.')
    table('Bảng dữ liệu',['i','x_i','y_i'],[[i,x[i],y[i]] for i in range(n+1)])
    heading('2. Tính khoảng mốc và hệ số góc')
    body(f'Chỉ số đoạn: i = 0, ..., n−1 (trong bài: 0 đến {n-1}).')
    formula(r'h_i=x_{i+1}-x_i,\quad \Delta y_i=y_{i+1}-y_i,\quad d_i=\frac{\Delta y_i}{h_i}')
    table('Các đại lượng trên từng đoạn',['i','h_i',r'\Delta y_i','d_i'],[[i,h[i],dy[i],d[i]] for i in range(n)])
    heading('3. Đặt ẩn và biểu diễn spline theo tài liệu')
    formula(r'\alpha_i=S^{\prime\prime}(x_i),\qquad i=0,\ldots,n')
    body(f'Mỗi nhánh S_i nằm trên [x_i, x_(i+1)]; các hệ số θ_i, v_i dùng i = 0, ..., n−1 (0 đến {n-1}).')
    formula(r'S_i(x)=-\frac{\alpha_i}{6h_i}(x-x_{i+1})^3+\frac{\alpha_{i+1}}{6h_i}(x-x_i)^3+\theta_i(x-x_i)+v_i')
    formula(r'\theta_i=d_i-\frac{h_i}{6}(\alpha_{i+1}-\alpha_i),\qquad v_i=y_i-\frac{\alpha_i h_i^2}{6}')
    heading('4. Lập hệ từ tính liên tục của đạo hàm bậc nhất')
    body(f'Chỉ xét nút trong: i = 1, ..., n−1 (trong bài: 1 đến {n-1}). Không thế i = 0 vào công thức này.')
    formula(r'h_{i-1}\alpha_{i-1}+2(h_{i-1}+h_i)\alpha_i+h_i\alpha_{i+1}=6(d_i-d_{i-1})')
    formula(r'i=1,\ldots,n-1')
    for i in range(1,n):formula(equation(matrix[i],rhs[i]))
    heading('5. Ghép điều kiện biên và giải hệ')
    if clamped:
        formula(rf'S^\prime({fmt(x[0])})={fmt(left)},\qquad S^\prime({fmt(x[-1])})={fmt(right)}')
        formula(r'2h_0\alpha_0+h_0\alpha_1=6(d_0-S^\prime(x_0))')
        formula(r'h_{n-1}\alpha_{n-1}+2h_{n-1}\alpha_n=6(S^\prime(x_n)-d_{n-1})')
    else:
        formula(rf'\alpha_0=S^{{\prime\prime}}({fmt(x[0])})={fmt(rhs[0])}')
        formula(rf'\alpha_{{{n}}}=S^{{\prime\prime}}({fmt(x[-1])})={fmt(rhs[-1])}')
    formula(equation(matrix[0],rhs[0]));formula(equation(matrix[-1],rhs[-1]))
    table('Hệ phương trình Aα = b',[rf'\alpha_{{{i}}}' for i in range(n+1)]+['b'],[list(row)+[b] for row,b in zip(matrix,rhs)])
    # Show the actual tridiagonal elimination used to recover the solution.
    diag=np.diag(matrix).copy();reduced=rhs.copy()
    body('Khử lần lượt phần tử dưới đường chéo, rồi thế ngược từ phương trình cuối.')
    for i in range(1,n+1):
        factor=matrix[i,i-1]/diag[i-1]
        diag[i]-=factor*matrix[i-1,i];reduced[i]-=factor*reduced[i-1]
    upper=np.diag(diag)+np.diag(np.diag(matrix,1),1)
    table('Hệ sau khử tiến',[rf'\alpha_{{{i}}}' for i in range(n+1)]+['b'],[list(row)+[b] for row,b in zip(upper,reduced)])
    for i in range(n,-1,-1):formula(rf'\alpha_{{{i}}}={fmt(alpha[i])}')
    heading('6. Tính θ_i và v_i')
    body(f'Tính theo từng đoạn: i = 0, ..., n−1 (trong bài: 0 đến {n-1}); bước này có dùng i = 0.')
    formula(rf'\theta_0=({fmt(d[0])})-\frac{{{fmt(h[0])}}}{{6}}[({fmt(alpha[1])})-({fmt(alpha[0])})]={fmt(theta[0])}')
    formula(rf'v_0={fmt(y[0])}-\frac{{({fmt(alpha[0])})({fmt(h[0])})^2}}{{6}}={fmt(v[0])}')
    table('Hệ số theo ký hiệu tài liệu',['i',r'\alpha_i',r'\alpha_{i+1}',r'\theta_i','v_i'],[[i,alpha[i],alpha[i+1],theta[i],v[i]] for i in range(n)])
    heading('7. Thế hệ số và viết từng nhánh spline')
    for i in range(n):
        body(f'Đoạn [{fmt(x[i])}, {fmt(x[i+1])}], đặt u_{i} = x - ({fmt(x[i])}).')
        formula(rf'S_{{{i}}}(x)=-\frac{{({fmt(alpha[i])})}}{{6({fmt(h[i])})}}(x-({fmt(x[i+1])}))^3+\frac{{({fmt(alpha[i+1])})}}{{6({fmt(h[i])})}}(x-({fmt(x[i])}))^3')
        formula(rf'\qquad +({fmt(theta[i])})(x-({fmt(x[i])}))+({fmt(v[i])})')
        for tex in split_equation(rf'S_{{{i}}}(x)',poly_latex(Polynomial(pp.c[:,i][::-1]),rf'u_{{{i}}}')):formula(tex)
    if target is not None:
        heading('8. Thay giá trị')
        formula(rf'S({fmt(target)})={fmt(pp(target))}')
        if not x[0]<=target<=x[-1]:body('Điểm tính ngoài bảng: ngoại suy bằng đoạn đầu/cuối.')
    report.presentation={'blocks':blocks}


def difference_model(method,x,y,a,center):
    central=method in ('Gauss I','Gauss II','Stirling','Bessel')
    n=len(x);offset=center if central else 0;selected=set()
    for k in range(1,n):
        if method=='Gauss I':selected.add((center-k//2,k))
        elif method=='Gauss II':selected.add((center-(k+1)//2,k))
        elif method=='Newton tiến' or method=='Newton bất kỳ':selected.add((0,k))
        elif method in ('Newton lùi','Newton bất kỳ lùi'):selected.add((n-k-1,k))
        elif method=='Stirling':
            if k%2:selected.update(((center-(k+1)//2,k),(center-(k-1)//2,k)))
            else:selected.add((center-k//2,k))
        elif method=='Bessel':
            if k%2:selected.add((center-k//2,k))
            else:selected.update(((center-k//2,k),(center-k//2+1,k)))
    order=[0]*n
    if method in ('Gauss I','Gauss II'):
        indices=[center]
        for r in range(1,n):
            indices.extend(j for j in ([center+r,center-r] if method=='Gauss I' else [center-r,center+r]) if 0<=j<n)
        for rank,j in enumerate(indices,1):order[j]=rank
    return dict(method=method,x=list(x),y=list(y),a=a,offset=offset,selected=selected,order=order)


def difference_image(model):
    """Render forward differences from their starting row, divided differences from their ending row."""
    x,y,a=model['x'],model['y'],model['a'];n=len(x);offset=model['offset'];divided=model['method'].startswith('Newton bất kỳ')
    # Bound rendered size; full numeric table remains available in export.
    count=min(n,18);degree=min(n-1,10);cw=145;rh=82;top=52
    widths=[58,140,140]+[cw]*degree
    # Measure both labels and numbers so long indices/values cannot overlap a cell.
    for k in range(1,degree+1):
        labels=[rf'f[x_i,\ldots,x_{{i+{k}}}]' if divided else rf'\Delta^{{{k}}}y_i']
        for j in range(min(count,n-k)):
            idx=j-offset
            labels.extend([rf'f[x_{{{idx}}},\ldots,x_{{{idx+k}}}]' if divided else rf'\Delta^{{{k}}}y_{{{idx}}}', fmt(a[j,k])])
        widths[k+2]=max(cw,max(math_image(v,13).width for v in labels)+32)
    widths[1]=max(140,max(math_image(fmt(v),13).width for v in x[:count])+32)
    widths[2]=max(140,max(math_image(fmt(v),13).width for v in y[:count])+32)
    image=Image.new('RGB',(sum(widths),top+rh*count+16),'white');draw=ImageDraw.Draw(image)
    xedges=np.cumsum([0]+widths).tolist()
    draw.rectangle((0,0,image.width,top),fill='#c3e1e3')
    for j in range(count):
        draw.rectangle((0,top+j*rh,image.width,top+(j+1)*rh),fill='#edf6f7' if j%2==0 else '#e2f0f2')
    for edge in xedges:draw.line((edge,0,edge,image.height),fill='white',width=2)
    def paste(tex,cx,cy,size=13,highlight=False):
        im=math_image(tex,size)
        xx=int(cx-im.width/2);yy=int(cy-im.height/2)
        if highlight:draw.rectangle((xx-6,yy-5,xx+im.width+6,yy+im.height+5),fill='#fff7db',outline='#186c67',width=2)
        image.paste(im,(xx,yy),im)
    ranked=model['method'] in ('Gauss I','Gauss II')
    headers=[r'\mathrm{STT}' if ranked else 'i','x_i','y_i']+[rf'f[x_i,\ldots,x_{{i+{k}}}]' if divided else (r'\Delta y_i' if k==1 else rf'\Delta^{{{k}}}y_i') for k in range(1,degree+1)]
    for c,tex in enumerate(headers):paste(tex,(xedges[c]+xedges[c+1])/2,top/2)
    for j in range(count):
        mid=top+(j+.5)*rh;idx=j-offset
        paste(str(model['order'][j] if ranked else idx),widths[0]/2,mid)
        paste(rf'x_{{{idx}}}',(xedges[1]+xedges[2])/2,mid-14)
        paste(fmt(x[j]),(xedges[1]+xedges[2])/2,mid+16)
        paste(rf'y_{{{idx}}}',(xedges[2]+xedges[3])/2,mid-14)
        paste(fmt(y[j]),(xedges[2]+xedges[3])/2,mid+16)
        for k in range(1,min(degree,n-j-1)+1):
            # The original PPS grid anchors a divided difference f[x_i,...,x_{i+k}]
            # at its final node x_{i+k}, making the triangle converge downward.
            # Finite forward differences deliberately keep their starting-node row.
            row=j+k if divided else j
            if row>=count:continue
            cx=(xedges[k+2]+xedges[k+3])/2;cy=top+(row+.5)*rh
            symbol=rf'f[x_{{{idx}}},\ldots,x_{{{idx+k}}}]' if divided else (rf'\Delta y_{{{idx}}}' if k==1 else rf'\Delta^{{{k}}}y_{{{idx}}}')
            paste(symbol,cx,cy-13,12, (j,k) in model['selected'])
            paste(fmt(a[j,k]),cx,cy+16,12)
    return image


def interpolation_presentation(result,method,x,y,a,q,p,origin,h,center,normalized,target,derivative):
    formulas=[]
    if normalized:
        formulas += [('Đổi biến theo tài liệu',[rf'x_0={fmt(origin)},\quad h={fmt(h)},\quad t=\frac{{x-x_0}}{{h}}=\frac{{x-({fmt(origin)})}}{{{fmt(h)}}}'])]
    templates={
      'Lagrange':r'P(x)=\sum_{i=0}^{n}y_i\prod_{j=0,\ j\ne i}^{n}\frac{x-x_j}{x_i-x_j}',
      'Newton bất kỳ':r'P(x)=y_0+\sum_{k=1}^{n}f[x_0,\ldots,x_k]\prod_{j=0}^{k-1}(x-x_j)',
      'Newton tiến':r'P(x_0+th)=y_0+\sum_{k=1}^{n}\frac{\Delta^k y_0}{k!}\prod_{j=0}^{k-1}(t-j)',
      'Newton lùi':r'P(x_n+th)=y_n+\sum_{k=1}^{n}\frac{\Delta^k y_{n-k}}{k!}\prod_{j=0}^{k-1}(t+j)',
      'Gauss I':r'P(x_0+th)=y_0+t\Delta y_0+\frac{t(t-1)}{2!}\Delta^2y_{-1}+\frac{t(t-1)(t+1)}{3!}\Delta^3y_{-1}+\cdots',
      'Gauss II':r'P(x_0+th)=y_0+t\Delta y_{-1}+\frac{t(t+1)}{2!}\Delta^2y_{-1}+\frac{t(t+1)(t-1)}{3!}\Delta^3y_{-2}+\cdots',
      'Stirling':r'P(x_0+th)=y_0+t\frac{\Delta y_{-1}+\Delta y_0}{2}+\frac{t^2}{2!}\Delta^2y_{-1}+\cdots',
      'Bessel':r'P(x_0+th)=\frac{y_0+y_1}{2}+\left(t-\frac{1}{2}\right)\Delta y_0+\frac{t(t-1)}{2!}\frac{\Delta^2y_{-1}+\Delta^2y_0}{2}+\cdots'
    }
    if method=='Newton lùi':
        formulas[0]=('Đổi biến theo tài liệu',[rf'x_n={fmt(origin)},\quad h={fmt(h)},\quad t=\frac{{x-x_n}}{{h}}=\frac{{x-({fmt(origin)})}}{{{fmt(h)}}}'])
    templates['Newton bất kỳ lùi']=r'P(x)=y_n+\sum_{k=1}^{n}f[x_{n-k},\ldots,x_n]\prod_{j=0}^{k-1}(x-x_{n-j})'
    template=templates[method]
    if method in ('Gauss I','Gauss II','Stirling','Bessel'):
        chunks=template.split('+\\frac')
        display=[chunks[0]]+[r'\qquad +\frac'+part for part in chunks[1:]]
    else:display=[template]
    formulas.append(('Công thức '+method,display))
    if normalized:formulas.append(('Đa thức theo biến t',split_equation('P(x_n+th)' if method=='Newton lùi' else 'P(x_0+th)',poly_latex(q,'t'))))
    if method not in ('Newton tiến','Newton lùi'):
        formulas.append(('Đa thức theo biến x',split_equation('P(x)',poly_latex(p))))
    if target is not None:
        val=result.value.deriv(derivative)(target)
        label='P' if derivative==0 else "P'" if derivative==1 else "P''" if derivative==2 else rf'P^{{({derivative})}}'
        formulas.append(('Thay giá trị',[rf'{label}({fmt(target)})={fmt(val)}']))
    result.presentation=dict(formulas=formulas,
                             omit_prefixes=('x_0 =','x_n =','P(x_0+th)','P(x_n+th)','P(x) =','P^('))
    if method != 'Lagrange':
        result.presentation['difference'] = difference_model(method,x,y,a,center)


def product_formula(report, method, products, center):
    """Show the unscaled basis from the table with factorials kept outside."""
    from numpy.polynomial import Polynomial
    from .common import polynomial_text
    def delta(k,j):
        return rf'y_{{{j}}}' if k==0 else rf'\Delta^{{{k}}}y_{{{j}}}'
    tex=[];plain=[]
    for k,coefficients in enumerate(products):
        if method=='Gauss I':j=-k//2 if k%2==0 else -(k//2)
        elif method=='Gauss II':j=-((k+1)//2)
        elif method=='Newton lùi':j=len(products)-1-k
        else:j=0
        if method=='Stirling':
            j=-(k//2)
            c=rf'\frac{{{delta(k,j-1)}+{delta(k,j)}}}{{2}}' if k%2 else delta(k,j)
        elif method=='Bessel':
            j=-(k//2)
            c=delta(k,j) if k%2 else rf'\frac{{{delta(k,j)}+{delta(k,j+1)}}}{{2}}'
        else:c=delta(k,j)
        poly=Polynomial(coefficients)
        term=rf'\left({c}\right)\frac{{{poly_latex(poly,"t")}}}{{{k}!}}'
        tex.append(('P(t) = ' if k==0 else r'\qquad + ')+term)
        # Plain text is deliberately explicit about averaging and factorials.
        def d(order,index):return f'y[{index}]' if order==0 else f'Delta^{order} y[{index}]'
        if method=='Stirling' and k%2:cp=f'({d(k,j-1)} + {d(k,j)})/2'
        elif method=='Bessel' and k%2==0:cp=f'({d(k,j)} + {d(k,j+1)})/2'
        else:cp=d(k,j)
        plain.append(('P(t) = ' if k==0 else '     + ')+f'({cp}) * ({polynomial_text(poly,"t")})/{k}!')
    note='P(t) ở đây là đa thức theo biến chuẩn hóa t; mỗi hàng bảng tích chưa chia giai thừa và chưa nhân sai phân.'
    report.presentation['after_tables']=[('Lập P(t) từ bảng tích — giữ giai thừa ở mẫu',tex)]
    report.presentation['after_tables_text']=[note,*plain]


def symbolic_latex(text,variables=('x',)):
    import sympy as sp
    from .common import expression
    expr,_=expression(text,variables)
    return sp.latex(expr)


def text_math(text):
    """Only convert known generated mathematical cells; never evaluate free prose."""
    s=str(text).strip()
    if not s:return None
    if re.fullmatch(r'[-+]?\d+(\.\d+)?[Ee][-+]?\d+',s):
        mantissa,exponent=re.split('[Ee]',s)
        return mantissa+r'\times10^{'+str(int(exponent))+'}'
    if re.fullmatch(r'[-+0-9.eE ]+',s):return s
    if re.fullmatch(r'[0-9xyzt()+*/^ .\-]+',s) and re.search('[xyzt]',s):
        s=re.sub(r'\^(-?\d+)',r'^{\1}',s).replace('*',r'\,')
        return s
    known={'x':'x_i','y':'y_i','k':'k','i':'i','t':'t','lambda':r'\lambda','hệ số':r'a_k',
           'hệ số lẻ':r'a_{2i+1}','hệ số chẵn':r'a_{2i}','trọng số':r'w_i',
           'f(x)':r'f(x_i)','t trước':r't_{k-1}','t sau':r't_k','độ lệch':r'|t_k-t_{k-1}|'}
    if s in known:return known[s]
    if s in ('x_i','y_i','t_j','x_j'):return s
    if re.fullmatch(r'[Δ∇]\^\d+ y_i',s):
        symbol=r'\nabla' if s.startswith('∇') else r'\Delta'
        return symbol+'^{'+s.split('^')[1].split()[0]+r'}y_i'
    if re.fullmatch(r'[a-z]_\{[0-9a-z+\-]+\}',s):return s
    if re.fullmatch(r'bậc \d+',s):return 't^{'+s.split()[-1]+'}'
    if re.fullmatch(r'[a-zA-Z0-9().+*/^ _\-]+',s):
        try:return symbolic_latex(s,('x','t','y','z'))
        except (ValueError,TypeError):pass
    return None


def table_display(report,index):
    name,headers,rows=report.tables[index]
    if name.startswith('Bảng tích:'):
        # The numeric engine stores ascending coefficients. Present highest degree first.
        variable='x' if report.title in ('Lagrange','Newton bất kỳ','Newton bất kỳ lùi') else 't'
        headers=[headers[0]]+[f'{variable}^{k}' if k>1 else variable if k==1 else '1' for k in range(len(headers)-2,-1,-1)]
        if 'tăng dần' in name:
            rows=[[r[0]]+list(reversed(r[1:])) for r in rows]
        if report.title == 'Lagrange':
            name='Bảng tích tử số cơ sở Lagrange — bậc giảm dần'
        elif report.title.startswith('Newton bất kỳ'):
            name='Bảng tích đa thức cơ sở Newton — bậc giảm dần'
        else:
            name='Bảng tích — chưa chia giai thừa, bậc giảm dần'
    elif name=='Các hệ số và tích phụ' and report.title=='Stirling':
        headers=['i','a_{2i+1}','a_{2i+2}','tích']
    return name,headers,rows


def table_image(report,index,limit=80):
    if index==0 and report.presentation.get('difference'):
        return difference_image(report.presentation['difference'])
    name,headers,rows=table_display(report,index);headers=headers[:12];rows=[r[:12] for r in rows[:limit]]
    font=ImageFont.truetype('C:/Windows/Fonts/segoeui.ttf',17)
    cells=[]
    for row in [headers]+rows:
        rr=[]
        for value in row:
            text=fmt(value);tex=text_math(text)
            if tex:
                rr.append(math_image(tex,12))
            else:
                box=font.getbbox(text or ' ');im=Image.new('RGBA',(max(1,box[2])+4,30),(255,255,255,0))
                ImageDraw.Draw(im).text((0,0),text,font=font,fill='#162d39');rr.append(im)
        cells.append(rr)
    widths=[max(95,max(row[j].width if j<len(row) else 0 for row in cells)+24) for j in range(len(headers))]
    heights=[max(44,max(im.height for im in row)+18) for row in cells]
    result=Image.new('RGB',(sum(widths),sum(heights)),'white');draw=ImageDraw.Draw(result);yy=0
    for i,row in enumerate(cells):
        xx=0
        for j,im in enumerate(row):
            draw.rectangle((xx,yy,xx+widths[j],yy+heights[i]),fill='#c3e1e3' if i==0 else '#edf6f7' if i%2 else 'white',outline='#d6e4e6')
            result.paste(im,(int(xx+(widths[j]-im.width)/2),int(yy+(heights[i]-im.height)/2)),im);xx+=widths[j]
        yy+=heights[i]
    return result


class MathDocument(ttk.Frame):
    def __init__(self,parent):
        super().__init__(parent)
        self.text=tk.Text(self,wrap='none',bg='white',fg='#162d39',font=('Segoe UI',11),padx=20,pady=18,borderwidth=0)
        vs=ttk.Scrollbar(self,command=self.text.yview);hs=ttk.Scrollbar(self,orient='horizontal',command=self.text.xview)
        self.text.configure(yscrollcommand=vs.set,xscrollcommand=hs.set)
        self.text.grid(row=0,column=0,sticky='nsew');vs.grid(row=0,column=1,sticky='ns');hs.grid(row=1,column=0,sticky='ew')
        self.rowconfigure(0,weight=1);self.columnconfigure(0,weight=1);self.images=[]
        self.text.tag_configure('title',font=('Segoe UI',18,'bold'),spacing3=12)
        self.text.tag_configure('heading',font=('Segoe UI',12,'bold'),spacing1=15,spacing3=8)
        self.text.tag_configure('body',spacing3=6)

    def begin(self):
        self.text.configure(state='normal');self.text.delete('1.0','end');self.images=[]

    def picture(self,im):
        photo=ImageTk.PhotoImage(im,master=self);self.images.append(photo)
        self.text.image_create('end',image=photo);self.text.insert('end','\n\n')

    def message(self,message):
        self.begin();self.text.insert('end',message,'body');self.text.configure(state='disabled')

    def show(self,report):
        self.begin();self.text.insert('end',report.title+'\n','title')
        if report.presentation.get('blocks'):
            for kind,value in report.presentation['blocks']:
                if kind=='formula':
                    im=math_image(value,report.presentation.get('math_size',14))
                    padded=Image.new('RGBA',(im.width+24,im.height+16),(255,255,255,0))
                    padded.paste(im,(12,8))
                    self.picture(padded)
                elif kind=='table':
                    self.text.insert('end',report.tables[value][0]+'\n','heading')
                    self.picture(table_image(report,value))
                else:self.text.insert('end',value+'\n','heading' if kind=='heading' else 'body')
            self.text.configure(state='disabled');self.text.yview_moveto(0);self.text.xview_moveto(0)
            return
        for title,formulas in report.presentation.get('formulas',[]):
            self.text.insert('end',title+'\n','heading')
            for formula in formulas:self.picture(math_image(formula))
        prefixes=report.presentation.get('omit_prefixes',())
        for line in report.lines:
            if not line.startswith(prefixes):self.text.insert('end',line+'\n','body')
        for i in range(len(report.tables)):
            name,headers,rows=table_display(report,i)
            self.text.insert('end',name+'\n','heading')
            if i==0 and report.presentation.get('difference'):
                if report.title.startswith('Newton bất kỳ'):
                    self.text.insert('end','Mốc x tăng từ trên xuống; tỷ sai phân f[x_i,\\ldots,x_{i+k}] cùng hàng với mốc cuối x_{i+k}, chụm xuống dưới.\nÔ viền xanh: tỷ sai phân dùng trong công thức Newton.\n','body')
                else:
                    self.text.insert('end','Mốc x tăng từ trên xuống; sai phân cùng hàng với mốc bắt đầu, chụm lên trên.\nÔ viền xanh: sai phân dùng trong công thức (hai ô khi lấy trung bình).\nVới Gauss, STT là thứ tự kết nạp mốc trong bài giảng.\n','body')
            self.picture(table_image(report,i))
            if len(rows)>80 or len(headers)>12 or (i==0 and report.presentation.get('difference') and len(rows)>11):self.text.insert('end','Bảng trình bày rút gọn; chọn Dữ liệu dạng lưới hoặc xuất Excel để xem đầy đủ.\n','body')
        for title,formulas in report.presentation.get('after_tables',[]):
            self.text.insert('end',title+'\n','heading')
            self.text.insert('end',report.presentation['after_tables_text'][0]+'\n','body')
            for formula in formulas:self.picture(math_image(formula))
        self.text.configure(state='disabled');self.text.yview_moveto(0);self.text.xview_moveto(0)


class MathTable(ttk.Frame):
    def __init__(self,parent):
        super().__init__(parent)
        self.canvas=tk.Canvas(self,bg='white',highlightthickness=0)
        vs=ttk.Scrollbar(self,command=self.canvas.yview);hs=ttk.Scrollbar(self,orient='horizontal',command=self.canvas.xview)
        self.canvas.configure(yscrollcommand=vs.set,xscrollcommand=hs.set)
        self.canvas.bind('<MouseWheel>',lambda e:self.canvas.yview_scroll(-int(e.delta/120),'units'))
        self.canvas.bind('<Shift-MouseWheel>',lambda e:self.canvas.xview_scroll(-int(e.delta/120),'units'))
        self.canvas.grid(row=0,column=0,sticky='nsew');vs.grid(row=0,column=1,sticky='ns');hs.grid(row=1,column=0,sticky='ew')
        self.columnconfigure(0,weight=1);self.rowconfigure(0,weight=1)

    def show(self,report,index):
        self.canvas.delete('all');self.photo=ImageTk.PhotoImage(table_image(report,index),master=self)
        self.canvas.create_image(16,16,image=self.photo,anchor='nw')
        self.canvas.configure(scrollregion=self.canvas.bbox('all'));self.canvas.xview_moveto(0);self.canvas.yview_moveto(0)


def export_html(report,path):
    def embed(im):
        out=io.BytesIO();im.save(out,format='PNG')
        return '<img alt="Bảng hoặc công thức toán" src="data:image/png;base64,'+base64.b64encode(out.getvalue()).decode()+'">'
    parts=['<!doctype html><html lang="vi"><meta charset="utf-8"><title>'+html.escape(report.title)+'</title>',
           '<style>body{font:17px Segoe UI,sans-serif;color:#162d39;max-width:1200px;margin:40px auto;padding:20px}h2{margin-top:28px}img{display:block;margin:12px 0}section{overflow:auto} @media print{img{max-width:100%;height:auto}section{overflow:visible}}</style>',
           '<h1>'+html.escape(report.title)+'</h1>']
    if report.presentation.get('blocks'):
        for kind,value in report.presentation['blocks']:
            if kind=='formula':parts.append('<section>'+embed(math_image(value))+'</section>')
            elif kind=='table':parts.append('<h2>'+html.escape(report.tables[value][0])+'</h2><section>'+embed(table_image(report,value))+'</section>')
            else:
                tag='h2' if kind=='heading' else 'p'
                parts.append('<'+tag+'>'+html.escape(value)+'</'+tag+'>')
        parts.append('</html>');path.write_text('\n'.join(parts),encoding='utf-8')
        return
    for title,formulas in report.presentation.get('formulas',[]):
        parts.append('<h2>'+html.escape(title)+'</h2><section>')
        parts.extend(embed(math_image(f)) for f in formulas);parts.append('</section>')
    for line in report.lines:
        if not line.startswith(report.presentation.get('omit_prefixes',())):parts.append('<p>'+html.escape(line)+'</p>')
    for i in range(len(report.tables)):
        name,_,rows=table_display(report,i);parts.append('<h2>'+html.escape(name)+'</h2><section>'+embed(table_image(report,i))+'</section>')
        if len(rows)>11 or len(report.tables[i][1])>12:parts.append('<p>Bảng trình bày có thể rút gọn; xuất Excel để lấy toàn bộ dữ liệu.</p>')
    for title,formulas in report.presentation.get('after_tables',[]):
        parts.append('<h2>'+html.escape(title)+'</h2><p>'+html.escape(report.presentation['after_tables_text'][0])+'</p><section>')
        parts.extend(embed(math_image(f)) for f in formulas);parts.append('</section>')
    parts.append('</html>');path.write_text('\n'.join(parts),encoding='utf-8')
