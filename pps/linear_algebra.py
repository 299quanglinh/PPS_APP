"""Bộ giải cơ sở cho hồi quy; giữ thang trọng số và cột thiết kế tách biệt."""
from decimal import Decimal, localcontext

import numpy as np


def sum_products(*vectors):
    """Tổng các tích hữu hạn, kể cả khi tổng/trung gian vượt miền float."""
    with localcontext() as ctx:
        ctx.prec=80
        terms=[]
        for values in zip(*vectors):
            term=Decimal(1)
            for value in values:term*=Decimal(str(float(value)))
            terms.append(term)
        nonzero=[v.adjusted() for v in terms if v]
        if nonzero:ctx.prec=max(80,max(nonzero)-min(nonzero)+30)
        return sum(terms,Decimal(0))


def fit_basis(theta, values, weights=None):
    """SVD sau đổi thang; QR Decimal khi cột quá gần phụ thuộc.

    Không lập hệ chuẩn để giải. Trọng số được chia cùng một hằng số,
    do đó không đổi bài toán cực tiểu.
    """
    theta=np.asarray(theta,float);values=np.asarray(values,float)
    if theta.ndim!=2 or values.shape!=(theta.shape[0],) or theta.shape[0]<theta.shape[1] or theta.shape[1]<1:
        raise ValueError('Số quan sát phải ít nhất bằng số cơ sở; ma trận và vế phải phải khớp cỡ.')
    if not np.all(np.isfinite(theta)) or not np.all(np.isfinite(values)):
        raise ValueError('Ma trận hồi quy và dữ liệu phải hữu hạn.')
    weights=np.ones(len(values)) if weights is None else np.asarray(weights,float)
    if weights.shape!=values.shape or not np.all(np.isfinite(weights)) or np.any(weights<=0):
        raise ValueError('Trọng số hồi quy phải dương, hữu hạn và cùng số quan sát.')
    sqrt_w=np.sqrt(weights/np.max(weights))
    with np.errstate(over='ignore',under='ignore',invalid='ignore'):
        scales=np.max(abs(theta),axis=0)
        if np.any(scales==0):raise ValueError('Có hàm cơ sở bằng 0 trên mọi quan sát.')
        design=(theta/scales)*sqrt_w[:,None]
        yscale=float(np.max(abs(values))) or 1.
        target=(values/yscale)*sqrt_w
    try:
        solution,_,rank,singular=np.linalg.lstsq(design,target,rcond=None)
    except np.linalg.LinAlgError as exc:
        raise ValueError('SVD hồi quy không hội tụ; kiểm tra cơ sở và tỉ lệ dữ liệu.') from exc
    high_precision=rank<theta.shape[1] or singular[-1]<singular[0]*1e-10 or np.any(sqrt_w==0)
    if not high_precision:
        with localcontext() as ctx:
            ctx.prec=80
            coefficients=np.array([float(Decimal(str(float(v)))*Decimal(str(yscale))/Decimal(str(float(s))))
                                   for v,s in zip(solution,scales)])
    else:
        all_values=np.r_[theta.ravel(),values,weights]
        nonzero=abs(all_values[all_values!=0])
        spread=int(np.log10(np.max(nonzero))-np.log10(np.min(nonzero))) if len(nonzero) else 0
        with localcontext() as ctx:
            ctx.prec=max(80,spread+60,4*theta.shape[1]+40)
            max_weight=Decimal(str(float(np.max(weights))))
            roots=[(Decimal(str(float(w)))/max_weight).sqrt() for w in weights]
            columns=[[Decimal(str(float(v)))*root for v,root in zip(column,roots)] for column in theta.T]
            target=[Decimal(str(float(v)))*root for v,root in zip(values,roots)]
            count=len(columns);q=[];r=[[Decimal(0)]*count for _ in range(count)]
            for j,column in enumerate(columns):
                v=column.copy();norm=ctx.sqrt(sum((a*a for a in v),Decimal(0)))
                # Hai lượt trực giao hóa giảm mất chữ số khi cột gần nhau.
                for _ in range(2):
                    for i,qi in enumerate(q):
                        projection=sum((a*b for a,b in zip(qi,v)),Decimal(0))
                        r[i][j]+=projection
                        v=[a-projection*b for a,b in zip(v,qi)]
                remaining=ctx.sqrt(sum((a*a for a in v),Decimal(0)))
                if remaining<=norm*Decimal(10)**(-ctx.prec+20):
                    raise ValueError('Các hàm cơ sở phụ thuộc tuyến tính trên tập dữ liệu.')
                r[j][j]=remaining;q.append([a/remaining for a in v])
            z=[sum((a*b for a,b in zip(qi,target)),Decimal(0)) for qi in q]
            solution=[Decimal(0)]*count
            for i in range(count-1,-1,-1):
                solution[i]=(z[i]-sum((r[i][j]*solution[j] for j in range(i+1,count)),Decimal(0)))/r[i][i]
            coefficients=np.array([float(v) for v in solution])
    if not np.all(np.isfinite(coefficients)):
        raise ValueError('Hệ số hồi quy vượt miền float; chọn cơ sở có đổi gốc hoặc đơn vị phù hợp.')
    return coefficients,high_precision


def exp_parameter(log_value, sign=1.):
    """Giữ tham số ở dạng exp(...) nếu tham số riêng lẻ không vừa float."""
    with np.errstate(over='ignore',under='ignore'):value=float(sign*np.exp(log_value))
    if np.isfinite(value) and value!=0:return value
    return ('-' if sign<0 else '')+f'exp({float(log_value):.17g})'
