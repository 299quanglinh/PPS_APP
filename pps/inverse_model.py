"""Dạng sai phân của phép lặp nội suy ngược, không tràn giai thừa."""
from decimal import Decimal, localcontext

from .common import spacing
from .polynomial_model import working_precision


def evaluate(coefficients, value):
    result=Decimal(0)
    for coefficient in reversed(coefficients):result=coefficient+value*result
    return result


def compose_interval(coefficients, center, radius):
    result=[coefficients[-1]]
    for coefficient in reversed(coefficients[:-1]):
        nxt=[Decimal(0)]*(len(result)+1)
        for i,v in enumerate(result):
            nxt[i]+=center*v;nxt[i+1]+=radius*v
        nxt[0]+=coefficient;result=nxt
    return result


class InverseMap:
    def __init__(self,x,y,target,backward=False,precision=50):
        exponents=[Decimal(str(abs(float(v)))).adjusted() for v in y if v!=0]
        spread=max(exponents)-min(exponents) if exponents else 0
        self.precision=max(working_precision(x),precision,spread+2*len(x)+30)
        with localcontext() as ctx:
            ctx.prec=self.precision
            self.h=Decimal(str(spacing(x)))
            self.origin=Decimal(str(float(x[-1] if backward else x[0])))
            self.target=Decimal(str(float(target)))
            column=[Decimal(str(float(v))) for v in y]
            self.y0=column[-1] if backward else column[0]
            coefficients=[];factorial=Decimal(1)
            for k in range(len(x)):
                if k:factorial*=k
                coefficients.append((column[-1] if backward else column[0])/factorial)
                column=[column[i+1]-column[i] for i in range(len(column)-1)]
            self.d1=coefficients[1]
            if not self.d1:raise ValueError('Sai phân cấp một tại gốc bằng 0; chọn nhóm mốc hoặc hướng khác.')
            while len(coefficients)>1 and coefficients[-1]==0:coefficients.pop()
            polynomial=[coefficients[-1]]
            for k in range(len(coefficients)-2,-1,-1):
                shift=Decimal(-k if backward else k)
                nxt=[Decimal(0)]*(len(polynomial)+1)
                for j,value in enumerate(polynomial):
                    nxt[j]-=shift*value;nxt[j+1]+=value
                nxt[0]+=coefficients[k];polynomial=nxt
            self.polynomial=polynomial
            phi=[-v/self.d1 for v in polynomial]
            phi[0]+=self.target/self.d1
            if len(phi)<2:phi.append(Decimal(0))
            phi[1]+=1
            while len(phi)>1 and phi[-1]==0:phi.pop()
            self.phi=phi

    def initial(self):return (self.target-self.y0)/self.d1

    def step(self,t):return evaluate(self.phi,t)

    def residual(self,t):return abs(evaluate(self.polynomial,t)-self.target)

    def location(self,t):return self.origin+self.h*t
