"""Đánh giá đa thức nội suy mà không bắt buộc khai triển theo lũy thừa x."""
from decimal import Decimal, localcontext

import numpy as np


def nonnegative_integer(value, label):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 0:
        raise ValueError(f'{label} phải là số nguyên không âm.')
    return int(value)


def working_precision(x):
    """Chọn số chữ số theo số mốc và độ chênh giữa toàn dải với bước nhỏ nhất."""
    with localcontext() as ctx:
        ctx.prec=1000
        ordered=sorted(Decimal(str(float(v))) for v in x)
        span=ordered[-1]-ordered[0]
        gap=min((ordered[i+1]-ordered[i] for i in range(len(ordered)-1)),default=Decimal(1))
        spread=max(0,span.adjusted()-gap.adjusted()+2) if span else 0
    return max(50,2*len(x)+30+spread)


def decimal_cell(value):
    """Giữ ô số bình thường; dùng chuỗi khoa học khi float tràn hoặc hụt số."""
    number=float(value)
    return number if np.isfinite(number) and (number!=0 or value==0) else format(value,'.17g')


class NewtonPolynomial:
    """Cùng đa thức qua mọi mốc; dùng số thập phân và Horner lồng.

    Số chữ số làm việc tăng theo số mốc để hạn chế triệt tiêu số khi bậc cao.
    Cách này không loại bỏ độ nhạy của nội suy đối với sai số trong dữ liệu.
    """

    def __init__(self, x, y):
        self.x = np.asarray(x, float).copy()
        self.y = np.asarray(y, float).copy()
        # Một bảng có cả mốc 1e-300 và 1e300 cần thêm chữ số để phép
        # chuẩn hóa không làm hai mốc phân biệt biến thành cùng một số.
        self.precision = working_precision(x)
        self.order = 0
        # Xếp từ giữa ra hai phía để giảm độ lớn các tích gần tâm bảng.
        center = (len(x) - 1) // 2
        indices = sorted(range(len(x)), key=lambda i: (abs(i-center), -i))
        with localcontext() as ctx:
            ctx.prec = self.precision
            raw = [Decimal(str(float(x[i]))) for i in indices]
            self.origin = raw[0]
            self.scale = max(abs(v-self.origin) for v in raw) or Decimal(1)
            self.nodes = tuple((v-self.origin)/self.scale for v in raw)
            work = [Decimal(str(float(y[i]))) for i in indices]
            coefficients = [work[0]]
            for k in range(1, len(raw)):
                work = [(work[j+1]-work[j])/(self.nodes[j+k]-self.nodes[j])
                        for j in range(len(work)-1)]
                coefficients.append(work[0])
            self.coefficients = tuple(coefficients)
            self.raw_nodes = tuple(raw)

    def deriv(self, m=1):
        m = nonnegative_integer(m, 'Cấp đạo hàm')
        view = object.__new__(type(self))
        view.__dict__ = dict(self.__dict__, order=self.order+m)
        return view

    def __call__(self, value):
        values = np.asarray(value, dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError('Điểm tính đa thức phải hữu hạn.')
        if self.order >= len(self.nodes):
            return np.zeros_like(values)[()] if values.ndim == 0 else np.zeros_like(values)
        answers = []
        with localcontext() as ctx:
            ctx.prec = self.precision
            for v in values.flat:
                z = (Decimal(str(float(v)))-self.origin)/self.scale
                # Horner đồng thời cho đạo hàm; cập nhật từ cấp cao xuống thấp.
                d = [self.coefficients[-1]] + [Decimal(0)] * self.order
                for k in range(len(self.nodes)-2, -1, -1):
                    factor = z-self.nodes[k]
                    for j in range(self.order, 0, -1):
                        d[j] = factor*d[j] + j*d[j-1]
                    d[0] = self.coefficients[k] + factor*d[0]
                answer = float(d[self.order]/self.scale**self.order)
                if not np.isfinite(answer):
                    raise ValueError('Giá trị đa thức/đạo hàm vượt miền số thực máy; chọn nhóm mốc nhỏ hơn hoặc đổi đơn vị dữ liệu.')
                answers.append(answer)
        result = np.asarray(answers).reshape(values.shape)
        return result[()] if values.ndim == 0 else result

    @property
    def coef(self):
        """Hệ số bậc tăng; chỉ chuyển sang float khi còn biểu diễn được."""
        if self.order>=len(self.nodes):return np.array([0.])
        with localcontext() as ctx:
            ctx.prec = self.precision
            c = [self.coefficients[-1]]
            for k in range(len(self.nodes)-2, -1, -1):
                root = self.raw_nodes[k]
                nxt = [Decimal(0)] * (len(c)+1)
                for j, v in enumerate(c):
                    nxt[j] -= root*v/self.scale
                    nxt[j+1] += v/self.scale
                nxt[0] += self.coefficients[k]
                c = nxt
            for _ in range(self.order):
                c = [j*c[j] for j in range(1, len(c))] or [Decimal(0)]
            result = np.asarray([float(v) for v in c])
        if not np.all(np.isfinite(result)):
            raise ValueError('Hệ số theo x vượt miền float; dùng dạng Newton lồng và gọi đa thức trực tiếp.')
        return result

    def degree(self):
        return max(0,len(self.nodes)-1-self.order)
