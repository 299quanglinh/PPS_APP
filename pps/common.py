"""Input validation, restricted mathematical expressions and reports."""
import ast
import math
from dataclasses import dataclass, field
import numpy as np
import sympy as sp


DISPLAY_DIGITS = 7

def set_display_digits(digits):
    global DISPLAY_DIGITS
    digits=int(digits)
    if not 0 <= digits <= 15: raise ValueError("Số chữ số thập phân phải từ 0 đến 15.")
    DISPLAY_DIGITS=digits

def fmt(value, digits=None):
    if digits is None: digits=DISPLAY_DIGITS
    if isinstance(value, str):
        return value
    if value is None:
        return ''
    v = float(value)
    if not math.isfinite(v):
        return str(v)
    if abs(v) >= 1e12 or (digits and 0 < abs(v) < 10.**(-digits)):
        return f'{v:.{max(1,digits)}g}'
    s = f'{v:.{digits}f}'
    if digits: s=s.rstrip('0').rstrip('.')
    return '0' if s in ('', '-0') else s


def numbers(text):
    try:
        if isinstance(text,str):
            values=np.array([float(v) for v in text.replace('−','-').replace(',', ' ').replace(';', ' ').split()])
        else:values=np.asarray(text,dtype=float)
    except (TypeError,ValueError) as exc:
        raise ValueError('Dãy phải chứa các số thực hữu hạn; dùng dấu chấm cho phần thập phân.') from exc
    if values.ndim!=1 or not len(values) or not np.all(np.isfinite(values)):
        raise ValueError('Nhập các số hữu hạn, cách nhau bằng dấu cách; dấu thập phân là dấu chấm.')
    return values


def data_pairs(rows):
    """Đọc hai cột x,y; chỉ bỏ dòng tiêu đề nhận diện được và dòng rỗng."""
    data=[];seen=False
    x_headers={'x','xi','x_i','t','t_i','year','năm'}
    y_headers={'y','yi','y_i','f(x)','f(t)','f','giá trị'}
    for index,row in enumerate(rows,1):
        if not row or all(v is None or str(v).strip()=='' for v in row):continue
        labels=[str(v).strip().casefold() for v in row[:2]]
        if not seen and len(labels)==2 and labels[0] in x_headers and labels[1] in y_headers:
            seen=True;continue
        seen=True
        try:
            pair=(float(row[0]),float(row[1]))
            if not np.all(np.isfinite(pair)):raise ValueError()
        except (ValueError,TypeError,IndexError) as exc:
            raise ValueError(f'Dòng {index} không có đủ hai số x,y hữu hạn ở hai cột đầu; không tự bỏ dòng dữ liệu lỗi.') from exc
        data.append(pair)
    if not data:raise ValueError('Không có dữ liệu. Hai cột đầu phải là x và y; tiêu đề có thể là x,y.')
    return data


def points(x, y, minimum=1):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) < minimum:
        raise ValueError(f'Dãy x và y phải cùng độ dài, có ít nhất {minimum} điểm.')
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError('Dữ liệu phải hữu hạn.')
    indices = np.argsort(x, kind='stable')
    x, y = x[indices], y[indices]
    if np.any(np.diff(x) == 0):
        raise ValueError('Mốc x trùng: nội suy/spline cần hoành độ phân biệt; hãy xử lý dữ liệu trước.')
    return x, y, 0


def spacing(x):
    if len(x) < 2:
        return 1.0
    gaps = np.diff(x)
    h = gaps[0]
    # Dung sai tuyệt đối cố định sẽ nhận nhầm lưới rất nhỏ nhưng không đều.
    if not np.isfinite(h) or h <= 0 or not np.all(np.isfinite(gaps)) or not np.allclose(gaps, h, rtol=1e-8, atol=0):
        raise ValueError('Phương pháp này yêu cầu mốc tăng dần, cách đều.')
    return float(h)


def expression(text, variables=('x',)):
    """Build SymPy expressions from a whitelist AST, never eval user input."""
    symbols = {v: sp.Symbol(v, real=True) for v in variables}
    names = {**symbols, 'pi': sp.pi, 'e': sp.E, 'E': sp.E}
    funcs = {k: getattr(sp, k) for k in ('sin', 'cos', 'tan', 'asin', 'acos', 'atan',
             'sinh', 'cosh', 'tanh', 'exp', 'sqrt', 'log', 'Abs')}
    funcs.update(ln=sp.log, abs=sp.Abs)
    if len(text) > 2000:
        raise ValueError('Biểu thức quá dài.')
    try:
        tree = ast.parse(text.strip().replace('^', '**'), mode='eval')
        if sum(1 for _ in ast.walk(tree)) > 300:
            raise ValueError('Biểu thức quá phức tạp.')
        def visit(node):
            if isinstance(node, ast.Constant) and type(node.value) in (int, float):
                if not math.isfinite(node.value):
                    raise ValueError('Hằng số ngoài miền hỗ trợ.')
                return sp.sympify(node.value)
            if isinstance(node, ast.Name) and node.id in names:
                return names[node.id]
            if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
                return visit(node.operand) * (-1 if isinstance(node.op, ast.USub) else 1)
            if isinstance(node, ast.BinOp):
                a, b = visit(node.left), visit(node.right)
                if isinstance(node.op, ast.Add): return a+b
                if isinstance(node.op, ast.Sub): return a-b
                if isinstance(node.op, ast.Mult): return a*b
                if isinstance(node.op, ast.Div): return a/b
                if isinstance(node.op, ast.Pow):
                    if b.is_number and abs(float(b)) > 1000:
                        raise ValueError('Số mũ quá lớn.')
                    return a**b
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                    and node.func.id in funcs and not node.keywords and len(node.args) == 1):
                return funcs[node.func.id](visit(node.args[0]))
            raise ValueError('Chỉ dùng biến, số, + - * / ^ và các hàm toán học được hỗ trợ.')
        expr = visit(tree.body)
        if expr.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
            raise ValueError('Biểu thức không hữu hạn.')
        return expr, tuple(symbols.values())
    except (SyntaxError, TypeError, OverflowError) as exc:
        raise ValueError('Sai cú pháp hàm. Ví dụ: sin(x) + x^2; phép nhân phải ghi *.') from exc


def function(text, variables=('x',)):
    expr, symbols = expression(text, variables)
    raw = sp.lambdify(symbols, expr, modules='numpy')
    def evaluate(*args):
        if any(not np.all(np.isfinite(arg)) for arg in args):
            raise ValueError('Đối số của hàm phải hữu hạn.')
        try:
            with np.errstate(all='ignore'):
                result = np.asarray(raw(*args))
        except (ZeroDivisionError, OverflowError, FloatingPointError) as exc:
            raise ValueError('Hàm không có giá trị thực hữu hạn tại điểm đang tính.') from exc
        if np.iscomplexobj(result):
            raise ValueError('Hàm không có giá trị thực hữu hạn tại điểm đang tính.')
        try:result=result.astype(float)
        except (TypeError,ValueError,OverflowError) as exc:
            raise ValueError('Giá trị biểu thức vượt miền số thực float.') from exc
        if not np.all(np.isfinite(result)):
            raise ValueError('Hàm không có giá trị thực hữu hạn tại điểm đang tính.')
        shape = np.broadcast_shapes(*(np.shape(a) for a in args))
        return np.broadcast_to(result, shape)
    return evaluate


def difference(x, y, divided=False):
    n = len(x)
    a = np.full((n, n), np.nan)
    a[:, 0] = y
    for k in range(1, n):
        a[:n-k, k] = np.diff(a[:n-k+1, k-1])
        if divided:
            a[:n-k, k] /= x[k:]-x[:n-k]
    return a


def polynomial_text(poly, variable='x'):
    terms = []
    for k, c in enumerate(poly.coef):
        s = fmt(abs(c))
        if s == '0': continue
        term = s if k == 0 else f'{s}*{variable}' + (f'^{k}' if k > 1 else '')
        terms.append(('-' if c < 0 else '+', term))
    if not terms: return '0'
    return ('-' if terms[0][0] == '-' else '') + terms[0][1] + ''.join(f' {s} {v}' for s,v in terms[1:])


@dataclass
class Report:
    title: str
    lines: list = field(default_factory=list)
    tables: list = field(default_factory=list)
    curves: list = field(default_factory=list)
    scatter: tuple | None = None
    value: object = None
    presentation: dict = field(default_factory=dict)

    def table(self, name, headers, rows):
        self.tables.append((name, list(headers), [list(r) for r in rows]))
        return self

    def text(self):
        if self.presentation.get('blocks'):
            out=[self.title,'='*64]
            for kind,value in self.presentation['blocks']:
                if kind=='table':
                    name,headers,rows=self.tables[value]
                    out.append(Report(name,tables=[(name,headers,rows)]).text())
                else:out.append(str(value))
            return '\n\n'.join(out)
        out = [self.title, '=' * 64, *self.lines]
        for name, headers, rows in self.tables:
            cells=[list(headers)]+[[fmt(v) for v in row] for row in rows[:600]]
            widths=[max(len(row[j]) if j<len(row) else 0 for row in cells) for j in range(len(headers))]
            out.extend(['',name])
            out.extend('  '.join(value.ljust(widths[j]) for j,value in enumerate(row)).rstrip() for row in cells)
            if len(rows) > 600:
                out.append(f'Hiển thị 600/{len(rows)} dòng; xuất Excel/CSV để xem toàn bộ.')
        out.extend(['',*self.presentation.get('after_tables_text',[])])
        return '\n'.join(out)
