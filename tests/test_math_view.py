import tempfile
import unittest
from pathlib import Path
import numpy as np
from pps.interpolation import interpolate
from pps.math_view import table_display,export_html,math_png,difference_image


class MathViewTests(unittest.TestCase):
    def test_equal_newton_presents_t_without_expanding_x(self):
        for method,origin in [('Newton tiến','x_0'),('Newton lùi','x_n')]:
            r=interpolate([9,10,11],[2,3,6],method)
            self.assertIn(f'P({origin}+th)',r.text())
            self.assertNotIn('P(x) =',r.text())
            self.assertNotIn('Đa thức theo biến x',[title for title,_ in r.presentation['formulas']])
            self.assertIn('Đa thức theo biến t',[title for title,_ in r.presentation['formulas']])
            np.testing.assert_allclose(r.value([9,10,11]),[2,3,6])

    def test_gauss_indices_and_selected_differences(self):
        x=np.arange(5.)
        r=interpolate(x,x**3,'Gauss I',1.2)
        m=r.presentation['difference']
        self.assertEqual(m['offset'],2)
        self.assertEqual(m['order'],[5,3,1,2,4])
        self.assertEqual(m['selected'],{(2,1),(1,2),(1,3),(0,4)})
        self.assertEqual(m['a'][2,1],19)
        im=difference_image(m);self.assertGreater(im.width,500)
        self.assertFalse(any('P^{(0)}' in f for _,fs in r.presentation['formulas'] for f in fs))

    def test_lagrange_table_follows_the_original_app_workflow(self):
        r=interpolate([0,1,2],[1,2,5],'Lagrange')
        self.assertEqual(len(r.tables),1)
        name,headers,rows=table_display(r,0)
        self.assertEqual(name,'Các bước lập đa thức Lagrange')
        self.assertEqual(headers,['Cột 0','Cột 1','Cột 2','Cột 3','Ghi chú'])
        notes=[row[-1] for row in rows]
        self.assertIn('Ma trận (x_j - x_i)',notes)
        self.assertIn('Bảng tích',notes)
        self.assertIn('Bảng thương W_i = w_(n+1)/(x-x_i)',notes)
        self.assertIn('Hệ số đa thức nội suy P(x)',notes)

    def test_two_tables_all_polynomial_methods(self):
        from pps.interpolation import METHODS
        for method in METHODS:
            x=np.arange(6. if method=='Bessel' else 5.)
            r=interpolate(x,x**3,method)
            self.assertEqual(len(r.tables),1 if method=='Lagrange' else 2,method)
            self.assertNotIn('Các hệ số và tích phụ',r.text())
            if method!='Lagrange':
                self.assertIn('bậc giảm dần',r.text())
                self.assertEqual(r.tables[1][1][1],f'bậc {len(x)-1}')
                self.assertEqual(len(r.tables[1][2]),len(x))
            np.testing.assert_allclose(r.value(x),x**3,atol=1e-10)

    def test_unscaled_products_and_factorial_formula(self):
        from numpy.polynomial import Polynomial
        import math
        for method in ('Newton tiến','Newton lùi','Gauss I','Gauss II','Stirling','Bessel'):
            x=np.arange(6. if method=='Bessel' else 5.)
            y=x**3-2*x+4
            r=interpolate(x,y,method)
            rows=r.tables[1][2]
            for k,row in enumerate(rows):
                self.assertEqual(row[1:][len(x)-1-k],1,method)
            for _,formulas in r.presentation['after_tables']:
                for k,f in enumerate(formulas):
                    self.assertIn(str(k)+'!',f)
                    self.assertTrue(math_png(f).startswith(b'\x89PNG'))
            self.assertIn('/2!',r.text())
            if method=='Gauss I':
                self.assertEqual(rows[4][1:],[1,-2,-1,2,0])
                a=r.presentation['difference']['a'];c=2
                built=sum((a[c-k//2,k]/math.factorial(k))*Polynomial(row[:0:-1]) for k,row in enumerate(rows))
                np.testing.assert_allclose(built(x-c),y,atol=1e-10)

    def test_math_parses_and_exports_offline(self):
        r=interpolate([0,1,2,3],[0,1,8,27],'Bessel',1.5)
        for _,forms in r.presentation['formulas']:
            for formula in forms:self.assertTrue(math_png(formula).startswith(b'\x89PNG'))
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'loi_giai.html';export_html(r,path)
            text=path.read_text(encoding='utf-8')
            self.assertIn('data:image/png;base64,',text)
            self.assertNotIn('<script',text)
            self.assertNotIn('src="http',text)


if __name__=='__main__':unittest.main()
