"""Các ca biên khi mang bảng dữ liệu đề thi vào GUI hoặc terminal."""
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from pps import interpolation as ip, calculus as ca, differential as de, chapter1_tools as ch1
from pps.common import fmt, function, spacing, numbers, data_pairs
from pps.gui import TOOLS, dispatch, fields
from pps.math_view import math_png, export_html


class LargeInterpolationTests(unittest.TestCase):
    def test_all_methods_use_more_than_100_nodes_and_analytic_derivatives(self):
        for method in ip.METHODS+['Newton bất kỳ lùi']:
            with self.subTest(method=method):
                x=np.arange(-50.,52. if method=='Bessel' else 51.)
                y=x**3-2*x+1
                r=ip.interpolate(x,y,method,.25,2)
                self.assertEqual(r.presentation['selection']['used_count'],len(x))
                assert_allclose(r.value([-.25,.25,1.5]),np.array([-.25,.25,1.5])**3-2*np.array([-.25,.25,1.5])+1,atol=1e-12)
                self.assertAlmostEqual(r.value.deriv(1)(.25),3*.25**2-2)
                self.assertAlmostEqual(r.value.deriv(2)(.25),1.5)
                self.assertAlmostEqual(r.value.deriv(3)(.25),6)
                for _, formulas in r.presentation['formulas']:
                    for formula in formulas:self.assertTrue(math_png(formula).startswith(b'\x89PNG'))

    def test_501_stirling_nodes_do_not_overflow_factorials(self):
        x=np.arange(-250.,251.)
        r=ip.interpolate(x,x*x,'Stirling',.5,1)
        self.assertEqual(r.presentation['selection']['used_count'],501)
        self.assertAlmostEqual(r.value(.5),.25)
        self.assertAlmostEqual(r.value.deriv()(.5),1)

    def test_large_noisy_polynomial_matches_independent_exact_rational_sum(self):
        # 35 mốc và dữ liệu không phải đa thức bậc thấp: không thể qua nhờ cắt bậc.
        from fractions import Fraction
        x=np.arange(-17.,18.)
        y=np.array([(k*k+3*k)%11-5 for k in range(35)],float)
        z=Fraction(1,3); expected=Fraction(0)
        for i in range(35):
            term=Fraction(int(y[i]))
            for j in range(35):
                if i!=j:term*= (z-int(x[j]))/int(x[i]-x[j])
            expected+=term
        # Cùng số thập phân mà API nhận, thay vì giả định float 1/3 là hữu tỉ đúng.
        r=ip.interpolate(x,y,'Stirling',float(z))
        self.assertAlmostEqual(r.value(float(z)),float(expected),places=11)

    def test_100_node_table_can_select_9_without_changing_data_silently(self):
        x=np.arange(100.)
        y=x*x
        with self.assertRaisesRegex(ValueError,'số mốc lẻ'):
            ip.interpolate(x,y,'Stirling',49.2)
        r=dispatch('Nội suy đa thức','Stirling',{'target':'49.2','count':'9','order':'1'},
                   ' '.join(map(str,x)),' '.join(map(str,y)))
        self.assertEqual(r.presentation['selection'],{'input_count':100,'used_count':9})
        assert_allclose(r.scatter[0],np.arange(45.,54.))
        self.assertAlmostEqual(r.value.deriv()(49.2),98.4)
        self.assertTrue(any('9/100' in line for line in r.lines))
        with self.assertRaises(ValueError):ip.interpolate(x,y,'Stirling',count=9)

    def test_even_gauss_tables_and_minimal_counts(self):
        for count in range(1,11):
            x=np.arange(float(count));y=x**min(3,count-1)
            for method in ('Gauss I','Gauss II'):
                r=ip.interpolate(x,y,method,.25)
                self.assertAlmostEqual(r.value(.25),.25**min(3,count-1))
                self.assertEqual(sorted(r.presentation['difference']['order']),list(range(1,count+1)))

    def test_shifted_and_tiny_coordinates_keep_value_and_derivative(self):
        for method in ip.METHODS:
            x=1e12+np.arange(6. if method=='Bessel' else 5.)
            r=ip.interpolate(x,(x-1e12)**3,method,1e12+1.5,1)
            self.assertAlmostEqual(r.value(1e12+1.5),3.375)
            self.assertAlmostEqual(r.value.deriv()(1e12+1.5),6.75)
        x=np.arange(5.)*1e-150
        r=ip.interpolate(x,x,'Stirling',1.5e-150,1)
        self.assertAlmostEqual(r.value.deriv()(1.5e-150),1.)
        self.assertGreater(r.value(1.5e-150),0.)
        # Chuẩn hóa cũng phải giữ được hai mốc cách nhau rất ít so với toàn dải.
        x=[0.,1e-300,1e300,1.0000000000000002e300,1.0000000000000003e300]
        r=ip.interpolate(x,x,'Newton bất kỳ',.25,1)
        self.assertAlmostEqual(r.value(.25),.25)
        self.assertAlmostEqual(r.value.deriv()(.25),1.)

    def test_invalid_targets_orders_counts_and_wrong_grids(self):
        for kwargs in ({'target':np.nan},{'target':np.inf},{'derivative':-1},{'derivative':.5},
                       {'count':0},{'count':6},{'count':2.5}):
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):
                ip.interpolate(range(5),range(5),'Stirling',**kwargs)
        with self.assertRaises(ValueError):spacing([0,1e-14,2.4e-14])
        with self.assertRaises(ValueError):ip.interpolate([0,1e-14,2.4e-14],[0,1,2],'Stirling')
        self.assertEqual(ip.interpolate([0],[7]).value.deriv(1000)(.2),0)

    def test_compact_report_exports_offline(self):
        x=np.arange(-17.,18.);r=ip.interpolate(x,x*x,'Stirling',.25)
        with tempfile.TemporaryDirectory() as folder:
            target=Path(folder)/'large.html';export_html(r,target)
            content=target.read_text(encoding='utf-8')
            self.assertIn('data:image/png;base64,',content)
            self.assertNotIn('src="http',content)


class NumericalBoundaryTests(unittest.TestCase):
    def test_horner_high_derivative_and_chebyshev_above_old_cap(self):
        r=ip.horner([1,0,0],2,1000)
        self.assertEqual(r.value,4)
        self.assertTrue(any('P^(1000)(2) = 0'==line for line in r.lines))
        x=ip.chebyshev(1200,-2,3).value
        self.assertEqual(len(x),1201)
        assert_allclose(np.cos(1201*np.arccos((x-.5)/2.5)),0,atol=3e-8)
        for args in ((1.5,0,1),(2,0,np.inf),(2,1,1)):
            with self.assertRaises(ValueError):ip.chebyshev(*args)

    def test_small_numbers_and_domain_errors_are_visible(self):
        self.assertEqual(float(fmt(1e-100)),1e-100)
        self.assertNotEqual(fmt(-1e-12),'0')
        for expr,z in [('1/x',0),('sqrt(x)',-1),('log(x)',0),('exp(x)',10000)]:
            with self.subTest(expr=expr),self.assertRaises(ValueError):function(expr)(z)

    def test_import_does_not_silently_drop_invalid_first_data_row(self):
        self.assertEqual(data_pairs([[],['x','y'],[0,1],[1,2]]),[(0.,1.),(1.,2.)])
        self.assertEqual(data_pairs([[0,1],[1,2]]),[(0.,1.),(1.,2.)])
        for rows in ([[0,'bad'],[1,2]],[['x','y'],[0,np.nan]],[[0,1],[2]],[]):
            with self.assertRaises(ValueError):data_pairs(rows)
        assert_allclose(numbers('−1; 2, 3'),[-1,2,3])
        assert_allclose(numbers([1,2]),[1,2])
        r=dispatch('Nội suy đa thức','Stirling',{'target':.5,'order':1},[0,1,2],[0,1,4])
        self.assertAlmostEqual(r.value.deriv()(.5),1)

    def test_newton_append_avoids_product_underflow(self):
        nodes=np.arange(1.,102.)*1e-5
        coefficients=np.zeros(101);coefficients[0]=2
        xx,cc=ch1.newton_append(nodes,coefficients,0.,2.)
        self.assertEqual(cc[-1],0.)
        self.assertEqual(ch1.newton_horner(xx,cc,.5),2.)
        with self.assertRaises(ValueError):ch1.chebyshev_value(1000,2.)

    def test_regression_plot_domain_does_not_discard_valid_fit(self):
        r=ip.least_squares([-2,-1,1,2],[1,1,1,1],basis=['1','sqrt(x^2-1)'])
        assert_allclose(r.value([-2,-1,1,2]),1,atol=1e-12)
        self.assertEqual(r.curves,[])
        self.assertTrue(any('Không vẽ' in line for line in r.lines))
        with self.assertRaises(ValueError):r.value(0)

    def test_bounds_above_factorial_float_limit_and_underflow(self):
        x=np.arange(-100.,101.)
        bound=ch1.interpolation_error_bound(x,.5,1).value
        self.assertGreater(bound,0)
        self.assertTrue(np.isfinite(bound))
        tiny=ch1.chebyshev_global_bound(501,0,1,1).value
        self.assertGreater(tiny,0)  # Không báo cận chính xác bằng 0 do underflow.
        self.assertEqual(ch1.interpolation_error_bound(x,0,1).value,0)

    def test_scaled_linear_systems_do_not_become_singular(self):
        for scale in (1e-100,1.,1e100):
            assert_allclose(ch1.gaussian_elimination(np.eye(2)*scale,np.array([2.,3.])*scale).value,[2,3])
            assert_allclose(ch1.thomas([0,0],np.array([1,2])*scale,[0,0],np.array([2,6])*scale).value,[2,3])

    def test_bisection_endpoints_and_unattainable_float_tolerance(self):
        r=ch1.bisect_model(lambda x:x*x,1,1,2)
        self.assertEqual(r.value,(1,(1,1)))
        a=1e12;b=np.nextafter(a,np.inf)
        with self.assertRaisesRegex(ValueError,'độ phân giải'):
            ch1.bisect_model(lambda x:-1 if x==a else 1,0,a,b,1e-20)
        with self.assertRaises(ValueError):ch1.bisect_model(lambda x:x,0,-1,1,max_steps=0)

    def test_quadrature_rejects_tiny_irregular_tables_and_nonfinite_tolerances(self):
        with self.assertRaises(ValueError):
            ca.integrate('Simpson',a=0,b=2e-14,x=[0,1.4e-14,2e-14],y=[0,1,2])
        for tolerance in (np.nan,np.inf,0,-1):
            with self.assertRaises(ValueError):ca.integrate(tolerance=tolerance)
            with self.assertRaises(ValueError):ca.double_integral('x',0,1,0,1,tolerance=tolerance)
        for fn in (lambda:ca.quadrature_weights(2.5,'Hình thang'),
                   lambda:ca.differentiate([0,1],[0,1],np.nan),
                   lambda:ca.differentiate([0,1],[0,1],.5,degree=1.5)):
            with self.assertRaises(ValueError):fn()

    def test_derivative_weights_resist_constant_cancellation_and_high_degree(self):
        x=np.array([0.,1e-7,2e-7])
        self.assertAlmostEqual(ca.differentiate(x,np.full(3,1e12),.5e-7,2,1).value,0,delta=1e-25)
        x=np.arange(-50.,51.)
        self.assertAlmostEqual(ca.differentiate(x,x**3,.25,100,2).value,1.5,places=12)
        # Đạo hàm cấp hai x^2 bằng 2 dù trọng số cấp hai vượt miền float.
        x=np.arange(3.)*1e-160
        self.assertAlmostEqual(ca.differentiate(x,x*x,1.5e-160,2,2).value,2,places=12)

    def test_ode_reverse_time_zero_interval_and_finite_stopping(self):
        for method in ('Euler hiện','Euler ẩn','Hình thang','RK4','Adams 3'):
            t,y=de.ode(['x'],[math.e],1,0,-.01,method).value
            self.assertEqual(t[-1],0)
            self.assertLess(abs(y[-1,0]-1),.01)
        t,y=de.ode(['x'],[7],2,2,.1).value
        assert_allclose(t,[2]);assert_allclose(y,[[7]])
        with self.assertRaises(ValueError):de.ode(['x'],[1],0,1,.1,tolerance=np.nan)
        with self.assertRaises(ValueError):de.grid(0,1,1e-320)
        with self.assertRaises(ValueError):de.tableau([[0]],1,[0])

    def test_neumann_singularity_and_spline_optional_target(self):
        with self.assertRaisesRegex(ValueError,'Neumann'):
            de.bvp('1','0','1',0,1,.1,(1,0,0),(1,0,0))
        r=dispatch('Spline','Cubic tự nhiên',{'target':''},'0 1 2','0 1 0')
        assert_allclose(r.value([0,1,2]),[0,1,0])
        with self.assertRaises(ValueError):ip.spline([0,1],[0,1],'Cubic kẹp',np.nan,0)


class OfflineDispatchExamples(unittest.TestCase):
    def test_every_method_has_a_working_example_without_creating_a_tk_window(self):
        for tool,(methods,_) in TOOLS.items():
            for method in methods:
                with self.subTest(tool=tool,method=method):
                    values={key:default for key,_,default,_ in fields(tool,method)}
                    x='0 1 2 3 4';y='0 1 4 9 16'
                    if method=='Bessel':x='0 1 2 3';y='0 1 4 9'
                    if method=='Cubic tuần hoàn':y='0 1 0 -1 0'
                    if tool=='Bình phương tối thiểu' and method.startswith('Lũy thừa'):
                        x='1 2 3 4 5';y='1 4 9 16 25'
                    if tool=='Bình phương tối thiểu' and method.startswith('Hàm mũ'):y='1 2 4 8 16'
                    if tool=='Nội suy ngược':x='1 1.1 1.2 1.3 1.4';y='1 1.21 1.44 1.69 1.96'
                    report=dispatch(tool,method,values,x,y)
                    self.assertTrue(report.title)
                    for _, formulas in report.presentation.get('formulas',[]):
                        for formula in formulas:self.assertTrue(math_png(formula).startswith(b'\x89PNG'))


if __name__=='__main__':unittest.main()
