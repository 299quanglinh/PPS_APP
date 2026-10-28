import math
import unittest

import numpy as np
from numpy.testing import assert_allclose
from numpy.polynomial import Polynomial

from pps import interpolation as ip
from pps.gui import dispatch
from pps.math_view import math_png
from pps.chapter1_tools import (bisect_model, chebyshev_value, fill_missing,
    chebyshev_global_bound, horner_division, interpolation_error_bound,
    gaussian_elimination, interpolation_global_bound, inverse_by_swap, inverse_controlled,
    monotone_branches, newton_append, newton_horner, original_scale_fit,
    scan_level, spline_at, spline_boundary, spline_degree, thomas, transformed_fit)


class ChapterOneAudit(unittest.TestCase):
    def test_polynomial_auxiliary_algorithms(self):
        self.assertAlmostEqual(chebyshev_value(3, .4), 4*.4**3-3*.4)
        q, remainder = horner_division([-1,-2,0,3], 2)
        assert_allclose(q, [10,6,3]); self.assertAlmostEqual(remainder,19)
        nodes, coeffs = newton_append([0,1],[1,1],2,5)
        assert_allclose(nodes,[0,1,2]); assert_allclose(coeffs,[1,1,1])
        self.assertAlmostEqual(newton_horner(nodes,coeffs,1.5),3.25)
        horner_report=ip.horner([3,0,-2,-1],2,3)
        assert_allclose([table[2][-1][1] for table in horner_report.tables],[19,34,18,3])
        self.assertIn('P^(2)(2) = 36',horner_report.lines)
        self.assertAlmostEqual(interpolation_error_bound([0,.5,1],.25,math.e).value,math.e/128)
        p=ip.interpolate([0,.5,1],[1,math.exp(.5),math.e]).value
        self.assertLessEqual(abs(math.exp(.25)-p(.25)),math.e/128)
        self.assertAlmostEqual(interpolation_global_bound([0,.5,1],math.e,'cách đều').value,
                               math.e*.5**3/12)
        self.assertAlmostEqual(chebyshev_global_bound(3,-1,1,math.e).value,math.e/24)
        # Với 2 mốc, dấu của D_i trong bảng Lagrange phải theo x_i-x_j.
        rows=ip.interpolate([0,1],[1,2],'Lagrange').tables[0][2]
        weights=next(row for row in rows if row[-1].startswith('Hệ số D_i'))
        assert_allclose(weights[:2],[-1,2])

    def test_missing_values_and_contradiction(self):
        assert_allclose(fill_missing([0,1,2,3],[0,None,4,9],2).value,[0,1,4,9])
        with self.assertRaises(ValueError):fill_missing([0,1,2,3,4],[0,1,4,10,None],2,tolerance=1e-12)
        with self.assertRaises(ValueError):fill_missing([0,1,2,3],[0,None,None,9],2)

    def test_inverse_models_are_distinct(self):
        self.assertAlmostEqual(inverse_by_swap([0,1,2],[0,1,4],2).value,5/3)
        root, bracket = bisect_model(Polynomial([0,0,1]),2,1,2,1e-12).value
        self.assertAlmostEqual(root,math.sqrt(2),places=11)
        self.assertLessEqual(bracket[0],math.sqrt(2)); self.assertGreaterEqual(bracket[1],math.sqrt(2))
        roots,intervals,flats=scan_level([0,1,2,3],[1,-1,0,2],0).value
        self.assertEqual(roots,[2]); self.assertEqual(intervals,[(0,1)]); self.assertEqual(flats,[])
        with self.assertRaises(ValueError):bisect_model(lambda x:x*x,2,0,1)
        self.assertEqual(monotone_branches([0,1,2,3,4],[0,2,1,1,3]).value,
                         ([(0,1),(1,2),(3,4)],[(2,3)]))

    def test_controlled_inverse_contracts_only_with_valid_interval(self):
        expected=(-1+math.sqrt(1.2))/.2
        report=inverse_controlled([0,1,2],[0,1.1,2.4],.5,(0,1),
                                  tolerance=1e-10,residual_tolerance=1e-10)
        self.assertAlmostEqual(report.value,expected,places=9)
        self.assertTrue(any('cận nghiệm mô hình' in line for line in report.lines))
        with self.assertRaises(ValueError):
            inverse_controlled([0,1,2],[0,1.1,2.4],.5,(0,.1))

    def test_spline_conditions_and_knot_derivatives(self):
        pp=spline_degree([0,1,2],[0,1,0],4,[1,0,0]).value
        assert_allclose(pp.c[:,0],[0,0,0,1,0])
        assert_allclose(pp.c[:,1],[-2,0,0,1,1])
        self.assertAlmostEqual(spline_at(pp,1,3),0)
        pp2=spline_boundary([0,1,2],[0,1,0],4,
                            [('trái',1,1),('trái',2,0),('trái',3,0)]).value
        assert_allclose(pp2.c,pp.c,atol=1e-11)
        quadratic=spline_boundary([0,1,2],[0,1,0],2,
                                  [([('trái',1,1),('phải',1,1)],-4)]).value
        assert_allclose(quadratic.c[:,0],[1,0,0],atol=1e-11)
        natural=ip.spline([0,1,2],[0,1,0]).value
        self.assertAlmostEqual(spline_at(natural,.5),.6875)
        self.assertAlmostEqual(spline_at(natural,.5,1),1.125)
        linear=ip.spline([0,1,2],[0,1,0],'Tuyến tính').value
        with self.assertRaises(ValueError):spline_at(linear,1,1)
        with self.assertRaises(ValueError):spline_at(linear,-.1)

    def test_linear_solvers_detect_bad_pivots(self):
        assert_allclose(thomas([0,-1,-1],[2,2,2],[-1,-1,0],[0,0,4]).value,[1,2,3])
        assert_allclose(gaussian_elimination([[0,1],[1,1]],[1,3]).value,[2,1])
        with self.assertRaises(ValueError):thomas([0,1],[0,1],[1,0],[1,1])
        with self.assertRaises(ValueError):gaussian_elimination([[1,2],[2,4]],[1,2])

    def test_weighted_regression_and_domain(self):
        x=np.array([0,1,2]); y=np.array([1,2,2])
        fit=ip.least_squares(x,y,weights=[1,2,1])
        assert_allclose(fit.value(np.array([0,1])),[1.25,1.75],atol=1e-12)
        self.assertTrue(any('MSE' in line for line in fit.lines))
        with self.assertRaises(ValueError):ip.least_squares(x,y,weights=[1,0,1])
        with self.assertRaises(ValueError):ip.least_squares([0,1,2],[1,-2,3],'Hàm mũ')

    def test_new_spline_modes_require_proper_data(self):
        pp=ip.spline([0,1,2,3],[0,1,0,-1],'Cubic not-a-knot').value
        assert_allclose(pp([0,1,2,3]),[0,1,0,-1],atol=1e-12)
        with self.assertRaises(ValueError):ip.spline([0,1,2],[0,1,0],'Cubic not-a-knot')
        with self.assertRaises(ValueError):ip.spline([0,1,2],[0,1,2],'Cubic tuần hoàn')

    def test_transformed_models_restore_original_parameters(self):
        x=np.array([1.,2.,3.,4.])
        cases=[('1/(a+b*x)',1/(2+3*x)),
               ('x/(a+b*x)',x/(2+3*x)),
               ('a+b*ln(x)',2+3*np.log(x)),
               ('A*b^x',2*3**x)]
        for model,y in cases:
            with self.subTest(model=model):
                fit=transformed_fit(x,y,model)
                assert_allclose(fit.value(x),y,rtol=1e-11,atol=1e-11)
        with self.assertRaises(ValueError):transformed_fit([0,1],[1,2],'a+b*ln(x)')

    def test_original_scale_objective_differs_from_log_fit(self):
        x=np.array([0.,1.,2.,3.]);y=np.array([2.,3.1,5.2,8.9])
        log_fit=ip.least_squares(x,y,'Hàm mũ',['x'])
        direct=original_scale_fit(x,y,'Hàm mũ',['x'])
        self.assertLessEqual(np.sum((direct.value(x)-y)**2),
                             np.sum((log_fit.value(x)-y)**2)+1e-10)
        self.assertGreater(np.sum((log_fit.value(x)-y)**2),
                           np.sum((direct.value(x)-y)**2)+1e-5)
        power=original_scale_fit([1,2,3],[2,8,18],'Lũy thừa')
        assert_allclose(power.value(np.array([1.,2.,3.])),[2,8,18],rtol=1e-10)

    def test_gui_dispatch_reaches_new_paths(self):
        r=dispatch('Nội suy ngược','Hàm ngược',{'target':'2'},'0 1 2','0 1 4')
        self.assertAlmostEqual(r.value,5/3)
        r=dispatch('Nội suy ngược','Chia đôi đa thức',
                   {'target':'2','a':'1','b':'2','tol':'1e-10'},'0 1 2','0 1 4')
        self.assertAlmostEqual(r.value[0],math.sqrt(2),places=8)
        r=dispatch('Nội suy ngược','Lặp có cận co',
                   {'target':'1.3225','a':'1.1','b':'1.2','tol':'1e-10','residual_tol':'1e-10','direction':'Tiến'},
                   '1 1.1 1.2 1.3 1.4','1 1.21 1.44 1.69 1.96')
        self.assertAlmostEqual(r.value,1.15,places=8)
        r=dispatch('Spline','Bậc 4 — đạo hàm trái',
                   {'derivatives':'1 0 0','target':'1.5'},'0 1 2','0 1 0')
        self.assertAlmostEqual(r.value(1.5),1.375)
        r=dispatch('Spline','Cubic tuần hoàn',{'target':'1.5'},
                   '0 1 2 3 4','0 1 0 -1 0')
        self.assertAlmostEqual(r.value(0),r.value(4))
        for _, formulas in r.presentation['formulas']:
            for formula in formulas:self.assertTrue(math_png(formula).startswith(b'\x89PNG'))
        r=dispatch('Bình phương tối thiểu','Tuyến tính theo cơ sở',
                   {'basis':'1;x','weights':'1 2 1'},'0 1 2','1 2 2')
        self.assertAlmostEqual(r.value(0),1.25)
        r=dispatch('Bình phương tối thiểu','Hàm mũ — SSE thang y',
                   {'basis':'1;x','initial':''},'0 1 2 3','2 3.1 5.2 8.9')
        self.assertTrue(np.isfinite(r.value(1.5)))
        r=dispatch('Bình phương tối thiểu','Lũy thừa — SSE thang y',
                   {'initial':''},'1 2 3','2 8 18')
        self.assertAlmostEqual(r.value(2),8,places=9)

    def test_course_and_exam_examples(self):
        # Trang 52: 8 mốc quanh 0.64, Bessel và đạo hàm theo x.
        x=np.round(np.arange(0,2.01,.1),10)
        y=x**3-2*x+1
        xx,yy=ip.select_points(x,y,.64,8).value
        assert_allclose(xx,np.arange(.3,1.01,.1),atol=1e-12)
        p=ip.interpolate(xx,yy,'Bessel').value
        self.assertAlmostEqual(p(.64),-.017856,places=11)
        self.assertAlmostEqual(p.deriv()(.64),-.7712,places=10)
        # Đề 20241 câu a: so với đạo hàm giải tích, kể cả hệ số x² bằng 0.
        c=2.57
        lines=ip.horner([3.5,-2.78,5.63,0,-15.28,1.54],c,3).lines
        self.assertAlmostEqual(float(next(s for s in lines if s.startswith('P^(2)')).split('=')[-1]),
                               70*c**3-33.36*c**2+33.78*c,places=8)
        self.assertAlmostEqual(float(next(s for s in lines if s.startswith('P^(3)')).split('=')[-1]),
                               210*c**2-66.72*c+33.78,places=8)


if __name__ == '__main__': unittest.main()
