import math
import unittest
import numpy as np
from numpy.testing import assert_allclose
from pps import interpolation as ip, calculus as ca, differential as de
from pps.common import function, points, fmt


class InterpolationTests(unittest.TestCase):
    def test_lagrange_decimal_input_keeps_exact_display_coefficients(self):
        x=[1.2,1.5,1.7,1.8,2.1,2.3]
        y=[.892,1.179,1.358,1.445,1.688,1.839]
        report=ip.interpolate(x,y,'Lagrange')
        self.assertAlmostEqual(report.value(0), -34011/5500, places=13)

    def test_all_polynomial_methods_reproduce_known_polynomial(self):
        for method in ip.METHODS:
            with self.subTest(method=method):
                x=np.linspace(-1,1,6 if method=='Bessel' else 7)
                y=2*x**4-3*x**3+.7*x-1
                result=ip.interpolate(x,y,method,target=.21,derivative=2)
                xx=np.linspace(-1,1,31)
                assert_allclose(result.value(xx),2*xx**4-3*xx**3+.7*xx-1,atol=1e-11)
                assert_allclose(result.value.deriv(2)(.21),24*.21**2-18*.21,atol=1e-10)

    def test_arbitrary_nodes(self):
        x=np.array([0,.2,.9,1.8,3.])
        for method in ('Lagrange','Newton bất kỳ'):
            p=ip.interpolate(x,x**3-x,method).value
            assert_allclose(p([.1,1.2,2.5]),np.array([.1,1.2,2.5])**3-[.1,1.2,2.5],atol=1e-11)

    def test_random_data_interpolates_every_node(self):
        rng=np.random.default_rng(42)
        for count in range(1,10):
            methods=['Bessel'] if count%2==0 else ['Gauss I','Gauss II','Stirling']
            x=np.arange(count)*.25+2; y=rng.normal(size=count)
            for method in methods:
                with self.subTest(n=count,method=method):
                    r=ip.interpolate(x,y,method)
                    assert_allclose(r.value(x),y,atol=5e-5)

    def test_edge_cases_and_validation(self):
        self.assertEqual(ip.interpolate([4],[7],'Stirling').value(100),7)
        assert_allclose(ip.interpolate([2,3],[5,7],'Bessel').value([2,2.5,3]),[5,6,7])
        for method,x in [('Stirling',[0,1]),('Bessel',[0,1,2]),('Newton tiến',[0,1,3])]:
            with self.assertRaises(ValueError):ip.interpolate(x,x,method)
        x,y,d=points([2,0,2,1],[4,0,99,1]);assert_allclose(y,[0,1,4]);self.assertEqual(d,1)

    def test_inverse(self):
        x=np.linspace(1,1.4,5)
        for backward in (False,True):
            self.assertAlmostEqual(ip.inverse(x,x*x,1.3225,backward).value,1.15,places=8)
        with self.assertRaises(ValueError):ip.inverse([0,1,2],[0,1,0],.5)

    def test_splines(self):
        x=np.arange(5.)
        for kind in ['Tuyến tính','Bậc hai — biên trái','Bậc hai — biên phải','Cubic tự nhiên','Cubic kẹp','Cubic đạo hàm bậc hai']:
            p=ip.spline(x,x*x,kind,0 if kind!='Cubic đạo hàm bậc hai' else 2,8 if kind!='Cubic đạo hàm bậc hai' else 2).value
            assert_allclose(p(x),x*x,atol=1e-12)
            if kind not in ('Tuyến tính','Cubic tự nhiên'):
                assert_allclose(p(np.linspace(0,4,21)),np.linspace(0,4,21)**2,atol=1e-12)

    def test_least_squares(self):
        x=np.linspace(1,3,10)
        for kind,y,basis in [('Tuyến tính theo cơ sở',1+2*x+3*x*x,['1','x','x^2']),
                              ('Lũy thừa',2*x**1.3,['1']),('Hàm mũ',-2*np.exp(.3*x),['x'])]:
            r=ip.least_squares(x,y,kind,basis);assert_allclose(r.value(x),y,rtol=1e-10)
        # Mixed signs must still produce a finite back-transformed fit.
        r=ip.least_squares([-1,0,1,2],[-2,-1,1,4],'Lũy thừa')
        self.assertTrue(np.all(np.isfinite(r.value(np.array([-1,0,1,2])))))

    def test_select_chebyshev_horner(self):
        r=ip.select_points(np.arange(7),np.arange(7)**2,3.2,3)
        assert_allclose(r.value[0],[2,3,4])
        r=ip.select_points([0,1,2,3,4],[0,1,2,1,0],1.5,reverse=True)
        self.assertEqual(len(r.tables),2)
        roots=ip.chebyshev(4,-1,1).value
        assert_allclose(np.cos(5*np.arccos(roots)),0,atol=1e-14)
        self.assertEqual(ip.horner([1,0,-2,1],2,3).value,5)


class CalculusTests(unittest.TestCase):
    def test_derivative_weights(self):
        x=np.array([0,.3,.8,1.2,2.])
        for order,truth in [(1,4*.7**3-2),(2,12*.7**2),(3,24*.7)]:
            self.assertAlmostEqual(ca.differentiate(x,x**4-2*x,.7,4,order).value,truth,places=10)

    def test_quadrature_exactness(self):
        for method in ('Hình thang','Trung điểm','Simpson','Newton-Cotes'):
            self.assertAlmostEqual(ca.integrate(method,'2*x+1',0,1,20).value,2,places=12)
        self.assertAlmostEqual(ca.integrate('Simpson','x^3',0,2,20).value,4,places=12)
        for order in range(1,11):
            self.assertAlmostEqual(ca.integrate('Newton-Cotes',f'x^{order}',0,1,order*2,order).value,1/(order+1),places=9)

    def test_tables_and_refinement(self):
        x=np.linspace(0,1,11)
        self.assertAlmostEqual(ca.integrate('Simpson',a=0,b=1,x=x,y=x*x).value,1/3,places=12)
        x=(np.arange(10)+.5)/10
        self.assertAlmostEqual(ca.integrate('Trung điểm',a=0,b=1,x=x,y=x).value,.5,places=12)
        self.assertAlmostEqual(ca.integrate('Simpson','sin(x)',0,np.pi,10,tolerance=1e-9).value,2,places=8)
        with self.assertRaises(ValueError):ca.integrate('Simpson',n=3)
        with self.assertRaises(ValueError):ca.integrate('Trung điểm',x=[0,1],y=[0,1])
        x=np.linspace(0,2,21)
        self.assertAlmostEqual(ca.integrate('Simpson',a=.4,b=1.4,x=x,y=x*x).value,(1.4**3-.4**3)/3,places=12)
        result=ca.integrate('Hình thang','x^2',0,1,tolerance=1e-5,bound=2)
        self.assertLess(abs(result.value-1/3),1e-5)

    def test_double(self):
        for method in ('Hình thang','Simpson'):
            self.assertAlmostEqual(ca.double_integral('x+y',0,1,0,2,10,method).value,3,places=11)
        self.assertAlmostEqual(ca.double_integral('x^2+y^2',0,1,0,1).value,2/3,places=12)


class DifferentialTests(unittest.TestCase):
    def test_all_ode_methods(self):
        limits={'Euler hiện':.03,'Euler ẩn':.03,'Hình thang':.0001,'RK2':.0002,'RK3':1e-6,'RK4':1e-8}
        limits.update({f'Adams {s}':1e-5 for s in range(2,6)})
        for method,err in limits.items():
            with self.subTest(method=method):
                ts,ys=de.ode(['x'],[1],0,1,.01,method).value
                self.assertLess(abs(ys[-1,0]-math.e),err)
                self.assertEqual(ts[-1],1)

    def test_ode_system_and_custom(self):
        t,y=de.ode(['y','-x'],[1,0],0,1,.01,'RK4').value
        assert_allclose(y[-1],[math.cos(1),-math.sin(1)],atol=1e-9)
        for custom in [de.custom_rk2('alpha2',.5),de.custom_rk2('r1',.5),de.custom_rk2('r2',.5),de.custom_rk2('beta11',1),de.custom_rk3(.5,1)]:
            _,y=de.ode(['x'],[1],0,1,.01,'RK tùy chỉnh',custom=custom).value
            self.assertLess(abs(y[-1,0]-math.e),.0002)
        with self.assertRaises(ValueError):de.grid(0,1,.3)
        with self.assertRaises(ValueError):de.tableau([[1]],[1],[1])

    def test_bvp_dirichlet_and_robin(self):
        x,y=de.bvp('1','0','2',0,1,.05).value
        assert_allclose(y,x*(1-x),atol=1e-12)
        # u=x^2, left u'=0; right u'+u=3.
        x,y=de.bvp('1','0','-2',0,1,.05,(1,0,0),(1,1,3)).value
        assert_allclose(y,x*x,atol=1e-12)

    def test_eigenvalues(self):
        vals,x,vec=de.eigen('1','0','1',0,1,.01,3).value
        expected=4/.01**2*np.sin(np.arange(1,4)*np.pi/200)**2
        assert_allclose(vals,expected,atol=1e-9)
        assert_allclose(vec[[0,-1]],0)
        self.assertLess(abs(vals[0]-np.pi**2),.001)


class InputTests(unittest.TestCase):
    def test_expressions(self):
        assert_allclose(function('sin(x)^2 + cos(x)^2')([0,1]),1)
        self.assertEqual(float(function('e^0 + ln(e)')(0)),2)
        for expr in ["__import__('os')",'x.__class__','open(x)','[x][0]','10^10000']:
            with self.assertRaises(ValueError):function(expr)
        with self.assertRaises(ValueError):function('sqrt(x)')(-1)
        self.assertEqual(fmt(-1e-12),'0')


if __name__=='__main__':unittest.main()
