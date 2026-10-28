"""Kiểm định độc lập các thuật toán ngoài nội suy đa thức."""
import math
import unittest
from decimal import Decimal
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose
from numpy.polynomial import Polynomial

from pps import calculus as ca, differential as de, interpolation as ip, chapter1_tools as ch1


class QuadratureAudit(unittest.TestCase):
    def test_each_rule_has_its_algebraic_degree(self):
        for method,degree in [('Hình thang',1),('Trung điểm',1),('Simpson',3)]:
            for k in range(degree+1):
                with self.subTest(method=method,k=k):
                    r=ca.integrate(method,f'x^{k}',-2,3,12)
                    self.assertAlmostEqual(r.value,(3**(k+1)-(-2)**(k+1))/(k+1),places=10)
        for degree in range(1,11):
            exact_degree=degree+1 if degree%2==0 else degree
            r=ca.integrate('Newton-Cotes',f'x^{exact_degree}',0,1,2*degree,degree)
            self.assertAlmostEqual(r.value,1/(exact_degree+1),places=10)

    def test_reversed_and_zero_intervals(self):
        for method in ('Hình thang','Trung điểm','Simpson','Newton-Cotes'):
            with self.subTest(method=method):
                self.assertAlmostEqual(ca.integrate(method,'2*x+1',2,-1,12).value,-6)
                self.assertEqual(ca.integrate(method,'x^2',2,2,12).value,0)
        self.assertAlmostEqual(ca.double_integral('x+y',1,0,0,2).value,-3)
        self.assertEqual(ca.double_integral('x+y',0,0,0,2).value,0)

    def test_small_interval_large_function_has_finite_integral(self):
        r=ca.integrate('Hình thang','1e308',0,1e-308,2)
        self.assertAlmostEqual(r.value,1,places=12)

    def test_theoretical_bound_does_not_underflow_before_multiplication(self):
        # L^5*M4/180 = 1e-50/180, không phải 0.
        r=ca.integrate('Simpson','0',0,1e-70,2,bound=1e300)
        line=next(s for s in r.lines if s.startswith('Cận lý thuyết'))
        value=float(line.split('≤')[1].split(', nếu')[0])
        self.assertGreater(value,0)
        self.assertAlmostEqual(value/(1e-50/(16*180)),1,places=10)

    def test_auto_bound_and_runge_reach_requested_accuracy(self):
        for method in ('Hình thang','Trung điểm','Simpson','Newton-Cotes'):
            r=ca.integrate(method,'exp(x)',0,1,12,4,tolerance=1e-7)
            self.assertLess(abs(r.value-math.expm1(1)),2e-7)
        r=ca.integrate('Simpson','x^4',0,1,tolerance=1e-8,bound=24)
        self.assertLessEqual(abs(r.value-.2),1e-8)

    def test_table_function_and_midpoint_weights_agree(self):
        for method in ('Hình thang','Trung điểm','Simpson','Newton-Cotes'):
            x=(np.arange(12)+.5)/12 if method=='Trung điểm' else np.linspace(0,1,13)
            table=ca.integrate(method,a=0,b=1,x=x,y=x**4,n=12)
            direct=ca.integrate(method,'x^4',0,1,12)
            self.assertAlmostEqual(table.value,direct.value,places=14)
        with self.assertRaises(ValueError):ca.integrate(x=[0,1],y=None)

    def test_double_integral_has_correct_tensor_weights(self):
        self.assertAlmostEqual(ca.double_integral('x^3*y^3',-1,2,0,3,n=4).value,
                               ((2**4-1)/4)*(3**4/4),places=10)
        with self.assertRaises(ValueError):ca.double_integral('x',0,1,0,1,n=3)

    def test_adaptive_runge_checks_a_non_nested_grid_for_oscillatory_alias(self):
        r=ca.integrate('Simpson','cos(40*pi*x)',0,1,10,tolerance=1e-8)
        self.assertLess(abs(r.value),1e-8)
        r=ca.double_integral('cos(40*pi*x)',0,1,0,1,n=10,tolerance=1e-8)
        self.assertLess(abs(r.value),1e-8)


class OdeAudit(unittest.TestCase):
    def test_rk_orders_on_nonautonomous_forcing(self):
        for method,order in [('Euler hiện',1),('RK2',2),('RK3',3),('RK4',4)]:
            errors=[]
            for h in (.1,.05):
                _,y=de.ode(['x+sin(t)'],[1],0,1,h,method).value
                exact=1.5*math.e-(math.sin(1)+math.cos(1))/2
                errors.append(abs(y[-1,0]-exact))
            self.assertGreater(errors[0]/errors[1],2**(order-.3),method)

    def test_implicit_stiff_steps_match_the_discrete_equation(self):
        _,y=de.ode(['-100*x'],[1],0,1,.1,'Euler ẩn').value
        assert_allclose(y[:,0],(1/11)**np.arange(11),rtol=1e-7,atol=1e-14)
        _,y=de.ode(['-100*x'],[1],0,1,.1,'Hình thang').value
        assert_allclose(y[:,0],(-2/3)**np.arange(11),rtol=1e-10,atol=1e-12)

    def test_singular_implicit_step_is_not_reported_as_converged(self):
        with self.assertRaises(ValueError):de.ode(['x'],[1e-12],0,1,1,'Euler ẩn',tolerance=1e-14)

    def test_unused_state_variables_and_unrepresentable_time_grid_rejected(self):
        with self.assertRaises(ValueError):de.ode(['y'],[1],0,1,.1)
        with self.assertRaises(ValueError):de.grid(1e16,1e16+4,1)

    def test_rk3_special_families_satisfy_order_conditions(self):
        for c2,c3 in ((2/3,2/3),(2/3,0),(.5,1),(1/3,2/3),(.5,.75)):
            with self.subTest(c2=c2,c3=c3):
                a,b,c=de.custom_rk3(c2,c3)
                assert_allclose([b@c,b@(c*c),b@a@c],[.5,1/3,1/6],atol=1e-13)
                _,y=de.ode(['x'],[1],0,1,.025,'RK tùy chỉnh',custom=(a,b,c)).value
                self.assertLess(abs(y[-1,0]-math.e),2e-6)

    def test_all_rk2_parameterizations_represent_the_same_method(self):
        tables=[de.custom_rk2(name,value) for name,value in
                [('alpha2',2/3),('beta11',2/3),('r1',.25),('r2',.75)]]
        for table in tables[1:]:
            for a,b in zip(table,tables[0]):assert_allclose(a,b)

    def test_adams_predictor_corrector_exact_polynomial_forcing(self):
        for order in range(2,6):
            for degree in range(order):
                self.assertAlmostEqual(sum(c*(-j)**degree for j,c in enumerate(de.AB[order])),1/(degree+1),places=12)
            for degree in range(order+1):
                self.assertAlmostEqual(sum(c*(1-j)**degree for j,c in enumerate(de.AM[order])),1/(degree+1),places=12)
            t,y=de.ode(['t^2'],[0],0,1,.1,f'Adams {order}').value
            assert_allclose(y[:,0],t**3/3,atol=2e-13)

    def test_three_component_system_and_custom_table(self):
        t,y=de.ode(['y','-x','t'],[1,0,2],0,1,.025).value
        assert_allclose(y[-1],[math.cos(1),-math.sin(1),2.5],atol=4e-9)
        with self.assertRaises(ValueError):de.tableau([[0,0],[1,1]],[.5,.5],[0,2])


class BoundaryAndEigenAudit(unittest.TestCase):
    def test_variable_p_robin_boundary_is_second_order(self):
        errors=[]
        for h in (.1,.05,.025):
            x,y=de.bvp('1+x','0','-2-4*x',0,1,h,(1,0,0),(1,1,3)).value
            errors.append(float(np.max(abs(y-x*x))))
        self.assertGreater(errors[0]/errors[1],3.8,errors)
        self.assertGreater(errors[1]/errors[2],3.8,errors)

    def test_dirichlet_source_and_robin_scaling(self):
        args=('1','2','2*x^2-2',0,1,.025)
        x,y=de.bvp(*args,(0,1,0),(0,1,1)).value
        assert_allclose(y,x*x,atol=1e-11)
        _,z=de.bvp(*args,(0,1e-100,0),(0,1e100,1e100)).value
        assert_allclose(z,y,atol=1e-12)

    def test_eigen_discrete_spectrum_under_extreme_weights(self):
        h=.1;k=np.arange(1,4)
        expected=4/h**2*np.sin(k*np.pi*h/2)**2
        for weight in (1e200,1e-200,1.):
            with self.subTest(weight=weight):
                eig,_,_=de.eigen('1','0',str(weight),0,1,h,3).value
                assert_allclose(eig*weight,expected,rtol=2e-12,atol=0)

    def test_eigen_generalized_residual_and_orthogonality(self):
        values,x,v=de.eigen('1+x','2','1+x^2',0,1,.05,3).value
        h=x[1]-x[0];weights=1+x[1:-1]**2;faces=1+(x[:-1]+x[1:])/2
        diag=(faces[:-1]+faces[1:])/h**2+2;off=-faces[1:-1]/h**2
        a=np.diag(diag)+np.diag(off,1)+np.diag(off,-1)
        for j,value in enumerate(values):
            residual=a@v[1:-1,j]-value*weights*v[1:-1,j]
            self.assertLess(np.max(abs(residual)),1e-9)
        g=v[1:-1].T@(weights[:,None]*v[1:-1])
        assert_allclose(g-np.diag(np.diag(g)),0,atol=1e-11)

    def test_single_interior_node_and_invalid_weight(self):
        eig,x,v=de.eigen('1','0','1',0,1,.5,1).value
        assert_allclose(eig,[8]);self.assertEqual(v.shape,(3,1))
        with self.assertRaises(ValueError):de.eigen('1','0','0',0,1,.1,1)


class SplinesAndInverseAudit(unittest.TestCase):
    def test_spline_degree_one_through_four_reproduce_polynomials(self):
        x=np.array([0.,.3,1.,1.4,2.])
        for degree in range(1,5):
            p=Polynomial(np.arange(1.,degree+2))
            derivatives=[p.deriv(k)(0) for k in range(1,degree)]
            pp=ch1.spline_degree(x,p(x),degree,derivatives).value
            assert_allclose(pp(np.linspace(0,2,21)),p(np.linspace(0,2,21)),atol=2e-10)
            right=ch1.spline_boundary(x,p(x),degree,[('phải',k,p.deriv(k)(x[-1])) for k in range(1,degree)]).value
            assert_allclose(right(np.linspace(0,2,21)),p(np.linspace(0,2,21)),atol=2e-10)
            for knot in x[1:-1]:
                for order in range(degree):
                    self.assertAlmostEqual(ch1.spline_at(pp,knot,order),p.deriv(order)(knot),places=8)

    def test_spline_two_endpoint_system_handles_small_units(self):
        x=np.array([0.,1.,2.])*1e-8
        pp=ch1.spline_boundary(x,x,4,[('trái',1,1),('trái',2,0),('trái',3,0)]).value
        assert_allclose(pp(x),x,rtol=1e-10,atol=1e-20)
        self.assertAlmostEqual(pp(1.5e-8,nu=1),1,places=9)

    def test_quadratic_splines_both_boundaries_and_cubic_variants(self):
        x=np.array([0.,.2,.7,1.,2.]);y=x*x
        for kind in ('Bậc hai — biên trái','Bậc hai — biên phải','Cubic kẹp','Cubic đạo hàm bậc hai','Cubic not-a-knot'):
            left,right=(2,2) if kind=='Cubic đạo hàm bậc hai' else (0,4)
            pp=ip.spline(x,y,kind,left,right).value
            assert_allclose(pp(np.linspace(0,2,25)),np.linspace(0,2,25)**2,atol=1e-12)
        x=np.linspace(0,2*math.pi,9);y=np.sin(x);y[-1]=y[0]
        pp=ip.spline(x,y,'Cubic tuần hoàn').value
        self.assertAlmostEqual(pp(0,1),pp(2*math.pi,1),places=13)
        self.assertAlmostEqual(pp(0,2),pp(2*math.pi,2),places=13)

    def test_inverse_increasing_decreasing_and_backward(self):
        x=np.linspace(1,1.4,5)
        for sign in (1,-1):
            for backward in (False,True):
                r=ip.inverse(x,sign*x*x,sign*1.3225,backward)
                self.assertAlmostEqual(r.value,1.15,places=8)
                r=ch1.inverse_controlled(x,sign*x*x,sign*1.3225,(1.1,1.2),backward)
                self.assertAlmostEqual(r.value,1.15,places=8)
        r=ch1.inverse_by_swap(x,-x, -1.15)
        self.assertAlmostEqual(r.value,1.15,places=12)

    def test_scan_level_and_missing_data_without_assuming_monotonic_function(self):
        roots,intervals,flat=ch1.scan_level([0,1,2,3],[0,1,1,0],1).value
        self.assertEqual(roots,[1,2]);self.assertEqual(intervals,[]);self.assertEqual(flat,[(1,2)])
        r=ch1.fill_missing([3,1,0,2],[9,None,0,4],2)
        assert_allclose(r.value,[9,1,0,4])

    def test_large_inverse_maps_and_small_certified_interval(self):
        x=np.arange(201.)
        for backward in (False,True):
            self.assertAlmostEqual(ip.inverse(x,2*x+1,1.5,backward).value,.25)
            self.assertAlmostEqual(ch1.inverse_controlled(x,2*x+1,1.5,(0,1),backward).value,.25)
        x=np.linspace(1,1.4,5)
        self.assertAlmostEqual(ch1.inverse_controlled(x,x*x,1.3225,(1.149,1.151)).value,1.15,places=9)

    def test_tiny_derivative_jump_is_not_hidden_by_absolute_tolerance(self):
        pp=ip.spline([0,1,2],[0,0,1e-12],'Tuyến tính').value
        with self.assertRaises(ValueError):ch1.spline_at(pp,1,1)
        with self.assertRaises(ValueError):ch1.fill_missing([0,1,2,3],[0,1e-12,5e-12,None],1)

    def test_large_splines_do_not_build_dense_presentation_systems(self):
        x=np.linspace(0,1,1001)
        with patch('pps.math_view.np.linalg.solve',side_effect=AssertionError('Không dựng hệ đặc để trình bày bảng lớn')):
            report=ip.spline(x,x,'Cubic tự nhiên')
        assert_allclose(report.value([.25,.75]),[.25,.75],atol=1e-12)
        self.assertLess(max(len(headers) for _,headers,_ in report.tables),10)
        pp=ch1.spline_boundary(x,x,2,[('trái',1,1)]).value
        assert_allclose(pp([.25,.75]),[.25,.75],atol=1e-12)
        pp=ch1.spline_degree(x,x,4,[1,0,0]).value
        assert_allclose(pp([.25,.75]),[.25,.75],atol=1e-12)
        for side in ('trái','phải'):
            pp=ch1.spline_boundary(x,x,4,[(side,1,1),(side,2,0),(side,3,0)]).value
            assert_allclose(pp([.25,.75]),[.25,.75],atol=1e-12)


class RegressionAudit(unittest.TestCase):
    def test_weight_scale_does_not_change_fit_or_destroy_normal_matrix(self):
        x=np.array([0.,1.,2.]);y=np.array([1.,2.,2.])
        for factor in (1e-300,1.,1e300):
            r=ip.least_squares(x,y,weights=np.array([1,2,1])*factor)
            assert_allclose(r.value([0,1]),[1.25,1.75],atol=1e-12)
            for name,_,rows in r.tables:
                if 'ΘᵀWΘ' in name:
                    self.assertTrue(all(np.isfinite(float(v)) for row in rows for v in row))

    def test_normal_matrix_can_exceed_float_without_losing_fit(self):
        x=np.arange(3.)*1e10
        r=ip.least_squares(x,[1,2,2],weights=[1e300,2e300,1e300])
        assert_allclose(r.value([0,1e10]),[1.25,1.75],atol=1e-12)
        rows=next(rows for name,_,rows in r.tables if 'ΘᵀWΘ' in name)
        self.assertTrue(all(Decimal(str(v)).is_finite() for row in rows for v in row))
        self.assertGreater(Decimal(str(rows[1][1])),Decimal('1e308'))

    def test_affine_basis_at_large_offset_remains_distinct(self):
        x=1e12+np.arange(4.)
        r=ip.least_squares(x,[1,2,3,4])
        assert_allclose(r.value(x),[1,2,3,4],atol=1e-12)
        self.assertAlmostEqual(r.value(1e12+.5),1.5)

    def test_independent_weighted_solution_and_residual_orthogonality(self):
        x=np.array([-.7,0,.4,1,2]);y=np.array([2.,1.,2.,2.,5.]);w=np.array([1.,3.,2.,4.,1.])
        r=ip.least_squares(x,y,basis=['1','x','x^2'],weights=w)
        a=np.column_stack((np.ones(5),x,x*x))
        expected=np.linalg.solve(a.T@(w[:,None]*a),a.T@(w*y))
        assert_allclose(r.value(x),a@expected,atol=1e-12)
        assert_allclose(a.T@(w*(r.value(x)-y)),0,atol=1e-11)

    def test_transformed_exponential_predicts_when_base_exceeds_float(self):
        x=np.array([0.,.001,.002,.003]);y=2*np.exp(1000*x)
        r=ch1.transformed_fit(x,y,'A*b^x')
        assert_allclose(r.value(x),y,rtol=1e-12)

    def test_original_scale_prediction_preserves_array_shape(self):
        x=np.array([1.,2.,3.,4.]);y=2*x*x
        r=ch1.original_scale_fit(x,y,'Lũy thừa')
        z=np.array([[1.,2.],[3.,4.]])
        assert_allclose(r.value(z),2*z*z,rtol=1e-10)

    def test_original_scale_negative_and_zero_data(self):
        for y in (-2*np.exp(.3*np.arange(4)),np.zeros(4)):
            r=ch1.original_scale_fit(np.arange(4),y)
            assert_allclose(r.value(np.arange(4)),y,atol=1e-10)

    def test_original_scale_optimum_is_independent_of_y_units(self):
        x=np.arange(4.);y=np.array([2.,3.1,5.2,8.9])
        base=ch1.original_scale_fit(x,y).value(x)
        tiny=ch1.original_scale_fit(x,y*1e-200).value(x)/1e-200
        assert_allclose(tiny,base,rtol=1e-8,atol=1e-8)

    def test_log_and_reciprocal_model_domains_and_dependent_basis(self):
        with self.assertRaises(ValueError):ip.least_squares([0,1,2],[1,2,3],basis=['1','2'])
        with self.assertRaises(ValueError):ch1.transformed_fit([1,2],[0,1],'1/(a+b*x)')
        with self.assertRaises(ValueError):ch1.transformed_fit([0,1],[1,2],'a+b*ln(x)')
        r=ch1.transformed_fit([1,2,3],1/(2+3*np.arange(1,4)),'1/(a+b*x)')
        assert_allclose(r.value([1,2,3]),1/(2+3*np.arange(1,4)),atol=1e-13)


class LinearSolversAudit(unittest.TestCase):
    def test_mixed_row_units_and_intermediate_cancellation(self):
        assert_allclose(ch1.gaussian_elimination([[1e-100,0],[0,1e100]],[2e-100,3e100]).value,[2,3])
        assert_allclose(ch1.thomas([0,0],[1e-100,1e100],[0,0],[2e-100,3e100]).value,[2,3])
        r=ch1.gaussian_elimination([[1e308,-1e308],[0,1]],[0,2])
        assert_allclose(r.value,[2,2])

    def test_both_solvers_recover_constructed_solutions(self):
        rng=np.random.default_rng(142)
        for n in range(1,12):
            lower=np.r_[0.,rng.uniform(-1,1,n-1)]
            upper=np.r_[rng.uniform(-1,1,n-1),0.]
            diagonal=2+abs(lower)+abs(upper)
            a=np.diag(diagonal)+np.diag(lower[1:],-1)+np.diag(upper[:-1],1)
            exact=rng.normal(size=n);rhs=a@exact
            assert_allclose(ch1.thomas(lower,diagonal,upper,rhs).value,exact,atol=1e-12)
            assert_allclose(ch1.gaussian_elimination(a,rhs).value,exact,atol=1e-12)


if __name__=='__main__':unittest.main()
