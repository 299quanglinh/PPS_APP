import unittest
import math
from pps.gui import dispatch, fields, TOOLS
from pps.inverse_workflow import iteration_workflow
from pps.math_view import math_png

class InverseWorkflowTests(unittest.TestCase):
    def test_user_data_and_independent_quadratic(self):
        # Mốc cục bộ từ bảng 20251, kèm một mốc ngoài bộ được chọn.
        r=iteration_workflow([2.38,2.495,2.61,2.725],[-1.5819,-1.5127,-1.3977,-1.3814],-1.43)
        self.assertEqual(len(r.value),1)
        result=r.value[0]
        self.assertEqual(result['status'],'DAT')
        expected=2.495+.115*(.16435-math.sqrt(.16435**2-4*.04935*.0827))/(2*.04935)
        self.assertLessEqual(abs(result['value']-expected),result['bound'])
        self.assertEqual(result['steps'],8)
        self.assertGreaterEqual(result['displayed_bound'],result['bound'])
        self.assertLessEqual(abs(float(f"{result['value']:.12f}")-expected),result['displayed_bound'])
        self.assertLess(result['residual'],1e-8)
        content=r.text()
        for heading in ['5. Đánh giá sai số','6. Kết luận nhánh','Chặn hậu nghiệm B_j','Chặn dùng kiểm tra E_j']:
            self.assertIn(heading,content)


    def test_pasted_pdf_example_auto_interval_and_manual_37(self):
        x=[1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8]
        y=[4.02,3.98,4.23,3.67,2.99,1.24,.87,2.08,4.54]
        r=iteration_workflow(x,y,3,tolerance=1e-6)
        self.assertEqual([a['status'] for a in r.value],['DAT','THIEU_MOC'])
        self.assertAlmostEqual(r.value[0]['value'],1.998348,places=5)
        manual=iteration_workflow(x,y,3,tolerance=1e-6,alpha=.985,beta=.997)
        self.assertEqual(manual.value[0]['status'],'DAT')

    def test_three_branches(self):
        y=[3.1130627,3.1318007,3.144226,3.1512136,3.1536928,3.1526381,
           3.149059,3.1439906,3.1384827,3.1335894,3.1303588,3.1298225,
           3.1329848,3.1408133,3.1542291,3.1740972,3.2012183,3.2363205,3.2800524]
        x=[2.4+i*.1 for i in range(len(y))]
        r=iteration_workflow(x,y,3.15,degree=5)
        self.assertEqual(len(r.value),3)
        self.assertTrue(all(v['status']=='DAT' for v in r.value))
        for a,b in zip(r.value,[2.67628868,2.97825566,3.77295959]):self.assertAlmostEqual(a['value'],b,places=7)

    def test_backward_and_node(self):
        r=iteration_workflow([1,2,3,4],[1,4,9,16],7,True)
        self.assertEqual(r.value[0]['status'],'DAT')
        self.assertAlmostEqual(r.value[0]['value'],math.sqrt(7),places=7)
        r=iteration_workflow([0,1,2],[0,1,4],1)
        self.assertEqual(r.value[0]['status'],'TRUNG_MOC')

    def test_failures_keep_branch_status(self):
        r=iteration_workflow([0,1,2],[0,1,4],.5,True)
        self.assertEqual(r.value[0]['status'],'THIEU_MOC')
        r=iteration_workflow([2,3,4],[4,9,16],7,max_steps=1)
        self.assertEqual(r.value[0]['status'],'CHUA_DAT')
        r=iteration_workflow([2,3,4],[4,9,16],7,alpha=0,beta=.1)
        self.assertEqual(r.value[0]['status'],'CHUA_CHUNG_NHAN')
        with self.assertRaises(ValueError):iteration_workflow([0,1,3],[0,1,9],.5)
        self.assertEqual(iteration_workflow([0,1,2],[0,1,4],10).value,[])

    def test_dispatch_and_all_formulas_render(self):
        self.assertEqual(TOOLS['Nội suy ngược'][0],['Lặp tiến','Lặp lùi'])
        self.assertIn('degree',[f[0] for f in fields('Nội suy ngược','Lặp tiến')])
        for method in ['Lặp tiến','Lặp lùi']:
            r=dispatch('Nội suy ngược',method,{'target':'7','tol':'1e-6','degree':'2'},'1 2 3 4','1 4 9 16')
            self.assertEqual(r.value[0]['status'],'DAT')
            for kind,value in r.presentation['blocks']:
                if kind=='formula':self.assertTrue(math_png(value).startswith(b'\x89PNG'))
            self.assertIn('Đầu vào',r.text());self.assertIn('Kiểm tra hội tụ',r.text())

    def test_method_specific_notation_and_backward_table(self):
        for backward in (False,True):
            r=iteration_workflow([1,2,3,4],[1,4,9,16],7,backward)
            formulas='\n'.join(v for kind,v in r.presentation['blocks'] if kind=='formula')
            for generic in ['x_g','y_g','c_1','F(']:self.assertNotIn(generic,formulas)
            if backward:
                self.assertIn(r'x_{k+1}',formulas)
                self.assertIn(r'\nabla y_{k+1}',formulas)
                self.assertIn(r'\psi(t)',formulas)
                _,headers,rows=next(t for t in r.tables if 'bảng sai phân lùi' in t[0])
                self.assertIn('∇^1 y_i',headers)
                self.assertEqual(rows[-1][1:],['9.000000000000','5.000000000000','2.000000000000'])
                self.assertEqual(rows[0][2:],["",""])
            else:
                self.assertIn('x_k',formulas)
                self.assertIn(r'\Delta y_k',formulas)
                self.assertIn(r'\varphi(t)',formulas)
