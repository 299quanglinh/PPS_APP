import unittest
import numpy as np
from pps.interpolation import least_squares

class RepeatedObservations(unittest.TestCase):
    def test_replicates_contribute_to_objective(self):
        x=np.array([0.,1.,1.,2.]);y=np.array([0.,1.,3.,2.])
        r=least_squares(x,y,basis=['1','x'])
        expected=np.linalg.lstsq(np.column_stack([np.ones(4),x]),y,rcond=None)[0]
        np.testing.assert_allclose(r.value(x),expected[0]+expected[1]*x)
        self.assertEqual(len(r.tables[-1][2]),4)

if __name__=='__main__':unittest.main()
