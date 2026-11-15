import unittest
import numpy as np
from check_nonconformal_flow_accuracy import metrics,difference

class AccuracyChecks(unittest.TestCase):
    def test_affine_gradient_on_oblique_tet(self):
        x=np.array([[.2,.3,.4],[1.2,.4,.7],[.1,1.8,.5],[.5,.2,2.1]])
        a=np.array([[2.,3.,4.],[-3.,-1.,2.],[1.,4.,.5]])
        u=np.column_stack([x@a,np.ones(4)])
        cells=np.array([[0,1,2,3]])
        m=metrics(x,u,cells)
        self.assertAlmostEqual(m['raw_cell_divergence_l1'],m['volume']*abs(np.trace(a)))
        self.assertAlmostEqual(m['raw_cell_divergence_l2'],np.sqrt(m['volume'])*abs(np.trace(a)))
        # Reversing orientation must not change an integral.
        self.assertEqual(m,metrics(x,u,cells[:,[0,2,1,3]]))

    def test_constant_difference_mass_integral(self):
        x=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        c=np.array([[0,1,2,3]])
        d=np.tile([1.,2.,2.,4.],(4,1)); z=np.zeros((4,4))
        result=difference((x,d,c),(x,z,c))
        self.assertAlmostEqual(result['velocity_interpolant_l2'],3/np.sqrt(6))
        self.assertAlmostEqual(result['pressure_interpolant_l2'],4/np.sqrt(6))
        with self.assertRaises(ValueError): difference((x,d,c),(np.concatenate([x,x]),np.concatenate([z,z]),c))

if __name__=='__main__': unittest.main()
