import unittest
from check_nonconformal_flow_fields import compare

class FieldChecks(unittest.TestCase):
    def test_difference_and_identity_are_not_hidden_by_rank_order(self):
        a={1:dict(x=0.,y=0.,z=0.,velocity_x=.1,velocity_y=0.,velocity_z=0.,pressure=.2)}
        self.assertTrue(compare(a,a)['passed'])
        b={1:dict(a[1],velocity_y=1e-4)}
        self.assertFalse(compare(a,b)['passed'])
        with self.assertRaises(ValueError):compare(a,{2:a[1]})
        with self.assertRaises(ValueError):compare(a,{1:dict(a[1],x=.1)})

if __name__=='__main__':unittest.main()
