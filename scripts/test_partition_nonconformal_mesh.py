import unittest
import tempfile
from pathlib import Path
import numpy as np
from test_prepare_nonconformal_mesh import fixture,prepare
from partition_nonconformal_mesh import partitions,read_prepared,write_partitions

class PartitionChecks(unittest.TestCase):
    def test_incident_cells_and_separate_interface_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp)/'test.exo'; fixture(f)
            arrays=prepare.prepare(f,flow_boundaries=True)
            source=Path(tmp)/'input.mnc'; prepare.write_mesh(source,arrays)
            a=read_prepared(source)
            for nr in (1,2,4,8):
                parts=list(partitions(a,nr))
                for rank,part in enumerate(parts):
                    nodes,elements,x,cells,k,bc,side,owners,tags,outer,peers=part
                    own=set(range(len(a[0])*rank//nr,len(a[0])*(rank+1)//nr))
                    wanted={e for e,cell in enumerate(a[1]) if own.intersection(cell)}
                    self.assertEqual(set(elements),wanted)
                    self.assertTrue(own.issubset(nodes))
                    np.testing.assert_array_equal(nodes[cells],a[1][elements])
                    np.testing.assert_array_equal(x,a[0][nodes])
                    for tag in tags:
                        self.assertIn(tuple([elements[tag[0]],*tag[1:]]),map(tuple,a[5]))
                    for peer in peers: self.assertIn(rank,parts[peer][-1])
                # Coincident coordinates retain independent field identities.
                self.assertEqual(len(set(np.concatenate([p[0] for p in parts]))),8)
            report=write_partitions(source,Path(tmp)/'parts',4)
            self.assertEqual(report['ranks'],4)
            with self.assertRaises(FileExistsError): write_partitions(source,Path(tmp)/'parts',4)
            source.write_bytes(source.read_bytes()[:-1])
            with self.assertRaises(ValueError): read_prepared(source)

if __name__=='__main__': unittest.main()
