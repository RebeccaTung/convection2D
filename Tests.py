''' Unittests for RedbackContinuation functions '''

import os, sys, logging, unittest

from CheckMooseOutput import checkMooseOutput, MooseException

STRING1 = '''
Postprocessor Values:
+----------------+----------------+----------------+----------------+----------------+----------------+----------------+----------------+
| time           | max_temp       | temp_pt_0      | temp_pt_1      | temp_pt_2      | temp_pt_3      | temp_pt_4      | temp_pt_5      |
+----------------+----------------+----------------+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
+----------------+----------------+----------------+----------------+----------------+----------------+----------------+----------------+

 0 Nonlinear |R| = [32m6.000000e-03[39m
      0 Linear |R| = [32m6.000000e-03[39m
      1 Linear |R| = [32m2.268915e-12[39m
 1 Nonlinear |R| = [32m4.310185e-06[39m
      0 Linear |R| = [32m4.310185e-06[39m
      1 Linear |R| = [32m3.157957e-15[39m
 2 Nonlinear |R| = [32m2.591767e-12[39m

Postprocessor Values:
+----------------+----------------+----------------+----------------+----------------+----------------+----------------+----------------+
| time           | max_temp       | temp_pt_0      | temp_pt_1      | temp_pt_2      | temp_pt_3      | temp_pt_4      | temp_pt_5      |
+----------------+----------------+----------------+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+00 |   5.219478e-03 |   0.000000e+00 |   1.868683e-03 |   3.330103e-03 |   4.378302e-03 |   5.008965e-03 |   5.219478e-03 |
+----------------+----------------+----------------+----------------+----------------+----------------+----------------+----------------+


'''
STRING2 = '''
      8 Linear |R| = [33m2.326106e+00[39m
      9 Linear |R| = [32m1.455026e-08[39m
Aborting as solve did not converge

'''

class NullHandler(logging.Handler):
  def emit(self, record):
    pass

class TestStringMethods(unittest.TestCase):

  def setUp(self):
    ''' Initial set up for all tests '''
    self.logger = logging.getLogger('tests')
    self.logger.setLevel(logging.INFO)
    h = NullHandler()
    self.logger.addHandler(h)
  
  def tearDown(self):
    ''' Clean up function called after each test '''
    return

  def test_checkMooseOutputOK1(self):
    self.assertEqual(checkMooseOutput(STRING1, self.logger), 0)

  def test_checkMooseOutputNoConvergence1(self):
    #self.assertEqual(checkMooseOutput(STRING2, self.logger), 1)
    self.assertRaises(MooseException, checkMooseOutput, STRING2, self.logger)

if __name__ == "__main__":
  #unittest.main()
  suite = unittest.TestLoader().loadTestsFromTestCase(TestStringMethods)
  unittest.TextTestRunner(verbosity=2).run(suite)

  print 'Finished'