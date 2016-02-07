''' Unittests for RedbackContinuation functions '''

import os, sys, logging, unittest

from CheckMooseOutput import checkMooseOutput, MooseException
from MooseInputFileRW import MooseInputFileRW

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

from Utils import NullHandler, getLogger

class TestStringMethods(unittest.TestCase):

  def setUp(self):
    ''' Initial set up for all tests '''
    self.logger = getLogger('log.txt', level=logging.INFO)
  
  def tearDown(self):
    ''' Clean up function called after each test '''
    return

  def test_checkMooseOutputOK1(self):
    self.assertEqual(checkMooseOutput(STRING1, self.logger), 0)

  def test_checkMooseOutputNoConvergence1(self):
    self.assertRaises(MooseException, checkMooseOutput, STRING2, self.logger)
  
  def test_mooseInputFileRW(self):
    handler = MooseInputFileRW()
    #handler.logger.setLevel(logging.DEBUG)
    filename_in = os.path.join('tests', 'bench1_a.i')
    filename_out = os.path.join('tests', 'bench1_a_recreated.i')
    data = handler.read(filename_in)
    handler.write(data, filename_out)
    with open(filename_in, 'r') as f_1:
      file_content_1 = f_1.read()
    with open(filename_out, 'r') as f_2:
      file_content_2 = f_2.read()
    self.assertMultiLineEqual(file_content_1, file_content_2, 
                              msg='Problem with MooseInputFileRW')
  
  def assertMultiLineEqual(self, first, second, msg=None):
    ''' Assert that two multi-line strings are equal.
        If they aren't, show a nice diff.
        See http://stackoverflow.com/questions/3942820/how-to-do-unit-testing-of-functions-writing-files-using-python-unittest
    '''
    self.assertTrue(isinstance(first, str), 'First argument is not a string')
    self.assertTrue(isinstance(second, str), 'Second argument is not a string')
    if first != second:
      message = ''.join(difflib.ndiff(first.splitlines(True),
                                      second.splitlines(True)))
      if msg:
        message += ' : ' + msg
      self.fail('Multi-line strings are unequal:\n' + message)


if __name__ == "__main__":
  #unittest.main()
  suite = unittest.TestLoader().loadTestsFromTestCase(TestStringMethods)
  unittest.TextTestRunner(verbosity=2).run(suite)

  print 'Finished'