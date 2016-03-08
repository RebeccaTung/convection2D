''' Unittests for RedbackContinuation functions '''

import os, sys, logging, unittest, multiprocessing, csv, math, difflib

from CheckMooseOutput import checkMooseOutput, MooseException
from MooseInputFileRW import MooseInputFileRW
from RedbackContinuation import runContinuation
from Utils import NullHandler, getLogger

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

class TestStringMethods(unittest.TestCase):

  def setUp(self):
    ''' Initial set up for all tests '''
    self.logger = getLogger('log.txt', level=logging.INFO)
    self.error_filename = os.path.join('tests', 'error.txt')
    self.nb_procs = multiprocessing.cpu_count()
  
  def tearDown(self):
    ''' Clean up function called after each test '''
    return

  def test_checkMooseOutputOK1(self):
    self.assertEqual(checkMooseOutput(STRING1, self.error_filename, self.logger), 0)

  def test_checkMooseOutputNoConvergence1(self):
    self.assertRaises(MooseException, checkMooseOutput, STRING2, self.error_filename, self.logger)
  
  def test_mooseInputFileRW(self):
    handler = MooseInputFileRW()
    #handler.logger.setLevel(logging.DEBUG)
    filename_in = os.path.join('tests', 'MooseInputFilesRW', 'bench1_a.i')
    filename_out = os.path.join('tests', 'MooseInputFilesRW', 'bench1_a_recreated.i')
    data = handler.read(filename_in)
    handler.write(data, filename_out)
    with open(filename_in, 'r') as f_1:
      file_content_1 = f_1.read()
    with open(filename_out, 'r') as f_2:
      file_content_2 = f_2.read()
    self.assertMultiLineEqual(file_content_1, file_content_2, 
                              msg='Problem with MooseInputFileRW')
  
  def test_SCurve_1D(self):
    ''' Test that 1D simulation with 1 parameter produces the correct S-curve '''
    # run simulation
    parameters = {
      'lambda_initial_1':1e-8,
      'lambda_initial_2':2e-8,
      'ds_initial':5e-2,
      's_max':0.5,
      # Rescaling factor
      'rescaling_factor':4.5399929762e-5, # to multiply continuation parameter
      # Numerical parameters
      'exec_loc':'~/projects/redback/redback-opt',
      'nb_threads':self.nb_procs,
      'input_file':'tests/S_curve_1D/bench1_a.i',
      'running_dir':'tests/running_tmp',
      'result_curve_csv':'tests/S_curve_1D/S_curve.csv',
      'error_filename':'tests/error_output.txt',
      'plot_s_curve':False,
      'ref_s_curve':'tests/S_curve_1D/ref.csv'
    }
    results = runContinuation(parameters, self.logger)
    # check results
    self._compareCsvFiles(parameters['result_curve_csv'], 'tests/S_curve_1D/S_curve_gold.csv', 
                          relative_error=1e-5)
  
  def _parseResultCsvFile(self, filename):
    ''' Parse csv file storing S-curve results and return data. Raises Exception if problem found.
        @param[in] filename - string, filename to parse
        @return data - dictionary of data, keys=column header, values=list of data 
    '''
    data = {}
    with open(filename, 'rb') as csvfile:
      csvreader = csv.reader(csvfile)
      line_i = 0 # line index
      for row in csvreader:
        if line_i == 0:
          # Headers
          headers = [elt.strip().lower() for elt in row]
          nb_col = len(headers)
          for label in headers:
            data[label] = []
          line_i += 1
          continue
        # data line
        assert len(row) == nb_col
        for i, elt in enumerate(row):
          data[headers[i]].append(float(elt))
    return data
  
  def _compareCsvFiles(self, filename1, filename2, abs_zero=1e-11, relative_error=1e-6):
    ''' Compare the contents of 2 csv files and raise Exception if contents differ by more than
        acceptable tolerance (for non zero entries), or by more than absolute tolerance for zero entries.
        @param[in] filename1 - string, filename of first file (simulation results)
        @param[in] filename2 - string, filename of second file (reference results)
        @param[in] abs_zero - float, threshold to determine when a value is ~zero
        @param[in] relative_error - float, precision required to match ratios of non zero values
    '''
    for filename in [filename1, filename2]:
      if not os.path.isfile(filename):
        raise Exception, 'File "{0}" not found'.format(filename)
    data1 = self._parseResultCsvFile(filename1)
    data2 = self._parseResultCsvFile(filename2)
    # Check headers
    headers1 = data1.keys()
    headers2 = data2.keys()
    nb_col = len(headers2)
    assert len(headers1) == nb_col, 'Different number of columns in files "{0}" and "{1}"'.\
      format(filename1, filename2)
    # Check data
    if not len(data1):
      # nothing to compare
      return
    nb_rows = len(data2[headers1[0]])
    assert len(data1[headers1[0]]) == nb_rows, 'File "{0}" contains {1} data lines instead of {2}'.\
      format(filename1, len(data1[headers1[0]]), nb_rows)
    for col_i in range(nb_col):
      header = headers2[col_i]
      for row_i in range(nb_rows):
        if math.fabs(data2[header][row_i]) < abs_zero:
          assert math.fabs(data1[header][row_i]) < abs_zero, 'Data for column {0} row {1} differ: {2} != {3} (ref)'.\
            format(col_i, row_i, data1[header][row_i], data2[header][row_i])
        else:
          measured_error = math.fabs((data1[header][row_i] - data2[header][row_i]) / data2[header][row_i]) 
          assert measured_error < relative_error, \
            'Data for column "{0}", row {1} differ: {2} != {3} (ref), relative error = {4} > {5}'.\
            format(header, row_i, data1[header][row_i], data2[header][row_i],
                   measured_error, relative_error)
    return
  
  def assertMultiLineEqual(self, first, second, msg=None):
    ''' Assert that two multi-line strings are equal.
        If they aren't, show a nice diff.
        See http://stackoverflow.com/questions/3942820/how-to-do-unit-testing-of-functions-writing-files-using-python-unittest
    '''
    self.assertTrue(isinstance(first, str), 'First argument is not a string')
    self.assertTrue(isinstance(second, str), 'Second argument is not a string')
    if first.strip() != second.strip():
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