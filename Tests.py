''' Unittests for RedbackContinuation functions '''

import os, sys, logging, unittest

from CheckMooseOutput import checkMooseOutput

class TestStringMethods(unittest.TestCase):

  def setUp(self):
    ''' Initial set up for all tests '''
    self.logger = logging.getLogger('tests')
    self.logger.setLevel(logging.INFO)
  
  def tearDown(self):
    ''' Clean up function called after each test '''
    return

  def test_checkMooseOutput1(self):
    stdout = 'No error here...'
    self.assertEqual(checkMooseOutput(stdout, self.logger), 0)
    # assertRaises(Exception, checkMooseOutput, [stdout, self.logger])

if __name__ == "__main__":
  #unittest.main()
  suite = unittest.TestLoader().loadTestsFromTestCase(TestStringMethods)
  unittest.TextTestRunner(verbosity=2).run(suite)

  print 'Finished'