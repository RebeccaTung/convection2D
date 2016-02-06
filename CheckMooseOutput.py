''' Function to check if Moose simulation converged or not based on the printed output '''

import os, sys, logging

def checkMooseOutput(stdout, logger):
  ''' Check string containing moose's stdout output from simulation.
      Raise Exception if string shows simulation failed.
      @param[in] stdout - string containing stdout from Moose simulation
      @param[in] logger - python logger instance  
  '''
  lines = stdout.split('\n')
  # Print all output in debug mode
  if logger.level <= logging.DEBUG:
    for line in lines:
      print line
  # Check for errors
  error = False
  logger.warning('TODO: checkMooseOutput to be implemented and tested')
  #import pdb;pdb.set_trace()
  if error:
    for line in lines:
      print line
    logger.error('Simulation failed\n' + '\n'.join(lines[-5:-1]))
    return 1
  return 0

if __name__ == "__main__":
  stdout = 'hello world'
  checkMooseOutput(stdout, logging.getLogger('test'))
  print 'Finished'