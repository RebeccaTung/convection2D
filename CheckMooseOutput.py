''' Function to check if Moose simulation converged or not based on the printed output '''

import os, sys, logging

class MooseException(Exception):
  pass

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
  ### Known error case 1: 'Aborting as solve did not converge'
  # Find last non-empty line
  index = len(lines) - 1
  last_line = ''
  while index > -1:
    if lines[index].strip() is not '':
      last_line = lines[index]
      break
    index -= 1
  if last_line is not '':
    if last_line.lower().strip() == 'aborting as solve did not converge':
      error = True
  
  if error:
    for line in lines:
      print line
    logger.error('Simulation failed\n' + '\n'.join(lines[-5:-1]))
    raise MooseException('Simulation did not converge')
  return 0

if __name__ == "__main__":
  stdout = 'hello world'
  checkMooseOutput(stdout, logging.getLogger('test'))
  print 'Finished'