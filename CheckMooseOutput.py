''' Function to check if Moose simulation converged or not based on the printed output '''

import os, sys, logging

class MooseException(Exception):
  pass

def checkMooseOutput(stdout, error_filename, logger):
  ''' Check string containing moose's stdout output from simulation.
      Raise Exception if string shows simulation failed.
      @param[in] stdout - string containing stdout from Moose simulation
      @param[in] logger - python logger instance  
  '''
  lines = stdout.split('\n')
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
    #for line in lines:
    #  print line
    # store full output to file
    root_dir, short_filename = os.path.split(error_filename)
    if not os.path.isdir(root_dir):
      os.makedirs(root_dir)
    if os.path.isfile(error_filename):
      os.remove(error_filename)
    with open(error_filename, 'wb') as error_file:
      error_file.writelines(lines)
    logger.debug('Wrote error file to "{0}"'.format(error_filename))
    raise MooseException('Simulation did not converge')
  else:
      # Print all output in debug mode
    if logger.level <= logging.DEBUG:
      for line in lines:
        print line
  
  return 0

if __name__ == "__main__":
  stdout = 'hello world'
  checkMooseOutput(stdout, logging.getLogger('test'))
  print 'Finished'