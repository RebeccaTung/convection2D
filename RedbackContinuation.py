''' Script to run Redback continuation. 
  TODO: everything...
'''

import os, sys, random, logging, subprocess, shutil, math, csv
from os.path import expanduser

from Utils import getLogger
from CheckMooseOutput import checkMooseOutput
from PlotSCurve import plotSCurve

def createRedbackFilesRequired(parameters, logger):
  ''' Create Redback files required to run continuation 
      @param[in] parameters - dictionary of input parameters
      @param[in] logger - python logger instance
  '''
  # Create output csv file for S-curve results
  with open(parameters['result_curve_csv'], 'wb') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['step', 'lambda', 'norm'])
  # Create other input files
  logger.warning('TODO: createRedbackFilesRequired not implemented')

def checkAndCleanInputParameters(parameters, logger):
  ''' Check input parameters provided by user 
      @param[in] parameters - dictionary of input parameters
      @param[in] logger - python logger instance  
      @return: found_error - boolean, True if any error was found
  '''
  found_error = False
  # Check all parameters are provided
  if type(parameters) != dict:
    logger.error('"parameters" must be of type dictionary, not {0}!'.format(type(parameters)))
    return True
  params_real = ['lambda_initial_1', 'lambda_initial_2', 'ds_initial', 's_max']
  params_int = ['nb_threads']
  params_str = ['exec_loc', 'input_dir', 'running_dir', 'result_curve_csv', 'ref_s_curve']
  params_bool = ['plot_s_curve']
  all_keys = params_real + params_int + params_str + params_bool
  missing_param = False
  for param in all_keys:
    if not param in parameters.keys():
      logger.error('Input parameter "{0}" is missing!'.format(param))
      missing_param = True
  if missing_param:
    return True
  # Ensure proper type of numerical parameters
  for key in all_keys:
    if key not in parameters.keys():
      logger.error('Input parameter "{0}={1}" should be an integer!'.format(key, parameters[key]))
      found_error = True
  if found_error:
    return Ture
  for param in params_real:
    if not type(parameters[param]) in [float, int]:
      logger.error('Input parameter "{0}={1}" should be a number!'.format(param, parameters[param]))
      found_error = True
  for param in params_int:
    if not type(parameters[param]) == int:
      logger.error('Input parameter "{0}={1}" should be an integer!'.format(param, parameters[param]))
      found_error = True
  for param in params_str:
    if not type(parameters[param]) == str:
      logger.error('Input parameter "{0}={1}" should be a string!'.format(param, parameters[param]))
      found_error = True
  for param in params_bool:
    if not type(parameters[param]) == bool:
      logger.error('Input parameter "{0}={1}" should be a boolean!'.format(param, parameters[param]))
      found_error = True
  # Ensure that first two continuation parameter are sorted
  if not os.path.isdir(parameters['running_dir']):
    os.mkdir(parameters['running_dir'])
    logger.info('Created running directory "{0}"'.format(os.path.realpath(parameters['running_dir'])))
  # Check that all input folders and files exist
  exec_loc = parameters['exec_loc']
  tmp = [elt.strip() for elt in exec_loc.split('/')]
  if '~' in tmp:
    home_dir = expanduser('~')
    index = tmp.index('~')
    tmp.remove('~')
    tmp.insert(index, home_dir)
    exec_loc = '/'.join(tmp)
    logger.debug('exec_loc="{0}"'.format(exec_loc))
  if not os.path.isfile(exec_loc):
    logger.error('Input parameter "exec_loc" does not point to an existing Redback executable')
    found_error = True
  if not os.path.isdir(parameters['input_dir']):
    logger.error('Input parameter "input_dir" does not point to an existing directory')
    found_error = True
  if parameters['ref_s_curve']:
    if not os.path.isfile(parameters['ref_s_curve']):
      logger.error('Input parameter "ref_s_curve" does not point to an existing file')
      found_error = True
    parameters['ref_s_curve'] = os.path.realpath(parameters['ref_s_curve'])
  # replace directories with full path as we're going to change working directory
  parameters['input_dir'] = os.path.realpath(parameters['input_dir'])
  parameters['running_dir'] = os.path.realpath(parameters['running_dir'])
  parameters['exec_loc'] = os.path.realpath(exec_loc)
  parameters['result_curve_csv'] = os.path.realpath(parameters['result_curve_csv'])
  # Create running directory if required
  if not os.path.isdir(parameters['running_dir']):
    os.mkdir(parameters['running_dir'])
    logger.info('Created running directory "{0}"'.format(os.path.realpath(parameters['running_dir'])))

  return found_error

def runInitialSimulation1(parameters, logger):
  ''' Run initial simulation (the very first one) with provided value of lambda 
      @param[in] parameters - dictionary of input parameters
      @param[in] logger - python logger instance
  '''
  input_file = os.path.join(parameters['input_dir'], 'extra_param_initial_guess1.i')
  command1 = '{exec_loc} --n-threads={nb_procs} '\
              '-i {input_i} Outputs/csv=true '\
              'Materials/adim_rock/gr={gr}'\
              .format(nb_procs=parameters['nb_threads'], exec_loc=parameters['exec_loc'],
                      input_i=input_file, gr=parameters['lambda_initial_1'])
  try:
    logger.debug(command1)
    stdout = subprocess.check_output(command1.split())
    checkMooseOutput(stdout, logger)
  except:
    logger.error('Execution failed! (First initialisation step)')
    sys.exit(1)

def runInitialSimulation2(parameters, logger):
  ''' Run initial simulation (the second one) with provided value of lambda 
      @param[in] parameters - dictionary of input parameters
      @param[in] logger - python logger instance
  '''
  input_file = os.path.join(parameters['input_dir'], 'extra_param_initial_guess2.i')
  command2 = '{exec_loc} --n-threads={nb_procs} '\
              '-i {input_i} Outputs/csv=true '\
              'Materials/adim_rock/gr={gr}'\
              .format(nb_procs=parameters['nb_threads'], exec_loc=parameters['exec_loc'],
                      input_i=input_file, gr=parameters['lambda_initial_2'])
  try:
    logger.debug(command2)
    stdout = subprocess.check_output(command2.split())
    checkMooseOutput(stdout, logger)
  except:
    logger.error('Execution failed! (Second initialisation step)')
    sys.exit(1)

def parseCsvFile(csv_filename):
  ''' Parse csv file to extract latest value of lambda and max_temp '''
  logger.debug('Parsing csv file "{0}"'.format(csv_filename))
  latest_lambda = None
  latest_max_temp = None
  index_column_lambda = None
  index_column_max_temp = None
  with open(csv_filename, 'rb') as csvfile:
    csvreader = csv.reader(csvfile)
    line_i = 0 # line index
    for row in csvreader:
      if line_i == 0:
        # Headers. Find the column "lambda"
        headers = row
        if 'lambda' in headers:
          index_column_lambda = headers.index('lambda')
        if 'max_temp' in headers:
          index_column_max_temp = headers.index('max_temp')
        line_i += 1
        continue
      # Data line
      if len(row) < len(headers):
          break # finished reading all data
      if index_column_lambda is not None:
        latest_lambda = float(row[index_column_lambda])
      if index_column_max_temp is not None:
        latest_max_temp = float(row[index_column_max_temp])
      line_i += 1
      continue # go to next data line
  return latest_lambda, latest_max_temp

def writeResultsToCsvFile(results, step_index):
  ''' Write results to S-curve csv file for give step 
      @param[in] results - dictionary of results returned by runContinuation()
      @step_index[in] - int, step index
  '''
  with open(parameters['result_curve_csv'], 'a') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow([step_index, '{0:3.16e}'.format(results[step_index][0]), 
                        '{0:3.16e}'.format(results[step_index][1])])

def parseScurveCsv(parameters, logger):
  ''' Parse S-curve csv file
      @param[in] parameters - dictionary of input parameters
      @param[in] logger - python logger instance
      @return lambda_vals - list of continuation values (float)
      @return max_temp_vals - list of norms (float)
  '''
  logger.debug('Parsing csv file "{0}"'.format(parameters['result_curve_csv']))
  lambda_vals = []
  max_temp_vals = []
  with open(parameters['result_curve_csv'], 'rb') as csvfile:
    csvreader = csv.reader(csvfile)
    line_i = 0 # line index
    for row in csvreader:
      if line_i == 0:
        # Headers
        line_i += 1
        continue
      # Data line
      if len(row) < 3:
          break # finished reading all data
      lambda_vals.append(float(row[1]))
      max_temp_vals.append(float(row[2]))
      line_i += 1
      continue # go to next data line
  return lambda_vals, max_temp_vals

def runContinuation(parameters, logger):
  ''' Master function to run pseudo arc-length continuation 
      @param[in] logger - python logger instance
  '''
  logger.info('='*20+' Starting continuation... '+'='*20)
  found_error = checkAndCleanInputParameters(parameters, logger)
  if found_error:
    return
  createRedbackFilesRequired(parameters, logger)
  initial_cwd = os.getcwd()
  os.chdir(parameters['running_dir'])
  # Pseudo arc-length continuation algorithm
  s = 0
  step_index = 0
  lambda_old = 0
  lambda_older = 0
  ds = parameters['ds_initial']
  results = {} # key=step_index, value=[lambda. max_temp]
  
  logger.info('Step {0} (first initial)'.format(step_index))
  runInitialSimulation1(parameters, logger)
  lambda_old = parameters['lambda_initial_1']
  dummy, max_temp = parseCsvFile('extra_param_initial_guess1.csv')
  results[step_index] = [lambda_old, max_temp]
  writeResultsToCsvFile(results, step_index)
  
  step_index += 1
  logger.info('Step {0} (second initial)'.format(step_index))
  runInitialSimulation2(parameters, logger)
  lambda_older = parameters['lambda_initial_1']
  lambda_old = parameters['lambda_initial_2']
  dummy, max_temp = parseCsvFile('extra_param_initial_guess2.csv')
  results[step_index] = [lambda_old, max_temp]
  writeResultsToCsvFile(results, step_index)
  
  finished = False
  #raw_input('About to start the first iterative step.\nPress any key to continue...')
  while not finished:
    step_index += 1
    logger.info('Step {0}, s={1}'.format(step_index, s))
    ds_old = ds 
    ds = parameters['ds_initial']*step_index
    # run simulation
    input_file = os.path.join(parameters['input_dir'], 'extra_param_iteration.i')
    logger.warning('TODO: fix nodes IDs in C++ code to get programmatically all the node IDs')
    lambda_ic = 2*lambda_old - lambda_older 
    previous_exodus_filename = 'extra_param_initial_guess2.e'
    if step_index > 2:
      previous_exodus_filename = 'extra_param_iteration.e'
    command1 = '{exec_loc} --n-threads={nb_procs} '\
                '-i {input_i} Outputs/csv=true Mesh/file={mesh_file} '\
                'Variables/lambda/initial_condition={lambda_IC} '\
                'GlobalParams/ds={ds} GlobalParams/ds_old={ds_old} '\
                'ScalarKernels/continuation_kernel/continuation_parameter_old={lambda_old_value} '\
                'ScalarKernels/continuation_kernel/continuation_parameter_older={lambda_older_value} '\
                'UserObjects/old_temp_UO/mesh={previous_exodus} '\
                'UserObjects/older_temp_UO/mesh={previous_exodus} '\
                .format(nb_procs=parameters['nb_threads'], exec_loc=parameters['exec_loc'],
                        input_i=input_file, mesh_file='extra_param_initial_guess2.e', # TODO: previous exodus file
                        ds=ds, ds_old=ds_old, lambda_old_value=lambda_old, lambda_older_value=lambda_older,
                        lambda_IC=lambda_ic, previous_exodus=previous_exodus_filename)
    try:
      logger.debug(command1)
      stdout = subprocess.check_output(command1.split())
      checkMooseOutput(stdout, logger)
    except:
      logger.error('Execution failed! (Iteration step={0})'.format(step_index))
      sys.exit(1)
    # update lambda_ic
    lambda_older = lambda_old
    lambda_old, max_temp = parseCsvFile('extra_param_iteration.csv')
    results[step_index] = [lambda_old, max_temp]
    writeResultsToCsvFile(results, step_index)
    
    # check if finished
    if (s > parameters['s_max']):
      finished = True
    # increment s
    ds = parameters['ds_initial']
    s += ds
    
    if parameters['plot_s_curve']:
      plotSCurve(parameters, logger)
    #raw_input('Press any key to continue...')
  
  # Finished, clean up
  os.chdir(initial_cwd)
  return results

if __name__ == "__main__":
  # User input
  outpud_dir = '.'
  ds = 1e-2
  parameters = {
    'lambda_initial_1':ds,
    'lambda_initial_2':2*ds,
    'ds_initial':ds,
    's_max':1,
    # Numerical parameters
    'exec_loc':'~/projects/redback/redback-opt',
    'nb_threads':1,
    'input_dir':'benchmark_1_T',
    'running_dir':'running_tmp',
    'result_curve_csv':'S_curve.csv',
    'plot_s_curve':False,
    'ref_s_curve':'benchmark_1_T/ref.csv'
  }
  logger = getLogger('sim', os.path.join(outpud_dir, 'log.txt'), logging.INFO)
  results = runContinuation(parameters, logger)
  print 'Finished'