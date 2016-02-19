''' Script to run Redback continuation.
  TODO: everything...
'''

import os, sys, random, logging, subprocess, shutil, math, csv, copy
from os.path import expanduser

from Utils import getLogger
from CheckMooseOutput import checkMooseOutput, MooseException
from PlotSCurve import plotSCurve
from MooseInputFileRW import MooseInputFileRW
from GenerateInputFiles import writeInitialGuessFile, writeIterationFile, \
  SIM_IG1_NAME, SIM_IG2_NAME, SIM_ITER_NAME, updateMultiplyingCoefficientsForInitialGuess

def createRedbackFilesRequired(parameters, handler, logger):
  ''' Create Redback files required to run continuation
      Updates parameters dictionary
      @param[in] parameters - dictionary of input parameters
      @param[in] handler - MooseInputFileRW instance
      @param[in] logger - python logger instance
      @return sim_i - dictionary representing simulation for iteration steps
  '''
  # Create output csv file for S-curve results
  if os.path.isfile(parameters['result_curve_csv']):
      os.remove(parameters['result_curve_csv'])
  with open(parameters['result_curve_csv'], 'wb') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['step', 'lambda', 'norm'])
  # Check that our reader/writer can reproduce that input file faithfully
  filename_in = parameters['input_file']
  root, ext = os.path.splitext(filename_in)
  assert ext == '.i'
  filename_out = os.path.join(parameters['running_dir'], 'input_recreated.i')
  if handler == None:
    handler = MooseInputFileRW()
  data_sim = handler.read(filename_in)
  handler.write(data_sim, filename_out)
  with open(filename_in, 'r') as f_1:
    file_content_1 = f_1.read()
  with open(filename_out, 'r') as f_2:
    file_content_2 = f_2.read()
  if file_content_1 != file_content_2:
    error_msg = 'Reader/writer did not reproduce exactly the input file "{0}". '.format(filename_in)
    error_msg += 'Check manually if this is acceptable and overwrite this error message if necessary.'
    raise Exception, error_msg
  else:
    logger.debug('Our reader/writer successfully read and wrote the input file')
  # Create 3 simulation files (IG1, IG1, iteration)
  parameters['input_IG1'] = os.path.join(parameters['running_dir'], '{0}.i'.format(SIM_IG1_NAME))
  parameters['input_IG2'] = os.path.join(parameters['running_dir'], '{0}.i'.format(SIM_IG2_NAME))
  parameters['input_iteration'] = os.path.join(parameters['running_dir'], '{0}.i'.format(SIM_ITER_NAME))
  # First initial guess
  if 1:
    sim_1 = writeInitialGuessFile(1, data_sim, parameters['input_IG1'], handler, logger)
    sim_2 = writeInitialGuessFile(2, data_sim, parameters['input_IG2'], handler, logger)
    sim_i = writeIterationFile(data_sim, parameters['input_iteration'], handler, logger)
  else:
    sim_i = None
  return sim_i

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
  params_str = ['exec_loc', 'input_file', 'running_dir', 'result_curve_csv', 'ref_s_curve', 
                'error_filename']
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
  if not os.path.isfile(parameters['input_file']):
    logger.error('Input parameter "input_file" does not point to an existing file')
    found_error = True
  if parameters['ref_s_curve']:
    if not os.path.isfile(parameters['ref_s_curve']):
      logger.error('Input parameter "ref_s_curve" does not point to an existing file')
      found_error = True
    parameters['ref_s_curve'] = os.path.realpath(parameters['ref_s_curve'])
  # replace directories with full path as we're going to change working directory
  parameters['input_file'] = os.path.realpath(parameters['input_file'])
  parameters['running_dir'] = os.path.realpath(parameters['running_dir'])
  parameters['exec_loc'] = os.path.realpath(exec_loc)
  parameters['result_curve_csv'] = os.path.realpath(parameters['result_curve_csv'])
  parameters['error_filename'] = os.path.realpath(parameters['error_filename'])
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
  input_file = os.path.join(parameters['running_dir'], '{0}.i'.format(SIM_IG1_NAME))
  command1 = '{exec_loc} --n-threads={nb_procs} '\
              '-i {input_i} Outputs/csv=true '\
              'Materials/adim_rock/gr={gr}'\
              .format(nb_procs=parameters['nb_threads'], exec_loc=parameters['exec_loc'],
                      input_i=input_file, gr=parameters['lambda_initial_1'])
  try:
    logger.debug(command1)
    stdout = subprocess.check_output(command1.split())
    checkMooseOutput(stdout, parameters['error_filename'], logger)
  except:
    logger.error('Execution failed! (First initialisation step)')
    sys.exit(1)

def runInitialSimulation2(parameters, logger):
  ''' Run initial simulation (the second one) with provided value of lambda
      @param[in] parameters - dictionary of input parameters
      @param[in] logger - python logger instance
  '''
  input_file = os.path.join(parameters['running_dir'], '{0}.i'.format(SIM_IG2_NAME))
  command2 = '{exec_loc} --n-threads={nb_procs} '\
              '-i {input_i} Outputs/csv=true '\
              'Materials/adim_rock/gr={gr}'\
              .format(nb_procs=parameters['nb_threads'], exec_loc=parameters['exec_loc'],
                      input_i=input_file, gr=parameters['lambda_initial_2'])
  try:
    logger.debug(command2)
    stdout = subprocess.check_output(command2.split())
    checkMooseOutput(stdout, parameters['error_filename'], logger)
  except:
    logger.error('Execution failed! (Second initialisation step)')
    sys.exit(1)

def parseCsvFile(csv_filename):
  ''' Parse csv file to extract latest value of lambda, max_temp
      and L2_norm_u_diff
      @param[in] csv_filename - string, name of csv file to parse
      @return data - dictionary of data with following keys:
        ['lambda','solution_norm','L2_norm_u_diff','nli','nnli']
  '''
  logger.debug('Parsing csv file "{0}"'.format(csv_filename))
  latest_lambda = None
  latest_sol_norm = None
  latest_L2norm = None
  latest_nli = None
  latest_nnli = None
  index_column_lambda = None
  index_column_sol_norm = None
  index_L2norm = None
  index_nli = None
  index_nnli = None
  with open(csv_filename, 'rb') as csvfile:
    csvreader = csv.reader(csvfile)
    line_i = 0 # line index
    for row in csvreader:
      if line_i == 0:
        # Headers. Find the column "lambda"
        headers = row
        if 'lambda' in headers:
          index_column_lambda = headers.index('lambda')
        if 'solution_norm' in headers:
          index_column_sol_norm = headers.index('solution_norm')
        if 'L2_norm_u_diff' in headers:
          index_L2norm = headers.index('L2_norm_u_diff')
        if 'nli' in headers:
          index_nli = headers.index('nli')
        if 'nnli' in headers:
          index_nnli = headers.index('nnli')
        line_i += 1
        continue
      # Data line
      if len(row) < len(headers):
          break # finished reading all data
      if index_column_lambda is not None:
        latest_lambda = float(row[index_column_lambda])
      if index_column_sol_norm is not None:
        latest_sol_norm = float(row[index_column_sol_norm])
      if index_L2norm is not None:
        latest_L2norm = float(row[index_L2norm])
      if index_nli is not None:
        latest_nli = float(row[index_nli])
      if index_nnli is not None:
        latest_nnli = float(row[index_nnli])
      line_i += 1
      continue # go to next data line
  data = {
    'lambda':latest_lambda,
    'solution_norm':latest_sol_norm,
    'L2_norm_u_diff':latest_L2norm,
    'nli':latest_nli,
    'nnli':latest_nnli
  }
  return data

def writeResultsToCsvFile(results, step_index):
  ''' Write results to S-curve csv file for give step
      @param[in] results - dictionary of results returned by parseCsvFile()
      @step_index[in] - int, step index
  '''
  with open(parameters['result_curve_csv'], 'a') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow([step_index, '{0:3.16e}'.format(results[step_index]['lambda']),
                        '{0:3.16e}'.format(results[step_index]['solution_norm'])])

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

def getInitialStepLength(ds_old, ds0, previous_nb_attempts, logger):
  ''' Calculates length of initial continuation step new_ds for coming iteration, not exceeding ds0.
      @param[in] ds_old - float, length of previous successful step for previous iteration
      @param[in] ds0 - float, initial length of continuation parameter provided by user
      @param[in] previous_nb_attempts - int, number of attempts it took to succeed in previous iteration
      @param[in] logger - python logger instance
      @return new_ds - float, length of coming iteration
  '''
  logger.debug('getInitialStepLength(ds_old={0}, ds0={1}, previous_nb_attempts={2})'\
    .format(ds_old, ds0, previous_nb_attempts))
  return ds0 # always start with value given by user
  mult_coeff = 1 # mutliplying coefficient to apply to ds_old
  if previous_nb_attempts in [0]:
    # last iteration was too easy, increase step a bit
    mult_coeff = 1.2
  elif previous_nb_attempts in range(1, 5):
    # last iteration required just a few attempt, start from the same step size
    mult_coeff = 1
  else:
    # last iteration took a lot of attempts, reduce step size directly
    mult_coeff = 0.8
  # ensure step size doesn't exceed user provided size
  new_ds = min(mult_coeff*ds_old, ds0)
  logger.debug('  Changing ds_init by {0}x --> new_ds={1}'.format(mult_coeff, new_ds))
  return new_ds

def runContinuation(parameters, logger):
  ''' Master function to run pseudo arc-length continuation
      @param[in] logger - python logger instance
  '''
  MAX_ATTEMPTS = 8 # Maximum attempts to try and solve any given iteration
  logger.info('='*20+' Starting continuation... '+'='*20)
  found_error = checkAndCleanInputParameters(parameters, logger)
  if found_error:
    return
  handler = MooseInputFileRW()
  sim_i = createRedbackFilesRequired(parameters, handler, logger)
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
  results[step_index] = parseCsvFile('{0}.csv'.format(SIM_IG1_NAME))
  results[step_index]['lambda'] = parameters['lambda_initial_1']
  writeResultsToCsvFile(results, step_index)

  step_index += 1
  logger.info('Step {0} (second initial)'.format(step_index))
  runInitialSimulation2(parameters, logger)
  lambda_older = parameters['lambda_initial_1']
  lambda_old = parameters['lambda_initial_2']
  results[step_index] = parseCsvFile('{0}.csv'.format(SIM_IG2_NAME))
  results[step_index]['lambda'] = parameters['lambda_initial_2']
  writeResultsToCsvFile(results, step_index)
  attempt_index = 1 # This second initialisation step succeeded in 1 attempt
  ds = math.sqrt(results[step_index]['L2_norm_u_diff']**2 + (lambda_old - lambda_older)**2)
  
  input_file = os.path.join(parameters['running_dir'], '{0}.i'.format(SIM_ITER_NAME))
  finished = False
  while not finished:
    #raw_input('About to start the first iterative step.\nPress enter to continue...')
    step_index += 1
    ds_old = ds
    ds = getInitialStepLength(ds_old, parameters['ds_initial'], attempt_index, logger)
    step_succeeded = False
    attempt_index = 0
    while not step_succeeded and attempt_index < MAX_ATTEMPTS:
      logger.info('Step {0} (attempt {1}), s={2}, ds={3}'\
                  .format(step_index, attempt_index, s, ds))
      # Calculate multiplying coefficients of old and older solutions for coming initial guess
      coeff_mult = ds/ds_old
      coeff_guess_old = 1 + coeff_mult
      coeff_guess_older = -coeff_mult
      lambda_ic = coeff_guess_old*lambda_old + coeff_guess_older*lambda_older
      sim_i = updateMultiplyingCoefficientsForInitialGuess\
        (sim_i, coeff_guess_old, coeff_guess_older)
      handler.write(sim_i, parameters['input_iteration'])
      # run simulation
      previous_exodus_filename = '{0}.e'.format(SIM_IG2_NAME)
      if step_index > 2:
        previous_exodus_filename = '{0}.e'.format(SIM_ITER_NAME)
      command1 = '{exec_loc} --n-threads={nb_procs} '\
                  '-i {input_i} Outputs/csv=true Mesh/file={previous_exodus} '\
                  'Variables/lambda/initial_condition={lambda_IC} '\
                  'GlobalParams/ds={ds} GlobalParams/ds_old={ds_old} '\
                  'ScalarKernels/continuation_kernel/continuation_parameter_old={lambda_old_value} '\
                  'ScalarKernels/continuation_kernel/continuation_parameter_older={lambda_older_value} '\
                  'UserObjects/old_temp_UO/mesh={previous_exodus} '\
                  'UserObjects/older_temp_UO/mesh={previous_exodus} '\
                  .format(nb_procs=parameters['nb_threads'], exec_loc=parameters['exec_loc'],
                          input_i=input_file, ds=ds, ds_old=ds_old,
                          lambda_old_value=lambda_old, lambda_older_value=lambda_older,
                          lambda_IC=lambda_ic, previous_exodus=previous_exodus_filename)
      try:
        logger.debug(command1)
        stdout = subprocess.check_output(command1.split())
        checkMooseOutput(stdout, parameters['error_filename'], logger)
        # if it passes that point with no Exception raised, then attempt succeeded
        step_succeeded = True
        continue
      except MooseException, e:
        logger.debug('Attempt failed Iteration step={0} with ds={1}'\
                    .format(step_index, ds))
        attempt_index += 1
        ds = ds/4.
    if attempt_index >= MAX_ATTEMPTS:
      logger.error('Execution failed after {0} attempts! (Iteration step={0} with ds={1})'\
                   .format(attempt_index, step_index, ds))
      sys.exit(1)
    # update lambda_ic
    lambda_older = lambda_old
    results[step_index] = parseCsvFile('{0}.csv'.format(SIM_ITER_NAME))
    lambda_old = results[step_index]['lambda']
    writeResultsToCsvFile(results, step_index)

    if 0:
      # Compute ds externally for comparison
      logger.debug('  ds (set) = {0}, sol_l2norm = {1}, ds_old={2}'.format(ds, results[step_index]['L2_norm_u_diff'], ds_old))

    # increment s
    s += ds
    # check if finished
    if (s > parameters['s_max']):
      finished = True

    if parameters['plot_s_curve']:
      plotSCurve(parameters, logger)
  
  # Finished, clean up
  os.chdir(initial_cwd)
  return results

if __name__ == "__main__":
  # User input
  outpud_dir = '.'
  ds = 1e-1
  parameters = {
    'lambda_initial_1':1e-8, #ds,
    'lambda_initial_2':2e-8, #2*ds,
    'ds_initial':ds,
    's_max':15,
    # Numerical parameters
    'exec_loc':'~/projects/redback/redback-opt',
    'nb_threads':1,
    'input_file':'benchmark_1_T/bench1_a.i',
    'running_dir':'running_tmp',
    'result_curve_csv':'S_curve.csv',
    'error_filename':'error_output.txt',
    'plot_s_curve':False,
    'ref_s_curve':'benchmark_1_T/ref.csv'
  }
  logger = getLogger('sim', os.path.join(outpud_dir, 'log.txt'), logging.INFO)
  results = runContinuation(parameters, logger)
  print 'Finished'
