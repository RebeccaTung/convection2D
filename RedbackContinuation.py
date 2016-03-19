''' Script to run Redback continuation.
'''

import os, sys, random, logging, subprocess, shutil, math, csv, copy
from os.path import expanduser
import matplotlib.pyplot as plt

from Utils import getLogger, getListFromString, getListOfActiveVariableNames
from CheckMooseOutput import checkMooseOutput, MooseException
from PlotSCurve import plotSCurve
from MooseInputFileRW import MooseInputFileRW
from GenerateInputFiles import writeInitialGuessFile, writeIterationFile, \
  CONT_PARAM_NAMES, SIM_IG1_NAME, SIM_IG2_NAME, SIM_ITER_NAME, \
  updateMultiplyingCoefficientsForInitialGuess

def createRedbackFilesRequired(parameters, handler, logger):
  ''' Create Redback files required to run continuation
      Updates parameters dictionary
      @param[in] parameters - dictionary of input parameters
      @param[in] handler - MooseInputFileRW instance
      @param[in] logger - python logger instance
      @return sim_1 - dictionary representing simulation for first initial guess
      @return sim_2 - dictionary representing simulation for first initial guess
      @return sim_i - dictionary representing simulation for iteration steps
      @return nb_vars - int, number of variables in simulation (not counting continuation parameters)
  '''
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
  if file_content_1.strip() != file_content_2.strip():
    error_msg = 'Reader/writer did not reproduce exactly the input file "{0}". '.format(filename_in)
    error_msg += 'Check manually if this is acceptable and overwrite this error message if necessary.'
    raise Exception, error_msg
  else:
    logger.debug('Our reader/writer successfully read and wrote the input file')
  # Count number of variables
  variable_names = getListOfActiveVariableNames(data_sim, logger)
  nb_vars = len(variable_names)
  # Create 3 simulation files (IG1, IG1, iteration)
  parameters['input_IG1'] = os.path.join(parameters['running_dir'], '{0}.i'.format(SIM_IG1_NAME))
  parameters['input_IG2'] = os.path.join(parameters['running_dir'], '{0}.i'.format(SIM_IG2_NAME))
  parameters['input_iteration'] = os.path.join(parameters['running_dir'], '{0}.i'.format(SIM_ITER_NAME))
  sim_1 = writeInitialGuessFile(1, data_sim, variable_names, parameters, handler, logger)
  sim_2 = writeInitialGuessFile(2, data_sim, variable_names, parameters, handler, logger)
  sim_i = writeIterationFile(data_sim, variable_names, parameters['input_iteration'], handler, logger, parameters)
  # Create output csv file for S-curve results
  if os.path.isfile(parameters['result_curve_csv']):
      os.remove(parameters['result_curve_csv'])
  with open(parameters['result_curve_csv'], 'wb') as csvfile:
    csvwriter = csv.writer(csvfile)
    row = ['step', 'lambda']
    for i, variable_name in enumerate(variable_names):
      row.append('norm_L_inf_u{0} ({1})'.format(i, variable_name))
      row.append('norm_L2_u{0} ({1})'.format(i, variable_name))
    csvwriter.writerow(row)
  
  return sim_1, sim_2, sim_i, nb_vars

def checkAndCleanInputParameters(parameters, logger):
  ''' Check input parameters provided by user and add default values where needed.
      @param[in] parameters - dictionary of input parameters
      @param[in] logger - python logger instance
      @return: found_error - boolean, True if any error was found
  '''
  continuation_variables = ['Gruntfest', 'Lewis'] # Acceptable continuation variables
  plot_norms = ['L_inf', 'L2'] # Acceptable norms to plot
  found_error = False
  # Check all parameters are provided
  if type(parameters) != dict:
    logger.error('"parameters" must be of type dictionary, not {0}!'.format(type(parameters)))
    return True
  # add defaults
  if 'continuation_variable' not in parameters:
    parameters['continuation_variable'] = 'Gruntfest'
  if 'plot_norm' not in parameters:
    parameters['plot_norm'] = 'L_inf'
  if 'plot_solution_index' not in parameters:
    parameters['plot_solution_index'] = 0
  if 'step_change_factor' not in parameters:
    parameters['step_change_factor'] = 0.25
  if 'ref_s_curve' not in parameters:
    parameters['ref_s_curve'] = ''
  if 'plot_s_curve' not in parameters:
    parameters['plot_s_curve'] = False
  if 'non_blocking' not in parameters:
    parameters['non_blocking'] = True
  # check entries
  params_real = ['lambda_initial_1', 'lambda_initial_2', 'ds_initial', 's_max', 'rescaling_factor',
                 'step_change_factor']
  params_int = ['nb_threads', 'plot_solution_index']
  params_str = ['exec_loc', 'input_file', 'running_dir', 'result_curve_csv', 'ref_s_curve', 
                'error_filename', 'continuation_variable', 'plot_norm']
  params_bool = ['plot_s_curve', 'non_blocking']
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
  # Check continuation variable
  if not parameters['continuation_variable'] in continuation_variables:
    logger.error('The continuation variable must be in {0}'.format(continuation_variables))
    found_error = True
  # Check plot options
  if not parameters['plot_norm'] in plot_norms:
    logger.error('The "plot_norm" parameter must be in {0}'.format(plot_norms))
    found_error = True
  if not parameters['plot_solution_index'] > -1:
    logger.error('The "plot_solution_index" parameter must be positive')
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

def getRedbackMaterialNames(sim_data, logger):
  ''' Get list of all RedbackMaterial names in this simulation
      @param[in] sim_data - python structure containing simulation data.
      @param[in] logger - python logger instance
      @return: redback_mat_names - list of material names
  '''
  top_block_names = [elt['name'] for elt in sim_data['children']]
  materials_index = top_block_names.index('Materials')
  materials = sim_data['children'][materials_index]
  redback_mat_names = []
  for material in materials['children']:
    material_name = material['name']
    for elt in material['attributes']:
      if elt['name']=='type' and elt['value']=='RedbackMaterial':
        redback_mat_names.append(material_name)
  return redback_mat_names

def runInitialSimulation1(parameters, sim_data, logger):
  ''' Run initial simulation (the very first one) with provided value of lambda
      @param[in] parameters - dictionary of input parameters
      @param[in] sim_data - python structure containing simulation data.
      @param[in] logger - python logger instance
  '''
  cont_param_name = CONT_PARAM_NAMES[parameters['continuation_variable']]
  input_file = os.path.join(parameters['running_dir'], '{0}.i'.format(SIM_IG1_NAME))
  # Find all materials
  redback_mat_names = getRedbackMaterialNames(sim_data, logger)
  # Form command
  command1 = '{exec_loc} --n-threads={nb_procs} '\
              '-i {input_i} Outputs/csv=true '.\
              format(nb_procs=parameters['nb_threads'], exec_loc=parameters['exec_loc'],
                     input_i=input_file)
  for material_name in redback_mat_names:
    command1 += 'Materials/{mat_name}/{attr_name}={value} '\
      .format(mat_name=material_name, attr_name=cont_param_name,
              value=parameters['lambda_initial_1']*parameters['rescaling_factor'])
  try:
    logger.debug(command1)
    stdout = subprocess.check_output(command1.split())
    checkMooseOutput(stdout, parameters['error_filename'], logger)
  except:
    logger.error('Execution failed! (First initialisation step)')
    sys.exit(1)

def runInitialSimulation2(parameters, sim_data, logger):
  ''' Run initial simulation (the second one) with provided value of lambda
      @param[in] parameters - dictionary of input parameters
      @param[in] sim_data - python structure containing simulation data.
      @param[in] logger - python logger instance
  '''
  cont_param_name = CONT_PARAM_NAMES[parameters['continuation_variable']]
  input_file = os.path.join(parameters['running_dir'], '{0}.i'.format(SIM_IG2_NAME))
  # Find all materials
  redback_mat_names = getRedbackMaterialNames(sim_data, logger)
  # Form command
  command2 = '{exec_loc} --n-threads={nb_procs} '\
              '-i {input_i} Outputs/csv=true '.\
              format(nb_procs=parameters['nb_threads'], exec_loc=parameters['exec_loc'],
                     input_i=input_file)
  for material_name in redback_mat_names:
    command2 += 'Materials/{mat_name}/{attr_name}={value} '\
      .format(mat_name=material_name, attr_name=cont_param_name,
              value=parameters['lambda_initial_2']*parameters['rescaling_factor'])
  try:
    logger.debug(command2)
    stdout = subprocess.check_output(command2.split())
    checkMooseOutput(stdout, parameters['error_filename'], logger)
  except:
    logger.error('Execution failed! (Second initialisation step)')
    sys.exit(1)

def parseCsvFile(csv_filename, nb_vars, logger):
  ''' Parse csv file to extract latest values from file, including
      lambda, L_inf_norm_u, L2_norm_u_diff, nli, nnli
      @param[in] csv_filename - string, name of csv file to parse
      @param[in] nb_vars - int, number of variables active in the simulation
      @param[in] logger - python logger instance
      @return data - dictionary of data with following keys:
        ['lambda','L_inf_norm_u','L2_norm_u','L2_norm_u_diff','nli','nnli']
  '''
  logger.debug('Parsing csv file "{0}"'.format(csv_filename))
  latest_lambda = None
  latest_L_inf_norm = [None]*nb_vars
  latest_L2_norm = [None]*nb_vars
  latest_L2norm_diff = [None]*nb_vars
  latest_nli = None
  latest_nnli = None
  index_column_lambda = None
  index_column_L_inf_norm = [None]*nb_vars
  index_column_L2_norm = [None]*nb_vars
  index_L2norm_diff = [None]*nb_vars
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
        for k in range(nb_vars):
          # solutions L_inf norms
          label = 'L_inf_norm_u{0}'.format(k)
          if label in headers:
            index_column_L_inf_norm[k] = headers.index(label)
          # solutions L2 norms
          label = 'L2_norm_u{0}'.format(k)
          if label in headers:
            index_column_L2_norm[k] = headers.index(label)
          # solutions diff norms
          label = 'L2_norm_u{0}_diff'.format(k)
          if label in headers:
            index_L2norm_diff[k] = headers.index(label)
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
      for k in range(nb_vars):
        # solutions L_inf norms
        if index_column_L_inf_norm[k] is not None:
          latest_L_inf_norm[k] = float(row[index_column_L_inf_norm[k]])
        # solutions L2 norms
        if index_column_L2_norm[k] is not None:
          latest_L2_norm[k] = float(row[index_column_L2_norm[k]])
        # solutions diff norms
        label = 'L2_norm_u{0}_diff'.format(k)
        if index_L2norm_diff[k] is not None:
          latest_L2norm_diff[k] = float(row[index_L2norm_diff[k]])
      if index_nli is not None:
        latest_nli = float(row[index_nli])
      if index_nnli is not None:
        latest_nnli = float(row[index_nnli])
      line_i += 1
      continue # go to next data line
  data = {
    'lambda':latest_lambda,
    'L_inf_norm_u':latest_L_inf_norm,
    'L2_norm_u':latest_L2_norm,
    'L2_norm_u_diff':latest_L2norm_diff,
    'nli':latest_nli,
    'nnli':latest_nnli
  }
  return data

def writeResultsToCsvFile(results, step_index, parameters, variable_names):
  ''' Write results to S-curve csv file for give step
      @param[in] results - dictionary of results returned by parseCsvFile()
      @param[in] step_index - int, step index
      @param[in] parameters - dictionary of input parameters
      @param[in] variable_names - list of strings, simulation variable names
  '''
  with open(parameters['result_curve_csv'], 'a') as csvfile:
    csvwriter = csv.writer(csvfile)
    row = [step_index, '{0:3.16e}'.format(results[step_index]['lambda'])]
    for i, variable_name in enumerate(variable_names):
      row.append('{0:3.16e}'.format(results[step_index]['L_inf_norm_u'][i]))
      row.append('{0:3.16e}'.format(results[step_index]['L2_norm_u'][i]))
    csvwriter.writerow(row)

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

def computeDsForPreviousStep(results, step_index, lambda_old, lambda_older, 
  nb_vars, logger):
  ''' Compute length "ds" for previous step.
      @param[in] results - dictionary of results as returned by parseCsvFile()
      @param[in] step_index - int, step index
      @param[in] lambda_old - float, old value of continuation parameter
      @param[in] lambda_older - float, older value of continuation parameter
      @param[in] nb_vars - int, number of variables in the simulation
      @param[in] logger - python logger instance
      @return ds - float, length of continuation parameter increment
  '''
  var_diff_squared = 0
  for k in range(nb_vars):
    var_diff_squared += results[step_index]['L2_norm_u_diff'][k]**2
  ds = math.sqrt(var_diff_squared + (lambda_old - lambda_older)**2)
  return ds

def runContinuation(parameters, logger):
  ''' Master function to run pseudo arc-length continuation
      @param[in] parameters - dictionary of user input
      @param[in] logger - python logger instance
  '''
  MAX_ATTEMPTS = 8 # Maximum attempts to try and solve any given iteration
  logger.info('='*20+' Starting continuation... '+'='*20)
  found_error = checkAndCleanInputParameters(parameters, logger) # adds missing defaults too
  if found_error:
    return
  handler = MooseInputFileRW()
  sim_1, sim_2, sim_i, nb_vars = createRedbackFilesRequired(parameters, handler, logger)
  variable_names = getListOfActiveVariableNames(sim_1, logger)
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
  runInitialSimulation1(parameters, sim_1, logger)
  lambda_old = parameters['lambda_initial_1']
  results[step_index] = parseCsvFile('{0}.csv'.format(SIM_IG1_NAME), nb_vars, logger)
  results[step_index]['lambda'] = parameters['lambda_initial_1']*parameters['rescaling_factor']
  writeResultsToCsvFile(results, step_index, parameters, variable_names)
  if parameters['plot_s_curve']:
    fig_name = plotSCurve(parameters, logger, figure_name=None)

  step_index += 1
  logger.info('Step {0} (second initial)'.format(step_index))
  runInitialSimulation2(parameters, sim_2, logger)
  lambda_older = parameters['lambda_initial_1']
  lambda_old = parameters['lambda_initial_2']
  results[step_index] = parseCsvFile('{0}.csv'.format(SIM_IG2_NAME), nb_vars, logger)
  results[step_index]['lambda'] = parameters['lambda_initial_2']*parameters['rescaling_factor']
  writeResultsToCsvFile(results, step_index, parameters, variable_names)
  attempt_index = 1 # This second initialisation step succeeded in 1 attempt
  ds = computeDsForPreviousStep(results, step_index, lambda_old, lambda_older, nb_vars, logger)
  if parameters['plot_s_curve']:
    plotSCurve(parameters, logger, figure_name=fig_name)

  input_file = os.path.join(parameters['running_dir'], '{0}.i'.format(SIM_ITER_NAME))
  finished = False
  while not finished:
    #raw_input('About to start iterative step.\nPress enter to continue...')
    #ds_old = computeDsForPreviousStep(results, step_index, lambda_old, lambda_older, nb_vars, logger)
    step_index += 1
    ds_old = ds
    ds = getInitialStepLength(ds_old, parameters['ds_initial'], attempt_index, logger)
    previous_exodus_filename = '{0}.e'.format(SIM_IG2_NAME)
    if step_index > 2:
      previous_exodus_filename = '{0}.e'.format(SIM_ITER_NAME)
    # save that previous exodus file in case the coming step requires multiple attempts
    root_name, ext = os.path.splitext(previous_exodus_filename)
    backup_exodus_filename = root_name + "_BACKUP" + ext
    shutil.copyfile(previous_exodus_filename, backup_exodus_filename)
    step_succeeded = False
    attempt_index = 0
    while not step_succeeded and attempt_index < MAX_ATTEMPTS:
      logger.info('Step {0} (attempt {1}), s={2}, ds={3}'\
                  .format(step_index, attempt_index, s, ds))
      # use backup file if not in the first attempt
      if attempt_index > 0:
        shutil.copyfile(backup_exodus_filename, previous_exodus_filename)
      # Calculate multiplying coefficients of old and older solutions for coming initial guess
      coeff_mult = ds/ds_old
      coeff_guess_old = 1 + coeff_mult
      coeff_guess_older = -coeff_mult
      lambda_ic = coeff_guess_old*lambda_old + coeff_guess_older*lambda_older
      sim_i = updateMultiplyingCoefficientsForInitialGuess\
        (sim_i, coeff_guess_old, coeff_guess_older, nb_vars)
      handler.write(sim_i, parameters['input_iteration'])
      # run simulation
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
        ds = ds*parameters['step_change_factor']
    if attempt_index >= MAX_ATTEMPTS:
      logger.error('Execution failed after {0} attempts!'\
                   .format(attempt_index))
      sys.exit(1)
    # update lambda_ic
    lambda_older = lambda_old
    results[step_index] = parseCsvFile('{0}.csv'.format(SIM_ITER_NAME), nb_vars, logger)
    lambda_old = results[step_index]['lambda']
    results[step_index]['lambda'] *= parameters['rescaling_factor']
    writeResultsToCsvFile(results, step_index, parameters, variable_names)

    # increment s
    s += math.fabs(ds)
    # check if finished
    if (s > parameters['s_max']):
      finished = True

    if parameters['plot_s_curve']:
      plotSCurve(parameters, logger, figure_name=fig_name)
  
  # Finished, clean up
  os.chdir(initial_cwd)
  return results

if __name__ == "__main__":
  # User input
  outpud_dir = '.'
  parameters = {
    'continuation_variable':'Gruntfest', # in ['Gruntfest', 'Lewis']
    'lambda_initial_1':1e-8, #ds,
    'lambda_initial_2':2e-8, #2*ds,
    'ds_initial':5e-2,
    's_max':0.5,
    'rescaling_factor':4.5399929762e-5, # to multiply continuation parameter
    # Rescaling factor
    'rescaling_factor':1, # to multiply continuation parameter
    # Numerical parameters
    'exec_loc':'~/projects/redback/redback-opt',
    'nb_threads':1,
    'input_file':'input_files/benchmark_1_T/bench1_a.i',#'benchmark_4_TH/bench_TH.i',
    'running_dir':'input_files/benchmark_1_T/running_tmp',
    'result_curve_csv':'input_files/benchmark_1_T/S_curve.csv',
    'error_filename':'error_output.txt',
    # step refinement
    'step_change_factor':0.25, # multiplying factor of step size when step fails
    # plot
    'plot_s_curve':True,
    'non_blocking':True,
    'ref_s_curve':'input_files/benchmark_1_T/ref.csv',
    'plot_norm':'L_inf', # in ['L2', 'L_inf']
    'plot_solution_index':0, # index of solution to plot
  }
  logger = getLogger('sim', os.path.join(outpud_dir, 'log.txt'), logging.INFO)
  results = runContinuation(parameters, logger)
  print 'Finished'
