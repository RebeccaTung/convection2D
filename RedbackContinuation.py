''' Script to run Redback continuation. 
  TODO: everything...
'''

import os, sys, random, logging, subprocess, shutil, math, csv

def getLogger(name, log_file='log.txt', level=logging.INFO):
    ''' Creates logger with given name and level.
        @param[in]   name     - string, logger name
        @param[in]   log_file - string, file name for log file
        @param[in]   level    - int, logging level
        @return[out] logger   - logger instance built
    '''
    logger_name = '{0}_{1}'.format(name, random.random())
    logger = logging.getLogger(logger_name)
    hdlr1 = logging.FileHandler(log_file)
    hdlr2 = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')
    hdlr1.setFormatter(formatter)
    hdlr2.setFormatter(formatter)
    logger.addHandler(hdlr1)
    logger.addHandler(hdlr2)
    logger.setLevel(level)
    return logger

def createRedbackFilesRequired(logger):
  ''' Create Redback files required to run continuation 
      @param[in] logger - python logger instance
  '''
  logger.error('TODO: createRedbackFilesRequired not implemented')

def checkInputParameters(parameters, logger):
  ''' Check input parameters provided by user '''
  logger.error('TODO: checkInputParameters not implemented')

def runInitialSimulation1():
  ''' Run initial simulation (the very first one) with provided value of lambda 
      @param[in] logger - python logger instance
  '''
  logger.error('TODO: input_file hardcoded in runInitialSimulation1')
  input_file = os.path.join(parameters['input_dir'], 'extra_param_initial_guess.i')
  command1 = '{exec_loc} --n-threads={nb_procs} '\
              '-i {input_i} Outputs/csv=true'\
              .format(nb_procs=parameters['nb_threads'], exec_loc=parameters['exec_loc'],
                      input_i=input_file)
  try:
    logger.info(command1)
    retcode = subprocess.call(command1, shell=True)
    if retcode < 0:
        error_msg = 'Child process was terminated by signal {0} (First initialisation step)'.format(retcode)
        logger.error(error_msg)
  except:
        logger.error('Execution failed! (First initialisation step)')

def runInitialSimulation2():
  ''' Run initial simulation (the second one) with provided value of lambda 
      @param[in] logger - python logger instance
  '''
  logger.error('TODO: input_file hardcoded in runInitialSimulation2')
  input_file = os.path.join(parameters['input_dir'], 'extra_param_initial_guess2.i')
  command2 = '{exec_loc} --n-threads={nb_procs} '\
              '-i {input_i} Outputs/csv=true'\
              .format(nb_procs=parameters['nb_threads'], exec_loc=parameters['exec_loc'],
                      input_i=input_file)
  try:
    logger.info(command2)
    retcode = subprocess.call(command2, shell=True)
    if retcode < 0:
        error_msg = 'Child process was terminated by signal {0} (Second initialisation step)'.format(retcode)
        logger.error(error_msg)
  except:
        logger.error('Execution failed! (Second initialisation step)')

def readLambdaFromLastCsv(csv_filename):
  ''' Parse csv file to extract latest value of lambda and max_temp '''
  latest_lambda = None
  with open(csv_filename, 'rb') as csvfile:
    csvreader = csv.reader(csvfile)
    line_i = 0 # line index
    for row in csvreader:
        if line_i == 0:
          # Headers. Find the column "lambda"
          headers = row
          index_column_lambda = headers.index('lambda')
          index_column_max_temp = headers.index('max_temp')
          continue
        # Data line
        if len(row) < len(headers):
            break # finished reading all data
        latest_lambda = row[index_column_lambda]
        latest_max_temp = row[index_column_max_temp]
        line_i += 1
        continue # go to next data line
  logger.info('Finished parsing csv file')
  return latest_lambda, latest_max_temp

def runContinuation(parameters, logger):
  ''' Master function to run pseudo arc-length continuation 
      @param[in] logger - python logger instance
  '''
  checkInputParameters(parameters, logger)
  createRedbackFilesRequired(logger)
  initial_cwd = os.getcwd()
  os.chdir(parameters['running_dir'])
  # Pseudo arc-length continuation algorithm
  s = 0
  step_index = 0
  lambda_old = 0
  lambda_older = 0
  ds = parameters['ds_initial']
  results = {} # key=step_index, value=[lambda. max_temp]
  
  runInitialSimulation1()
  lambda_old = parameters['lambda_initial_1']
  
  step_index = 1
  runInitialSimulation2()
  lambda_older = parameters['lambda_initial_1']
  lambda_old = parameters['lambda_initial_2']
  
  finished = False
  while not finished:
    step_index += 1
    ds_old = ds 
    logger.error('TODO: what is the real value of ds for the first time???')
    ds = parameters['ds_initial']
    # run simulation
    logger.error('TODO: input_file hardcoded in runContinuation')
    input_file = os.path.join(parameters['input_dir'], 'extra_param_iteration.i')
    logger.error('TODO: fix nodes IDs in C++ code to get programmatically all the node IDs')
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
      logger.info(command1)
      retcode = subprocess.call(command1, shell=True)
      if retcode < 0:
          error_msg = 'Child process was terminated by signal {0} (First initialisation step)'.format(retcode)
          logger.error(error_msg)
    except:
          logger.error('Execution failed! (First initialisation step)')
    
    # update lambda_ic
    lambda_older = lambda_old
    logger.error('TODO: get max temp for first 3 simulations')
    if step_index > 2:
      lambda_old, max_temp = readLambdaFromLastCsv('extra_param_iteration.csv')
      results[step_index] = [lambda_old, max_temp]
    else:
      assert step_index == 2
      lambda_old = parameters['lambda_initial_2']

    # check if finished
    if (s > parameters['s_max']):
      finished = True
    # increment s
    ds = parameters['ds_initial']
    s += ds
  # Finished, clean up
  os.chdir(initial_cwd)
  return results

if __name__ == "__main__":
  # User input
  outpud_dir = '.'
  parameters = {
    'lambda_initial_1':1e-6,
    'lambda_initial_2':1.7e-6,
    'ds_initial':0.001,
    's_max':1,
    # Numerical parameters
    'exec_loc':'~/projects/redback/redback-opt',
    'nb_threads':1,
    'input_dir':'~/projects/redback/tests/benchmark_extra_param',
    'running_dir':'running_tmp'
  }
  
  
  logger = getLogger('sim', os.path.join(outpud_dir, 'log.txt'), logging.INFO)
  results = runContinuation(parameters, logger)
  import pdb;pdb.set_trace()
  print 'Finished'