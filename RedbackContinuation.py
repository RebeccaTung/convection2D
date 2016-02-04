''' Script to run Redback continuation. 
  TODO: everything...
'''

import os, sys, random, logging, subprocess, shutil, math

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
  ds = parameters['ds_initial']
  runInitialSimulation1()
  runInitialSimulation2()
  finished = False
  while not finished:
    ds_old = ds 
    logger.error('TODO: what is the real value of ds for the first time???')
    ds = parameters['ds_initial']
    # run simulation
    logger.error('TODO: input_file hardcoded in runContinuation')
    input_file = os.path.join(parameters['input_dir'], 'extra_param_iteration.i')
    logger.error('TODO: need mechanism INSIDE redback input file to calculate this properly: 2*lambda_old - older_lambda')
    logger.error('TODO: need mechanism INSIDE redback input file to calculate this properly: 2*old_temp - older_temp')
    logger.error('TODO: fix nodes IDs in C++ code to get programmatically all the node IDs')
    logger.error('TODO: get proper values of lambda_old_value and lambda_older_value')
    #import pdb;pdb.set_trace()
    command1 = '{exec_loc} --n-threads={nb_procs} '\
                '-i {input_i} Outputs/csv=true Mesh/file={mesh_file} '\
                'Variables/temp/initial_condition=0.03276648 '\
                'Variables/lambda/initial_condition=2.4e-6 '\
                'GlobalParams/ds={ds} GlobalParams/ds_old={ds_old} '\
                'ScalarKernels/continuation_kernel/continuation_parameter_old={lambda_old_value} '\
                'ScalarKernels/continuation_kernel/continuation_parameter_older={lambda_older_value} '\
                .format(nb_procs=parameters['nb_threads'], exec_loc=parameters['exec_loc'],
                        input_i=input_file, mesh_file='extra_param_initial_guess2.e', # TODO: previous exodus file
                        ds=ds, ds_old=ds_old, lambda_old_value=1.7e-6, lambda_older_value=1e-6)
    try:
      logger.info(command1)
      retcode = subprocess.call(command1, shell=True)
      if retcode < 0:
          error_msg = 'Child process was terminated by signal {0} (First initialisation step)'.format(retcode)
          logger.error(error_msg)
    except:
          logger.error('Execution failed! (First initialisation step)')
    
    logger.error('Breaking loop')
    return
    # check if finished
    if (s > parameters['s_max']):
      finished = True
    # increment s
    ds = parameters['ds_initial']
    s += ds
  # Finished, clean up
  os.chdir(initial_cwd)
  return

if __name__ == "__main__":
  # User input
  outpud_dir = '.'
  parameters = {
    'lambda_initial_1':1e-6,
    'lambda_initial_2':1.7e-6,
    'ds_initial':0.1,
    's_max':1,
    # Numerical parameters
    'exec_loc':'~/projects/redback/redback-opt',
    'nb_threads':1,
    'input_dir':'~/projects/redback/tests/benchmark_extra_param',
    'running_dir':'running_tmp'
  }
  
  
  logger = getLogger('sim', os.path.join(outpud_dir, 'log.txt'), logging.INFO)
  runContinuation(parameters, logger)
  print 'Finished'