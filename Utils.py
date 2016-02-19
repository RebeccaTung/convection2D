''' Utilities '''

import os, sys, logging, random

class NullHandler(logging.Handler):
  def emit(self, record):
    pass

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
  #h = NullHandler()
  #logger.addHandler(h)
  logger.setLevel(level)
  return logger

def getListFromString(text):
  ''' Get list of strings from MOOSE text representing a list of items
      @param[in] text - string, text to parse
      @return list - python list of strings
  '''
  tmp = text.strip()
  if tmp.startswith("'") and tmp.endswith("'"):
    tmp = tmp[1:-1]
  list = tmp.split(' ')
  return list

def getListOfActiveVariableNames(sim_data, logger):
  ''' Returns list of active variable names in simulation.
      @param[in] sim_data - dict, simulation dictionary as returned 
        by RedbackContinuation.createRedbackFilesRequired()
      @param[in] logger - python logger instance
      @return nb_vars - int, number of active variables
  '''
  top_block_names = [elt['name'] for elt in sim_data['children']]
  variables_index = top_block_names.index('Variables')
  variables = sim_data['children'][variables_index]
  attribute_names = [attr['name'] for attr in variables['attributes']]
  active_index = None
  if 'active' in attribute_names:
    active_index = attribute_names.index('active')
  all_variables_names = [elt['name'] for elt in variables['children']]
  # find active variables
  if active_index is None:
    active_index = 0
    active_variables_names = all_variables_names
  else:
    active_variables_names = getListFromString(variables['attributes'][active_index]['value'])
  return active_variables_names

