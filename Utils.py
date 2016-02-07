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

