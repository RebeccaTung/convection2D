''' Script to run simulation producing S-curve for 1D THC example '''

import os, sys, logging

# Add directory containing RedbackContinuation.py to the path
sys.path.append(os.path.join('..', '..'))

from Utils import getLogger
from RedbackContinuation import runContinuation

if __name__ == "__main__":
  # User input
  outpud_dir = '.'
  parameters = {
    'continuation_variable':'Lewis', # in ['Gruntfest', 'Lewis']
    'lambda_initial_1':1e-7,
    'lambda_initial_2':9e-8,
    'ds_initial':-1e-9,
    's_max':1e-7,
    # Rescaling factor
    'rescaling_factor':1, # to multiply continuation parameter
    # Numerical parameters
    'exec_loc':'~/projects/redback/redback-opt',
    'nb_threads':1,
    'input_file':'convection2D.i',
    'running_dir':'running_tmp',
    'result_curve_csv':'S_curve.csv',
    'error_filename':'error_output.txt',
    'plot_s_curve':False,
    'ref_s_curve':'',
  }
  logger = getLogger('sim', os.path.join(outpud_dir, 'log.txt'), logging.INFO)
  results = runContinuation(parameters, logger)
  print 'Finished'
