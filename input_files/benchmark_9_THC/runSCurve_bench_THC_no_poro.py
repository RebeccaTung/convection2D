''' Script to run simulation producing S-curve for 1D THC example '''

import os, sys, logging

# Add directory containing RedbackContinuation.py to the path
sys.path.append(os.path.join('..', '..'))

from Utils import getLogger
from RedbackContinuation import runContinuation

if __name__ == "__main__":
  # User input
  outpud_dir = '.'
  ds = 5e-2
  parameters = {
    'lambda_initial_1':1e-8, #ds,
    'lambda_initial_2':2e-8, #2*ds,
    'ds_initial':ds,
    's_max':2.66,
    # Rescaling factor
    'rescaling_factor':4.5399929762e-5, # to multiply continuation parameter
    # Numerical parameters
    'exec_loc':'~/projects/redback/redback-opt',
    'nb_threads':1,
    'input_file':'bench_THC_no_poro.i',
    'running_dir':'running_tmp',
    'result_curve_csv':'S_curve.csv',
    'error_filename':'error_output.txt',
    'plot_s_curve':False,
    'ref_s_curve':'',
  }
  logger = getLogger('sim', os.path.join(outpud_dir, 'log.txt'), logging.INFO)
  results = runContinuation(parameters, logger)
  print 'Finished'
