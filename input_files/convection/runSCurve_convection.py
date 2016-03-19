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
    'lambda_initial_1':4,
    'lambda_initial_2':4.1,
    'ds_initial':1e-1,
    's_max':4,
    # Rescaling factor
    'rescaling_factor':1e-8, # to multiply continuation parameter
    # Numerical parameters
    'exec_loc':'~/projects/redback/redback-opt',
    'nb_threads':1,
    'input_file':'convection2D.i',
    'running_dir':'running_tmp',
    'result_curve_csv':'S_curve.csv',
    'error_filename':'error_output.txt',
    'plot_s_curve':True,
    'non_blocking':True,
    'plot_post_processor':'Nusselt_number', # post-processor to plot. If empty, then norm is used.
    'plot_norm':'L2', # in ['L2', 'L_inf']
    'plot_solution_index':1, # index of solution to plot
    'ref_s_curve':'',
    'postprocessors_exported':['Nusselt_number'],
  }
  logger = getLogger('sim', os.path.join(outpud_dir, 'log.txt'), logging.INFO)
  results = runContinuation(parameters, logger)
  raw_input('press enter to finish')
  print 'Finished'
