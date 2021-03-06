''' Script to plot S-curve '''

import os, sys, logging

# Add directory containing PlotSCurve.py to the path
sys.path.append(os.path.join('..', '..'))

from Utils import getLogger
from PlotSCurve import plotSCurve


if __name__ == "__main__":
  parameters = {
    'result_curve_csv':'S_curve.csv',
    'ref_s_curve':'',
    'plot_post_processor':'Nusselt_number', # post-processor to plot. If empty, then norm is used.
    'plot_norm':'L2', # in ['L2', 'L_inf'], only if no post-processor to plot is provided
    'plot_solution_index':1, # index of solution to plot
  }
  plotSCurve(parameters, getLogger('plotSCurve', level=logging.INFO))
  print 'Finished'
