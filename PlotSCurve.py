''' Script to plot S-curve '''

import os, sys, random, logging, subprocess, shutil, math, csv
from os.path import expanduser
import numpy as np
import pylab as P
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

from Utils import getLogger

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

def plotSCurve(parameters, logger):
  ''' Plot S-Curve
      @param[in] parameters - dictionary of input parameters
      @param[in] logger - python logger instance
  '''
  # parse csv file with values to plot
  lambda_vals, max_temp_vals = parseScurveCsv(parameters, logger)
  # plot
  fig = plt.figure()
  ax_reload = plt.axes([0.02, 0.02, 0.1, 0.075])
  reload_button = Button(ax_reload, 'Reload')
  ax = fig.add_subplot(1,1,1)
  plt.subplots_adjust(left=0.1, bottom=0.15, right=0.97, top=0.95,
                      wspace=None, hspace=None)

  # Plot reference S-curve
  ref_x_values = []
  ref_y_values = []
  if 'ref_s_curve' in parameters and parameters['ref_s_curve']:
    #import pdb;pdb.set_trace()
    with open(parameters['ref_s_curve'], 'rb') as csvfile:
      csvreader = csv.reader(csvfile)
      line_i = 0 # line index
      for row in csvreader:
        if len(row) < 2:
            break # finished reading all data
        ref_x_values.append(float(row[0]))
        ref_y_values.append(float(row[1]))
        line_i += 1
        continue # go to next data line

  plt.xlabel('Continuation parameter', fontsize=20)
  plt.ylabel('Norm of the solution', fontsize=20)

  plt.plot(ref_x_values, ref_y_values,'r-x')
  plt.hold(True)
  x_values = np.array(lambda_vals)
  y_values = np.array(max_temp_vals)
  plt.plot(x_values, y_values,'-x', markerfacecolor='black', markeredgecolor='black',
           markersize=12)

  def reload(event):
    ''' Function to reload curves '''
    lambda_vals, max_temp_vals = parseScurveCsv(parameters, logger)
    x_values = np.array(lambda_vals)
    y_values = np.array(max_temp_vals)
    P.hold(True)
    [xmin, xmax, ymin, ymax] = plt.axis()
    plt.plot(x_values, y_values,'b-x', markerfacecolor='black', markeredgecolor='black',
             markersize=12)
    P.axis([xmin, xmax, ymin, ymax])
    P.hold(False)
    plt.draw()
    print 'Reloaded {0} with {1} data points'.format(parameters['result_curve_csv'], len(x_values))

  reload_button.on_clicked(reload)
  #plt.axes().set_aspect('equal', 'datalim')
  plt.show()

if __name__ == "__main__":
  parameters = {
    'result_curve_csv':'S_curve.csv',
    'ref_s_curve':'benchmark_1_T/ref.csv'
  }
  plotSCurve(parameters, getLogger('plotSCurve', level=logging.INFO))
  print 'Finished'
