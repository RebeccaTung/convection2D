''' Functions to generate input files to run continuation
'''

import os, sys, random, logging, subprocess, shutil, math, csv, copy
from os.path import expanduser

from Utils import getLogger
from MooseInputFileRW import MooseInputFileRW

SIM_IG1_NAME = 'extra_param_initial_guess1'
SIM_IG2_NAME = 'extra_param_initial_guess2'
SIM_ITER_NAME = 'extra_param_iteration'

CONT_PARAM_NAME = 'gr' # hardcoded for now, should disappear (TODO)

def writeInitialGuessFile(index_ig, sim_data, out_filename, handler, logger):
  ''' Write first or second initial guess file
      @param[in] index_ig - int, 1 or 2 (to distinguish first and second initial guess simulation)
      @param[in] sim_data - python structure containing simulation data.
      @param[in] out_filename - string, filename to write
      @param[in] handler - instance of MooseInputFileRW
      @param[in] logger - python logger instance
      @return sim_data - modified python structure
  '''
  assert index_ig in [1, 2]
  sim = copy.deepcopy(sim_data)
  top_block_names = [elt['name'] for elt in sim['children']]

  ### Mesh
  if index_ig == 2:
    sim = __setMeshBlockFromFile(sim, '{0}.e'.format(SIM_IG1_NAME))

  ### AuxVariables
  if index_ig == 2:
    aux_vars_index = top_block_names.index('AuxVariables')
    aux_vars = sim['children'][aux_vars_index]
    all_aux_var_names = [elt['name'] for elt in aux_vars['children']]
    assert 'old_temp' not in all_aux_var_names
    assert 'u_diff_auxvar' not in all_aux_var_names
    aux_vars['children'].append({
      'name':'old_temp',
      'children':[],
      'comments':[],
      'attributes':[{'name':'initial_from_file_var','value':"temp", 'comment':''}],
    })
    aux_vars['children'].append({
      'name':'u_diff_auxvar',
      'children':[],
      'comments':[],
      'attributes':[],
    })
    #TODO: get proper name of variable

  ### Kernels
  # Disable timeDerivative kernels
  sim = __disableTimeDerivativeKernels(sim)

  ### Materials
  # Change material (just to let user see more easily that this parameter will get overwritten)
  materials_index = top_block_names.index('Materials')
  materials = sim['children'][materials_index]
  for material in materials['children']:
    material_name = material['name']
    for elt in material['attributes']:
      if elt['name']==CONT_PARAM_NAME:
        elt['value'] = '9999'
        elt['comment'] = 'will be overwritten by continuation wrapper'

  ### AuxKernels
  if index_ig == 2:
    auxkernels_index = top_block_names.index('AuxKernels')
    auxkernels = sim['children'][auxkernels_index]
    all_auxkernel_names = [elt['name'] for elt in auxkernels['children']]
    assert 'u_diff_auxkernel' not in all_auxkernel_names
    auxkernels['children'].append({
      'name':'u_diff_auxkernel',
      'children':[],
      'comments':[],
      'attributes':[
        {'name':'type','value':'RedbackDiffVarsAux', 'comment':''},
        {'name':'variable','value':'u_diff_auxvar', 'comment':''},
        {'name':'variable_2','value':'old_temp', 'comment':''},
        {'name':'variable_1','value':'temp', 'comment':''},
        {'name':'execute_on','value':'timestep_end', 'comment':''}
        ],
    })

  ### Postprocessors
  sim = __addPostProcessors(sim, index_ig==2)

  ### Executioner
  sim = __setExecutionerSteady(sim)

  ### Outputs
  base_filename = SIM_IG1_NAME
  if index_ig == 2:
    base_filename = SIM_IG2_NAME
  sim = __updateOutputs(sim, base_filename)
  outputs_index = top_block_names.index('Outputs')

  # write to file
  handler.write(sim, out_filename)
  return sim

def writeIterationFile(sim_data, out_filename, handler, logger, coeff_mutlipliers=(2,-1)):
  ''' Write iteration file
      @param[in] sim_data - python structure containing simulation data
      @param[in] filename - string, filename to write
      @param[in] logger - python logger instance
      @param[in] coeff_mutlipliers - list of 2 multipliers (mult_old, mult_older) 
        to apply to u_old and u_older such that the initial guess for the solution is
        u_guess = mult_old*u_old + mult_older*u_older
      @return sim_data - modified python structure
  '''
  sim = copy.deepcopy(sim_data)
  top_block_names = [elt['name'] for elt in sim['children']]

  ### Mesh
  sim = __setMeshBlockFromFile(sim, '{0}.e'.format(SIM_IG2_NAME))

  ### Variables
  variables_index = top_block_names.index('Variables')
  variables = sim['children'][variables_index]
  attribute_names = [attr['name'] for attr in variables['attributes']]
  active_index = None
  if 'active' in attribute_names:
    active_index = attribute_names.index('active')
  all_variables_names = [elt['name'] for elt in variables['children']]
  # find active variables
  if active_index is None:
    variables['attributes'].insert(0, {'name':'active','value':'', 'comment':''})
    active_index = 0
    active_variables_names = all_variables_names
  else:
    active_variables_names = __getListFromString(variables['attributes'][active_index]['value'])
  # Add initial conditions for all active variables
  for variable in variables['children']:
    if variable['name'] in active_variables_names:
      index_ic = None
      for index_child, child in enumerate(variable['children']):
        if child['name'] == 'InitialCondition':
          index_ic = index_child
      if index_ic is None:
        variable['children'].insert(0, {
          'name':'InitialCondition',
          'children':[],
          'comments':[],
          'attributes':[],
        })
        index_ic = 0
      variable['children'][index_ic]['attributes'] = \
        [{'name':'function','value':'initial_solution', 'comment':''},
         {'name':'type','value':'FunctionIC', 'comment':''},
         {'name':'variable','value':'{0}'.format(variable['name']), 'comment':''}]
  # add continuation variable
  cont_var_name = 'lambda'
  while cont_var_name in all_variables_names:
    cont_var_name += 'x'
  variables['children'].append({
    'name':'{0}'.format(cont_var_name),
    'children':[],
    'comments':[],
    'attributes':[{'name':'family','value':'SCALAR', 'comment':''},
                  {'name':'initial_condition','value':'9999', 'comment':'to be set up by the continuation wrapper'},],
  })
  # Update list of active variables
  variables['attributes'][active_index]['value'] = "'{0}'"\
    .format(' '.join(active_variables_names + [cont_var_name]))

  ### ICs
  # Remove all ICs
  if 'ICs' in top_block_names:
    ics_index = top_block_names.index('ICs')
    sim['children'].pop(ics_index)
    top_block_names.pop(ics_index)

  ### Global parameteres
  if 'GlobalParams' in top_block_names:
    global_params_index = top_block_names.index('GlobalParams')
  else:
    global_params_index = variables_index + 1
    sim['children'].insert(global_params_index, {
      'name':'GlobalParams',
      'children':[],
      'comments':[],
      'attributes':[],
    })
    top_block_names.insert(global_params_index, 'GlobalParams')
  global_params = sim['children'][global_params_index]
  for label in ['ds', 'ds_old']:
    assert label not in [elt['name'] for elt in global_params['attributes']]
  global_params['attributes'].extend([{'name':'ds','value':"9999", 'comment':'will be overwritten by continuation wrapper'},
                                      {'name':'ds_old','value':"9999", 'comment':'will be overwritten by continuation wrapper'}])

  ### AuxVariables
  if 'AuxVariables' in top_block_names:
    aux_vars_index = top_block_names.index('AuxVariables')
  else:
    aux_vars_index = kernels_index + 1
    sim['children'].insert(aux_vars_index, {
      'name':'AuxVariables',
      'children':[],
      'comments':[],
      'attributes':[],
    })
    top_block_names.insert(aux_vars_index, 'AuxVariables')
  aux_vars = sim['children'][aux_vars_index]
  all_aux_var_names = [elt['name'] for elt in aux_vars['children']]
  for label in ['old_temp', 'older_temp', 'directional_derivative', 'u_diff_auxvar']:
    assert label not in all_aux_var_names
  aux_vars['children'].append({
    'name':'old_temp',
    'children':[],
    'comments':[],
    'attributes':[{'name':'initial_from_file_var','value':"temp", 'comment':''}],
  })
  aux_vars['children'].append({
    'name':'older_temp',
    'children':[],
    'comments':[],
    'attributes':[{'name':'initial_from_file_var','value':"old_temp", 'comment':''}],
  })
  aux_vars['children'].append({
    'name':'directional_derivative',
    'children':[],
    'comments':[],
    'attributes':[{'name':'family','value':"SCALAR", 'comment':''}],
  })
  aux_vars['children'].append({
    'name':'u_diff_auxvar',
    'children':[],
    'comments':[],
    'attributes':[],
  })

  ### Functions
  if 'Functions' in top_block_names:
    functions_index = top_block_names.index('Functions')
  else:
    functions_index = aux_vars_index + 1
    sim['children'].insert(functions_index, {
      'name':'Functions',
      'children':[],
      'comments':[],
      'attributes':[],
    })
    top_block_names.insert(functions_index, 'Functions')
  functions = sim['children'][functions_index]
  all_functions_names = [elt['name'] for elt in functions['children']]
  for label in ['u_old', 'u_older', 'initial_solution']:
    assert label not in all_functions_names
  functions['children'].append({
    'name':'u_old',
    'children':[],
    'comments':[],
    'attributes':[{'name':'type','value':"SolutionFunction", 'comment':''},
                  {'name':'solution','value':"old_temp_UO", 'comment':''}],
  })
  functions['children'].append({
    'name':'u_older',
    'children':[],
    'comments':[],
    'attributes':[{'name':'type','value':"SolutionFunction", 'comment':''},
                  {'name':'solution','value':"older_temp_UO", 'comment':''}],
  })
  functions['children'].append({
    'name':'initial_solution',
    'children':[],
    'comments':[],
    'attributes':[{'name':'type','value':"LinearCombinationFunction", 'comment':''},
                  {'name':'functions','value':"'u_old u_older'", 'comment':''},
                  {'name':'w','value':"'{} {}'".\
                   format(coeff_mutlipliers[0], coeff_mutlipliers[1]), 'comment':''}],
  })

  ### Kernels
  # Disable timeDerivative kernels
  kernels_index = top_block_names.index('Kernels')
  sim = __disableTimeDerivativeKernels(sim)

  ### ScalarKernels
  if 'ScalarKernels' in top_block_names:
    scalarkernels_index = top_block_names.index('ScalarKernels')
  else:
    scalarkernels_index = kernels_index + 1
    sim['children'].insert(scalarkernels_index, {
      'name':'ScalarKernels',
      'children':[],
      'comments':[],
      'attributes':[],
    })
    top_block_names.insert(scalarkernels_index, 'ScalarKernels')
  scalarkernels = sim['children'][scalarkernels_index]
  all_scalarkernels_names = [elt['name'] for elt in scalarkernels['children']]
  assert 'continuation_kernel' not in all_scalarkernels_names
  scalarkernels['children'].append({
    'name':'continuation_kernel',
    'children':[],
    'comments':[],
    'attributes':[{'name':'type','value':'RedbackContinuation', 'comment':''},
                  {'name':'continuation_parameter_old','value':'9999', 'comment':'will be overwritten by continuation wrapper'},
                  {'name':'continuation_parameter_older','value':'9999', 'comment':'will be overwritten by continuation wrapper'},
                  {'name':'directional_derivative','value':'directional_derivative', 'comment':''},
                  {'name':'variable','value':'{0}'.format(cont_var_name), 'comment':''}],
  })

  ### AuxKernels
  if 'AuxKernels' in top_block_names:
    auxkernels_index = top_block_names.index('AuxKernels')
  else:
    auxkernels_index = scalarkernels_index + 1
    sim['children'].insert(auxkernels_index, {
      'name':'AuxKernels',
      'children':[],
      'comments':[],
      'attributes':[],
    })
    top_block_names.insert(auxkernels_index, 'AuxKernels')
  auxkernels = sim['children'][auxkernels_index]
  all_auxkernel_names = [elt['name'] for elt in auxkernels['children']]
  assert 'u_diff_auxkernel' not in all_auxkernel_names
  auxkernels['children'].append({
    'name':'u_diff_auxkernel',
    'children':[],
    'comments':[],
    'attributes':[
      {'name':'type','value':'RedbackDiffVarsAux', 'comment':''},
      {'name':'variable','value':'u_diff_auxvar', 'comment':''},
      {'name':'variable_2','value':'old_temp', 'comment':''},
      {'name':'variable_1','value':'temp', 'comment':''},
      {'name':'execute_on','value':'timestep_end', 'comment':''}
      ],
  })

  ### AuxScalarKernels
  if 'AuxScalarKernels' in top_block_names:
    auxscalarkernels_index = top_block_names.index('AuxScalarKernels')
  else:
    auxscalarkernels_index = auxkernels_index + 1
    sim['children'].insert(auxscalarkernels_index, {
      'name':'AuxScalarKernels',
      'children':[],
      'comments':[],
      'attributes':[],
    })
    top_block_names.insert(auxscalarkernels_index, 'AuxScalarKernels')
  auxscalarkernels = sim['children'][auxscalarkernels_index]
  assert 'directional_derivative' not in all_auxkernel_names
  auxscalarkernels['children'].append({
    'name':'directional_derivative',
    'children':[],
    'comments':[],
    'attributes':[
      {'name':'type','value':'RedbackContinuationTangentAux', 'comment':''},
      {'name':'variable','value':'directional_derivative', 'comment':''},
      {'name':'nodes','value':'9999', 'comment':'irrelevant value, overwritten by C++ code'},
      {'name':'sum_var_1','value':'temp', 'comment':''},
      {'name':'sum_var_old_1','value':'old_temp', 'comment':''},
      {'name':'sum_var_older_1','value':'older_temp', 'comment':''}
      ],
  })

  ### Materials
  materials_index = top_block_names.index('Materials')
  materials = sim['children'][materials_index]
  for material in materials['children']:
    mat_attr_names = [elt['name'] for elt in material['attributes']]
    type_index = mat_attr_names.index('type')
    if material['attributes'][type_index]['value'] == 'RedbackMaterial':
      # add continuation parameter
      if 'continuation_parameter' in mat_attr_names:
        index_cont_param = mat_attr_names.index('continuation_parameter')
      else:
        index_cont_param = type_index + 1
        material['attributes'].insert(index_cont_param, {'name':'continuation_parameter','value':'', 'comment':''})
        mat_attr_names.insert(index_cont_param, 'continuation_parameter')
      material['attributes'][index_cont_param]['value'] = cont_var_name
      # set value of cont_param_name to 1
      if CONT_PARAM_NAME in mat_attr_names:
        index_cont_param = mat_attr_names.index(CONT_PARAM_NAME)
      else:
        index_cont_param = index_cont_param + 1
        material['attributes'].insert\
          (index_cont_param, {'name':CONT_PARAM_NAME,'value':'', 'comment':''})
      material['attributes'][index_cont_param]['value'] = '1'
      material['attributes'][index_cont_param]['comment'] = \
        'Gets multiplied by value of scalar variable {0}'.format(cont_var_name)

  ### PostProcessors
  sim = __addPostProcessors(sim, add_l2_norm_diff=True)
  pps_index = top_block_names.index('Postprocessors')

  ### UserObjects
  if 'UserObjects' in top_block_names:
    userobjects_index = top_block_names.index('UserObjects')
  else:
    userobjects_index = pps_index + 1
    sim['children'].insert(userobjects_index, {
      'name':'UserObjects',
      'children':[],
      'comments':[],
      'attributes':[],
    })
    top_block_names.insert(userobjects_index, 'UserObjects')
  userobjects = sim['children'][userobjects_index]
  all_userobjects_names = [elt['name'] for elt in userobjects['children']]
  for label in ['old_temp_UO', 'older_temp_UO']:
    assert label not in all_userobjects_names
  userobjects['children'].append({
    'name':'old_temp_UO',
    'children':[],
    'comments':[],
    'attributes':[
      {'name':'type','value':'SolutionUserObject', 'comment':''},
      {'name':'timestep','value':'LATEST', 'comment':''},
      {'name':'system_variables','value':'temp', 'comment':''},
      {'name':'mesh','value':'extra_param_initial_guess2.e', 'comment':'will be overwritten by continuation wrapper'},
      {'name':'execute_on','value':'initial', 'comment':''}
      ],
  })
  userobjects['children'].append({
    'name':'older_temp_UO',
    'children':[],
    'comments':[],
    'attributes':[
      {'name':'type','value':'SolutionUserObject', 'comment':''},
      {'name':'timestep','value':'LATEST', 'comment':''},
      {'name':'system_variables','value':'old_temp', 'comment':''},
      {'name':'mesh','value':'extra_param_initial_guess2.e', 'comment':'will be overwritten by continuation wrapper'},
      {'name':'execute_on','value':'initial', 'comment':''}
      ],
  })

  ### Executioner
  sim = __setExecutionerSteady(sim)

  ### Outputs
  sim = __updateOutputs(sim, base_filename=SIM_ITER_NAME)

  # write to file
  handler.write(sim, out_filename)
  return sim

def updateMultiplyingCoefficientsForInitialGuess(sim_data, coeff_guess_old, coeff_guess_older):
  ''' Update multiplying coefficients for initial guess of the solution
      @param[in,out] sim_data - python structure with simulation data. Gets modified
      @return sim_data - updated python structure
  '''
  top_block_names = [elt['name'] for elt in sim_data['children']]
  functions_index = top_block_names.index('Functions')
  functions = sim_data['children'][functions_index]
  all_functions_names = [elt['name'] for elt in functions['children']]
  index_ig = all_functions_names.index('initial_solution')
  function_ig = functions['children'][index_ig]
  attr_names = [attr['name'] for attr in function_ig['attributes']]
  index_w = attr_names.index('w')
  function_ig['attributes'][index_w]['value'] = "'{0} {1}'".format(coeff_guess_old, coeff_guess_older)
  return sim_data

def __setMeshBlockFromFile(sim_data, mesh_filename):
  ''' Update parameters to set mesh from file with given filename
      @param[in,out] sim_data - python structure with simulation data. Gets modified
      @return sim_data - updated python structure
  '''
  mesh_index = [elt['name'] for elt in sim_data['children']].index('Mesh')
  mesh = sim_data['children'][mesh_index]
  index_file = -1
  for i_attr, attr in enumerate(mesh['attributes']):
    if attr['name'] == 'type':
      attr['value'] = 'FileMesh'
    elif attr['name'] == 'file':
      index_file = i_attr
  if index_file < 0:
    mesh['attributes'].insert(1, {'name':'file','value':'', 'comment':''})
    index_file = 1
  mesh['attributes'][index_file]['value'] = '{0}'.format(mesh_filename)
  # remove useless attributes
  attributes_to_remove = ['nx', 'ny', 'nz', 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax']
  indices_to_remove = []
  for i_attr, attr in enumerate(mesh['attributes']):
    if attr['name'] in attributes_to_remove:
      indices_to_remove.append(i_attr)
  indices_to_remove.reverse()
  for i in indices_to_remove:
    mesh['attributes'].pop(i)
  return sim_data

def __disableTimeDerivativeKernels(sim_data):
  ''' Disable time derivative kernels
      @param[in,out] sim_data - python structure with simulation data. Gets modified
      @return sim_data - updated python structure
  '''
  kernels_index = [elt['name'] for elt in sim_data['children']].index('Kernels')
  kernels_to_keep = []
  kernels = sim_data['children'][kernels_index]
  for kernel in kernels['children']:
    kernel_name = kernel['name']
    kernel_type = [elt['value'] for elt in kernel['attributes'] if elt['name']=='type'][0]
    if kernel_type != 'TimeDerivative':
      kernels_to_keep.append(kernel_name)
  attribute_names = [attr['name'] for attr in kernels['attributes']]
  active_index = None
  if 'active' in attribute_names:
    active_index = attribute_names.index('active')
  if active_index is None:
    kernels['attributes'].insert(0, {'name':'active','value':"''", 'comment':''})
    active_index = 0
  kernels['attributes'][active_index]['value'] = "'{0}'".format(' '.join(kernels_to_keep))
  return sim_data

def __addPostProcessors(sim_data, add_l2_norm_diff):
  ''' Add PostProcessors required for continuation algorithm
      @param[in,out] sim_data - python structure with simulation data. Gets modified
      @param[in] add_l2_norm_diff - boolean, True if we want to write L2_norm_diff postprocessor
      @return sim_data - updated python structure
  '''
  # Check if PostProcessors block exists
  top_block_names = [elt['name'] for elt in sim_data['children']]
  if 'Postprocessors' in top_block_names:
    pps_index = top_block_names.index('Postprocessors')
  else:
    materials_index = top_block_names.index('Materials')
    pps_index = materials_index + 1
    sim['children'].insert(pps_index, {
      'name':'Postprocessors',
      'children':[],
      'comments':[],
      'attributes':[],
    })
  pps = sim_data['children'][pps_index]
  index_pp_solution_norm = -1
  index_pp_l2norm_u = -1
  index_pp_l2norm_u_diff = -1
  all_pps_names = []
  for i_pp, pp in enumerate(pps['children']):
    all_pps_names.append(pp['name'])
    if pp['name'] == 'solution_norm':
      index_pp_solution_norm = i_pp
    if pp['name'] == 'L2_norm_u':
      index_pp_l2norm_u = i_pp
    if pp['name'] == 'L2_norm_u_diff':
      index_pp_l2norm_u_diff = i_pp
  # Add 'solution_norm'
  if index_pp_solution_norm < 0:
    pps['children'].insert(0, {'name':'solution_norm', 'children':[],
                               'attributes':[{'name':'type','value':'NodalMaxValue', 'comment':''},
                                             {'name':'variable','value':'temp', 'comment':''}],
                               # TODO: 'temp' should be replaced by actual variable name
                               'comments':['L_inf norm of solution for continuation algorithm']})
  else:
    # there is already a pp with that name, let's overwrite it just to be sure...
    pp_solution_norm = pps['children'][index_pp_solution_norm]
    pp_solution_norm['attributes'] = \
      [{'name':'type','value':'NodalMaxValue', 'comment':''},
       {'name':'variable','value':'temp', 'comment':''}]
    pp_solution_norm['comments'] = ['L_inf norm of solution for continuation algorithm']
  # Add 'L2_norm_u'
  if index_pp_l2norm_u < 0:
    pps['children'].insert(0, {'name':'L2_norm_u', 'children':[],
                               'attributes':[{'name':'type','value':'NodalL2Norm', 'comment':''},
                                             {'name':'variable','value':'temp', 'comment':''}],
                               # TODO: 'temp' should be replaced by actual variable name
                               'comments':['L2 norm of solution for continuation algorithm']})
  else:
    # there is already a pp with that name, let's overwrite it just to be sure...
    pp_sol_l2norm = pps['children'][index_pp_l2norm_u]
    pp_sol_l2norm['attributes'] = \
      [{'name':'type','value':'NodalL2Norm', 'comment':''},
       {'name':'variable','value':'temp', 'comment':''}]
    pp_sol_l2norm['comments'] = ['L2 norm of solution for continuation algorithm']
  # Add 'L2_norm_u_diff'
  if add_l2_norm_diff:
    if index_pp_l2norm_u_diff < 0:
      pps['children'].insert(0, {'name':'L2_norm_u_diff', 'children':[],
                                 'attributes':[{'name':'type','value':'NodalL2Norm', 'comment':''},
                                               {'name':'variable','value':'u_diff_auxvar', 'comment':''}],
                                 # TODO: 'temp' should be replaced by actual variable name
                                 'comments':['L2 norm of solution for continuation algorithm']})
    else:
      # there is already a pp with that name, let's overwrite it just to be sure...
      pp_sol_l2norm = pps['children'][index_pp_l2norm_u_diff]
      pp_sol_l2norm['attributes'] = \
        [{'name':'type','value':'NodalL2Norm', 'comment':''},
         {'name':'variable','value':'u_diff_auxvar', 'comment':''}]
      pp_sol_l2norm['comments'] = ['L2 norm of delta solution for continuation algorithm']
  # find active postprocessors
  attribute_names = [attr['name'] for attr in pps['attributes']]
  active_index = None
  if 'active' in attribute_names:
    active_index = attribute_names.index('active')
  all_pps_names = [elt['name'] for elt in pps['children']]
  if active_index is None:
    pps['attributes'].insert(0, {'name':'active','value':'', 'comment':''})
    active_index = 0
    active_pps_names = all_pps_names
  else:
    active_pps_names = __getListFromString(pps['attributes'][active_index]['value'])
  # update/create list of active post processors
  for label in ['solution_norm', 'L2_norm_u']:
    if label not in active_pps_names:
      active_pps_names.append(label)
  if add_l2_norm_diff:
    active_pps_names.append('L2_norm_u_diff')
  text = ' '.join(active_pps_names)
  if len(active_pps_names) > 1:
    text = "'{0}'".format(text)
  pps['attributes'][active_index]['value'] = text
  return sim_data

def __setExecutionerSteady(sim_data):
  ''' Set executioner type to steady.
      @param[in,out] sim_data - python structure with simulation data. Gets modified
      @return sim_data - updated python structure
  '''
  top_block_names = [elt['name'] for elt in sim_data['children']]
  exec_index = top_block_names.index('Executioner')
  executioner = sim_data['children'][exec_index]
  do_change = False
  for i_attr, attr in enumerate(executioner['attributes']):
    if attr['name'] == 'type':
      if attr['value'] == 'Steady':
        pass
      elif attr['value'] == 'Transient':
        do_change = True
        attr['value'] = 'Steady'
        break
  if do_change:
    attributes_to_remove = ['num_steps', 'ss_check_tol', 'end_time', 'dtmax', 'scheme'] # not exhaustive for now...
    children_to_remove = ['TimeStepper'] # not exhaustive for now...
    indices_to_remove = []
    for i_attr, attr in enumerate(executioner['attributes']):
      if attr['name'] in attributes_to_remove:
        indices_to_remove.append(i_attr)
    indices_to_remove.reverse()
    for i in indices_to_remove:
      executioner['attributes'].pop(i)
    indices_to_remove = []
    for i_child, child in enumerate(executioner['children']):
      if child['name'] in children_to_remove:
        indices_to_remove.append(i_child)
    indices_to_remove.reverse()
    for i in indices_to_remove:
      executioner['children'].pop(i)
  return sim_data

def __updateOutputs(sim_data, base_filename):
  ''' Update outputs
      @param[in,out] sim_data - python structure with simulation data. Gets modified
      @param[in] base_filename - string, base filename for outputs
      @return sim_data - updated python structure
  '''
  top_block_names = [elt['name'] for elt in sim_data['children']]
  outputs_index = top_block_names.index('Outputs')
  outputs = sim_data['children'][outputs_index]
  found_file_base = False
  found_exodus = False
  found_csv = False
  found_execute_on = False
  for i_attr, attr in enumerate(outputs['attributes']):
    if attr['name'] == 'file_base':
      attr['value'] = base_filename
      found_file_base = True
    elif attr['name'] == 'exodus':
      attr['value'] = 'true'
      found_exodus = True
    elif attr['name'] == 'csv':
      attr['value'] = 'true'
      found_csv = True
    elif attr['name'] == 'execute_on':
      attr['value'] = "'initial timestep_end'"
      found_execute_on = True
  if not found_file_base:
    outputs['attributes'].append({'name':'file_base','value':base_filename, 'comment':''})
  if not found_exodus:
    outputs['attributes'].append({'name':'exodus','value':'true', 'comment':''})
  if not found_csv:
    outputs['attributes'].append({'name':'csv','value':'true', 'comment':''})
  if not found_execute_on:
    outputs['attributes'].append({'name':'execute_on','value':"'initial timestep_end'", 'comment':''})
  return sim_data

def __getListFromString(text):
  ''' Get list of strings from text representing a list of items
      @param[in] text - string, text to parse
      @return list - python list of strings
  '''
  tmp = text.strip()
  if tmp.startswith("'") and tmp.endswith("'"):
    tmp = tmp[1:-1]
  list = tmp.split(' ')
  return list

if __name__ == "__main__":
  logger = getLogger('sim', 'running_tmp/log.txt', logging.INFO)
  handler = MooseInputFileRW()
  data_sim = handler.read('benchmark_1_T/bench1_a.i')
  sim_1 = writeInitialGuessFile(1, data_sim, 'running_tmp/extra_param_initial_guess1.i', handler, logger)
  sim_2 = writeInitialGuessFile(2, data_sim, 'running_tmp/extra_param_initial_guess2.i', handler, logger)
  sim_i = writeIterationFile(data_sim, 'running_tmp/extra_param_iteration.i', handler, logger)
  print 'Finished'
