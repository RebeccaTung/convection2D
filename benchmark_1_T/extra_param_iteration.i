[Mesh]
  type = FileMesh
  file = extra_param_initial_guess2.e
  dim = 1
[]

[Variables]
  [./temp]
    # 0.03276648
    [./InitialCondition]
      function = initial_solution
      type = FunctionIC
      variable = temp
    [../]
  [../]
  [./lambda]
    family = SCALAR
    initial_condition = 9999 # to be set up by the continuation wrapper
  [../]
[]

[GlobalParams]
  ds = 99999 # will be overwritten by continuation wrapper
  ds_old = 99999 # will be overwritten by continuation wrapper
[]

[AuxVariables]
  [./total_porosity]
    order = FIRST
    family = MONOMIAL
  [../]
  [./old_temp]
    initial_from_file_var = temp
  [../]
  [./older_temp]
    initial_from_file_var = old_temp
  [../]
  [./directional_derivative]
    family = SCALAR
  [../]
[]

[Functions]
  [./u_old]
    type = SolutionFunction
    solution = old_temp_UO
  [../]
  [./u_older]
    type = SolutionFunction
    solution = older_temp_UO
  [../]
  [./initial_solution]
    type = LinearCombinationFunction
    functions = 'u_old u_older'
    w = '2  -1'
  [../]
[]

[Kernels]
  active = 'diff_temp mh_temp'
  [./td_temp]
    type = TimeDerivative
    variable = temp
  [../]
  [./diff_temp]
    type = Diffusion
    variable = temp
  [../]
  [./mh_temp]
    type = RedbackMechDissip
    variable = temp
  [../]
[]

[AuxKernels]
  [./total_porosity]
    type = RedbackTotalPorosityAux
    variable = total_porosity
  [../]
[]

[AuxScalarKernels]
  [./directional_derivative]
    type = RedbackContinuationTangentAux
    variable = directional_derivative
    nodes = '0 1 2 3 4 5 6 7 8 9 10'
    sum_var_1 = temp
    sum_var_old_1 = old_temp
    sum_var_older_1 = older_temp
  [../]
[]

[BCs]
  active = 'left_temp right_temp'
  [./left_temp]
    type = DirichletBC
    variable = temp
    boundary = left
    value = 0
  [../]
  [./right_temp]
    type = DirichletBC
    variable = temp
    boundary = right
    value = 0
  [../]
  [./disp_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'left right'
    value = 0
  [../]
  [./disp_x_left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 1
  [../]
  [./disp_x_rigth]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0
  [../]
[]

[Materials]
  [./adim_rock]
    type = RedbackMaterial
    block = 0
    alpha_2 = 1
    ar = 10
    temperature = temp
    ar_F = 40
    ar_R = 1
    phi0 = 0.1
    da_endo = 1
    gr = 1
    total_porosity = total_porosity
    continuation_parameter = lambda
  [../]
[]

[Postprocessors]
  active = 'temp_pt_1 temp_pt_0 temp_pt_3 temp_pt_2 temp_pt_5 temp_pt_4 max_temp'
  [./middle_temp]
    type = PointValue
    variable = temp
    point = '0 0 0'
  [../]
  [./strain]
    type = StrainRatePoint
    variable = temp
    point = '0 0 0'
  [../]
  [./temp_pt_0]
    type = PointValue
    variable = temp
    point = '-1 0 0'
  [../]
  [./temp_pt_1]
    type = PointValue
    variable = temp
    point = '-0.8 0 0'
  [../]
  [./temp_pt_2]
    type = PointValue
    variable = temp
    point = '-0.6 0 0'
  [../]
  [./temp_pt_3]
    type = PointValue
    variable = temp
    point = '-0.4 0 0'
  [../]
  [./temp_pt_4]
    type = PointValue
    variable = temp
    point = '-0.2 0 0'
  [../]
  [./temp_pt_5]
    type = PointValue
    variable = temp
    point = '0 0 0'
  [../]
  [./max_temp]
    type = NodalMaxValue
    variable = temp
  [../]
[]

[UserObjects]
  [./old_temp_UO]
    type = SolutionUserObject
    timestep = LATEST
    system_variables = temp
    mesh = extra_param_initial_guess2.e
    execute_on = initial
  [../]
  [./older_temp_UO]
    type = SolutionUserObject
    timestep = LATEST
    system_variables = old_temp
    mesh = extra_param_initial_guess2.e
    execute_on = initial
  [../]
[]

[Executioner]
  type = Steady
  #petsc_options_iname = '-pc_type -pc_hypre_type'
  #petsc_options_value = 'hypre boomeramg'

  l_tol = 1e-5
  l_max_its = 1000
  nl_rel_tol = 1e-5
  nl_max_its = 50
  nl_abs_tol = 1e-8 # 1e-50


  #l_abs_step_tol = 1e-10 # -1
  #nl_rel_step_tol = 1e-5 # 1e-50
  #nl_abs_step_tol = 1e-10 # 1e-50
[]

[Preconditioning]
  active = 'SMP'
  [./SMP]
    petsc_options = '-snes_monitor -snes_linesearch_monitor -ksp_monitor'
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -ksp_max_it -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres     asm      1E-14        1E-14      200         500          lu          NONZERO'
    type = SMP
    full = true
  [../]
  [./smp_precond]
    type = SMP
    full = true
    #solve_type = NEWTON
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -ksp_atol -ksp_rtol'
    petsc_options_value = 'gmres     bjacobi   1E-4      1E-5       10           1E-4       1E-5'
  [../]
[]

[Outputs]
  file_base = extra_param_iteration
  exodus = true
  csv = true
  execute_on = 'initial timestep_end'
[]

[ScalarKernels]
  [./continuation_kernel]
    continuation_parameter_old = 99999 # overwritten by continuation wrapper
    continuation_parameter_older = 99999 # overwritten by continuation wrapper
    directional_derivative = directional_derivative
    variable = lambda
    type = RedbackContinuation
  [../]
[]
