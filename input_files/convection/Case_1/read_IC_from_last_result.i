[GlobalParams]
  time_factor = 1.e-3
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 10
  xmin = 0
  xmax = 1
  ymin = -0.5
  ymax = 0
[]

[MeshModifiers]
  [./top_left]
    type = AddExtraNodeset
    new_boundary = 98
    coord = '0 0'
  [../]
  [./top_right]
    type = AddExtraNodeset
    new_boundary = 99
    coord = '1 0'
  [../]
[]

[Variables]
  [./temp]
    [./InitialCondition]
      function = init_temp
      type = FunctionIC
      variable = temp
    [../]
  [../]
  [./pore_pressure]
    [./InitialCondition]
      function = init_p
      type = FunctionIC
      variable = pore_pressure
    [../]
  [../]
[]

[Materials]
  [./redback_nomech]
    type = RedbackMaterial
    block = 0
    temperature = temp
    pore_pres = pore_pressure
    Aphi = 1
    ar = 10
    ar_F = 1
    ar_R = 1
    are_convective_terms_on = true
    delta = 0.0333333333333
    eta1 = 1e3
    fluid_compressibility = 0.002
    fluid_thermal_expansion = 0.0007
    gravity = '0 -1.96 0'
    Peclet_number = 1.0
    phi0 = 0.3
    pressurization_coefficient = 0.166923076923
    ref_lewis_nb = 2.8e-08
    solid_compressibility = 0.001
    solid_thermal_expansion = 1e-05
    total_porosity = total_porosity
  [../]
[]

[Functions]
  [./timestep_function]
    type = ParsedFunction
    value = 'min(max(1e-15, dt*max(0.2, 1-0.05*(n_li-50))), (1e-1)*50.0/max(abs(v_max), abs(v_min)))'
    vals = 'num_li num_nli min_fluid_vel_y max_fluid_vel_y dt'
    vars = 'n_li n_nli v_min v_max dt'
  [../]
  [./init_p]
    type = SolutionFunction
    from_variable = pore_pressure
    solution = UO_p
  [../]
  [./init_temp]
    type = SolutionFunction
    from_variable = temp
    solution = UO_temp
  [../]
[]

[BCs]
  active = 'temperature_bottom temperature_top top_corners_p'
  [./temperature_top]
    type = DirichletBC
    variable = temp
    boundary = top
    value = 0.0
  [../]
  [./temperature_bottom]
    type = DirichletBC
    variable = temp
    boundary = bottom
    value = 1.0
  [../]
  [./pore_pressure_top]
    type = DirichletBC
    variable = pore_pressure
    boundary = top
    value = 0.02
  [../]
  [./top_corners_p]
    type = DirichletBC
    variable = pore_pressure
    boundary = '98 99'
    value = 0
  [../]
[]

[AuxVariables]
  [./total_porosity]
    family = MONOMIAL
  [../]
  [./Lewis_number]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fluid_vel_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./grad_temp]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  active = 'pres_conv temp_diff temp_conv press_diff'
  [./temp_td]
    type = TimeDerivative
    variable = temp
  [../]
  [./temp_diff]
    type = RedbackThermalDiffusion
    variable = temp
  [../]
  [./temp_conv]
    type = RedbackThermalConvection
    variable = temp
  [../]
  [./press_td]
    type = TimeDerivative
    variable = pore_pressure
  [../]
  [./press_diff]
    type = RedbackMassDiffusion
    variable = pore_pressure
  [../]
  [./pres_conv]
    type = RedbackMassConvection
    variable = pore_pressure
    temperature = temp
  [../]
  [./press_thermPress]
    type = RedbackThermalPressurization
    variable = pore_pressure
    temperature = temp
  [../]
[]

[AuxKernels]
  [./total_porosity]
    type = RedbackTotalPorosityAux
    variable = total_porosity
    mechanical_porosity = 0
    execute_on = linear
  [../]
  [./Lewis_number]
    type = MaterialRealAux
    variable = Lewis_number
    property = lewis_number
  [../]
  [./fluid_vel_y]
    type = MaterialRealVectorValueAux
    component = 1
    variable = fluid_vel_y
    property = fluid_velocity
  [../]
  [./grad_temp]
    type = VariableGradientComponent
    variable = grad_temp
    component = y
    block = 0
    execute_on = timestep_end
    boundary = bottom
    gradient_variable = temp
  [../]
[]

[UserObjects]
  [./UO_temp]
    type = SolutionUserObject
    timestep = LATEST
    system_variables = temp
    mesh = last_result.e
    execute_on = initial
  [../]
  [./UO_p]
    type = SolutionUserObject
    timestep = LATEST
    system_variables = pore_pressure
    mesh = last_result.e
    execute_on = initial
  [../]
[]

[Preconditioning]
  # active = ''
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./middle_temp]
    type = PointValue
    variable = temp
    point = '0.5 -0.25 0'
  [../]
  [./middle_press]
    type = PointValue
    variable = pore_pressure
    point = '0.5 -0.25 0'
  [../]
  [./middle_porosity]
    type = PointValue
    variable = total_porosity
    point = '0.5 -0.25 0'
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./num_li]
    type = NumLinearIterations
  [../]
  [./num_nli]
    type = NumNonlinearIterations
  [../]
  [./new_timestep]
    type = FunctionValuePostprocessor
    function = timestep_function
  [../]
  [./max_fluid_vel_y]
    type = ElementExtremeValue
    variable = fluid_vel_y
    execute_on = nonlinear
    value_type = max
  [../]
  [./min_fluid_vel_y]
    type = ElementExtremeValue
    variable = fluid_vel_y
    execute_on = nonlinear
    value_type = min
  [../]
  [./Nusselt_number]
    type = DifferencePostprocessor
    value1 = 1
    value2 = one_minus_Nusselt
  [../]
  [./one_minus_Nusselt]
    type = SideIntegralVariablePostprocessor
    variable = grad_temp
    boundary = bottom
  [../]
[]

[Executioner]
  type = Transient
  l_max_its = 200
  nl_max_its = 10
  solve_type = PJFNK
  petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -ksp_gmres_restart '
  petsc_options_value = 'gmres asm lu 201'
  nl_abs_tol = 1e-9 # 1e-10 to begin with
  line_search = basic
[]

[Outputs]
  file_base = bex_test
  exodus = true
  csv = true
  execute_on = 'initial timestep_end'
  print_linear_residuals = false
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]

