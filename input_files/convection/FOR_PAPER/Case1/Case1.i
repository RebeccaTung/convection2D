[GlobalParams]
  time_factor = 1.e-3
[]

[Mesh]
  type = FileMesh
  file = 2d.msh
  dim = 2
[]

[Variables]
  [./temp]
  [../]
  [./pore_pressure]
  [../]
[]

[Materials]
  [./redback_materialA]
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
    ref_lewis_nb = 3.30e-08
    solid_compressibility = 0.001
    solid_thermal_expansion = 1e-05
    total_porosity = total_porosity
  [../]
[]

[Functions]
  [./init_gradient_T]
    # (0.0-y*(1.0-0.0)*1000/500) + 0.2*1/2*(cos(pi*(2*y-(-0.5)-0)/(0-(-0.5)))+1)*cos(pi*(2*x-0-1)/(1-0))
    # 
    # T_max + (y-y_min)*(T_min-T_max)/(y_max-y_min)
    # + amplitude*1/2*(cos(pi*(2*y-y_min-y_max)/(y_max-y_min))+1)*cos(pi*(2*x-x_min-x_max)/(x_max-x_min))
    # 
    type = ParsedFunction
    value = 'T_max + (y-y_min)*(T_min-T_max)/(y_max-y_min) + amplitude*1/2*(cos(pi*(2*y-y_min-y_max)/(y_max-y_min))+1)*cos(pi*(2*x-x_min-x_max)/(x_max-x_min))'
    vals = '-1            0          0          2          0           1           0.2'
    vars = 'y_min y_max x_min x_max  T_min T_max amplitude'
  [../]
  [./init_gradient_P]
    type = ParsedFunction
    value = (0.02-1.96*y)
  [../]
  [./timestep_function]
    type = ParsedFunction
    value = 'min(max(1e-15, dt*max(0.2, 1-0.05*(n_li-50))), (1e-1)*50.0/max(abs(v_max), abs(v_min)))'
    vals = 'num_li num_nli min_fluid_vel_y max_fluid_vel_y dt'
    vars = 'n_li n_nli v_min v_max dt'
  [../]
  [./New_Nusselt_calc_fn]
    # sqrt(max(abs(x_max),abs(x_min))^2+max(abs(y_max),abs(y_min))^2)
    # max_gradT_x max_gradT_y
    type = ParsedFunction
    value = sqrt(max(abs(x_max),abs(x_min))^2+max(abs(y_max),abs(y_min))^2)
    vals = 'max_gradT_x max_gradT_y min_gradT_x min_gradT_y'
    vars = 'x_max y_max x_min y_min'
  [../]
  [./grad_for_Nusselt]
    type = ParsedFunction
    value = (max_norm_grad_T)/((T_bottom-T_top)/height)
    vals = 'max_norm_grad_T  1               0        1'
    vars = 'max_norm_grad_T T_bottom T_top height'
  [../]
[]

[BCs]
  active = 'temperature_bottom temperature_top top_corners_p'
  [./temperature_top]
    type = DirichletBC
    variable = temp
    boundary = 2
    value = 0.0
  [../]
  [./temperature_bottom]
    type = DirichletBC
    variable = temp
    boundary = 0
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
    boundary = '21 22'
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
  [./fluid_vel_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./grad_temp_x_var]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./grad_temp_y_var]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./inv_Le_perturb]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./norm_grad_T]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  active = 'pres_conv press_td temp_diff temp_td temp_conv press_diff'
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
  [./fluid_vel_x]
    type = MaterialRealVectorValueAux
    variable = fluid_vel_x
    property = fluid_velocity
  [../]
  [./grad_temp_x_kernel]
    type = VariableGradientComponent
    variable = grad_temp_x_var
    component = x
    gradient_variable = temp
  [../]
  [./grad_temp_y_kernel]
    type = VariableGradientComponent
    variable = grad_temp_y_var
    component = y
    gradient_variable = temp
  [../]
  [./grad_T]
    type = VectorMagnitudeAux
    variable = norm_grad_T
    y = grad_temp_y_var
    x = grad_temp_x_var
    execute_on = timestep_end
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
  active = 'num_nli max_gradT_y max_gradT_x min_fluid_vel_y num_li min_gradT_y min_gradT_x max_norm_grad_T Bec_Nusselt_new middle_porosity middle_press new_timestep max_fluid_vel_y norm_grad_T dt New_Nusselt_postproc middle_temp'
  [./middle_temp]
    type = PointValue
    variable = temp
    point = '1.0 -0.50 0'
  [../]
  [./middle_press]
    type = PointValue
    variable = pore_pressure
    point = '1.0 -0.50 0'
  [../]
  [./middle_porosity]
    type = PointValue
    variable = total_porosity
    point = '1.0 -0.50 0'
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
    variable = 'grad_temp_x_var grad_temp_y_var'
    boundary = 4
  [../]
  [./max_gradT_x]
    type = ElementExtremeValue
    variable = grad_temp_x_var
  [../]
  [./max_gradT_y]
    type = ElementExtremeValue
    variable = grad_temp_y_var
  [../]
  [./New_Nusselt_postproc]
    type = FunctionValuePostprocessor
    function = New_Nusselt_calc_fn
  [../]
  [./min_gradT_x]
    type = ElementExtremeValue
    variable = grad_temp_x_var
    value_type = min
  [../]
  [./min_gradT_y]
    type = ElementExtremeValue
    variable = grad_temp_y_var
    value_type = min
  [../]
  [./Bec_Nusselt_new]
    type = FunctionValuePostprocessor
    function = grad_for_Nusselt
  [../]
  [./norm_grad_T]
    type = ElementL2Norm
    variable = norm_grad_T
    execute_on = timestep_end
  [../]
  [./max_norm_grad_T]
    type = ElementExtremeValue
    variable = norm_grad_T
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 4000
  l_max_its = 200
  nl_max_its = 10
  solve_type = PJFNK
  petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -ksp_gmres_restart '
  petsc_options_value = 'gmres asm lu 201'
  nl_abs_tol = 1e-09 # 1e-10 to begin with
  nl_rel_tol = 1e-07 # 1e-10 to begin with
  line_search = basic
[]

[Outputs]
  file_base = case1
  print_linear_residuals = false
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
  [./export_all_in_exo]
    output_material_properties = true
    file_base = case1
    type = Exodus
    elemental_as_nodal = true
  [../]
[]

[ICs]
  active = 'IC_temp IC_pressure'
  [./IC_temp]
    function = init_gradient_T
    variable = temp
    type = FunctionIC
  [../]
  [./IC_pressure]
    function = init_gradient_P
    variable = pore_pressure
    type = FunctionIC
  [../]
  [./inv_Le_randomIC]
    # Inverse of Lewis number perturbation such that
    # 1/Le = 1/Le + 1/Le_perturb
    variable = inv_Le_perturb
    standard_deviation = 0.5
    type = FunctionLogNormalDistributionIC
    mean = 2e+07
  [../]
[]

