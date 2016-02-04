function [data, dim, nb_vars, bifurc_param_name, error] = parseInputFile(filename)
    % Parse input parameter (XML) file and return data. The number of
    % variables is defined by the number of equations read.
    % * data: containers.Map containing:
    %   * params: Map with parameters (name:value)
    %   * numerical_params: Map with numerical parameters
    %   * equations: list of equations as strings (nb=nb_vars)
    %   * ics: list of initial conditions as strings (nb=nb_vars)
    %   * bcs: list of boundary conditions as strings (nb=nb_vars)
    %   * ppcs: list of postprocessing commands as strings
    % * dim: int, dimensionality of the problem
    % * nb_vars: int, number of variables (/equations)
    % * bifurc_param_name: string, name of bifurcation parameter
    % * error: int, flag equal to 1 if problem, 0 if parsing went OK 
    data = containers.Map;
    data('params') = containers.Map;
    data('numerical_params') = containers.Map;
    data('equations') = {};
    data('ics') = {};
    data('bcs') = {};
    data('ppcs') = {};
    dim = 0;
    nb_vars = 0;
    error = 1;
    % Check that input file exists
    %fprintf('Parsing "%s"\n',which(filename));
    if exist(filename) ~= 2
        fprintf('Error, file "%s" not found. Bye!\n',filename);
        return
    end
    %%% PARSE XML FILE %%%
    xdoc = xml2struct(filename);
    % Equation parameters
    bifurc_param_name = xdoc.Model.EquationParameters.bifurcationParameter.Text;
    nb_params = length(xdoc.Model.EquationParameters.param);
    params = containers.Map;
    if nb_params == 1
        name = strtrim(xdoc.Model.EquationParameters.param.name.Text);
        value = strtrim(xdoc.Model.EquationParameters.param.value.Text); % as string
        params(name) = value;
    else
        for i=1:nb_params
            name = strtrim(xdoc.Model.EquationParameters.param{i}.name.Text);
            value = strtrim(xdoc.Model.EquationParameters.param{i}.value.Text); % as string
            params(name) = value;
        end
    end
    data('params') = params;
    % Numerical parameters
    nb_numparams = length(xdoc.Model.NumericalParameters.param);
    numparams = containers.Map;
    if nb_numparams == 1
        name = strtrim(xdoc.Model.NumericalParameters.param.name.Text);
        value = strtrim(xdoc.Model.NumericalParameters.param.value.Text); % as string
        numparams(name) = value;
    else
        for i=1:nb_numparams
            name = strtrim(xdoc.Model.NumericalParameters.param{i}.name.Text);
            value = strtrim(xdoc.Model.NumericalParameters.param{i}.value.Text); % as string
            numparams(name) = value;
        end
    end
    data('numerical_params') = numparams;
    % Equations
    dim = str2num(xdoc.Model.Equations.dim.Text);
    nb_vars = length(xdoc.Model.Equations.EQ);
    var_ids = [];
    eqs = {};
    if nb_vars == 1
        var_id = str2num(xdoc.Model.Equations.EQ.Attributes.var);
        var_ids = [var_ids; var_id];
        eqs{1} = strtrim(xdoc.Model.Equations.EQ.Text);
    else
        for i=1:nb_vars
            var_id = str2num(xdoc.Model.Equations.EQ{i}.Attributes.var);
            var_ids = [var_ids; var_id];
            eqs{i} = strtrim(xdoc.Model.Equations.EQ{i}.Text);
        end
    end
    [dummy ordered_ids] = sort(var_ids);
    data('equations') = eqs(ordered_ids);
    % Boundary conditions
    nb_bcs = length(xdoc.Model.BoundaryConditions.BC);
    var_ids = [];
    bcs = {};
    if nb_bcs == 1
        var_id = str2num(xdoc.Model.BoundaryConditions.BC.Attributes.var);
        var_ids = [var_ids; var_id];
        bcs{1} = strtrim(xdoc.Model.BoundaryConditions.BC.Text);
    else
        for i=1:nb_bcs
            var_id = str2num(xdoc.Model.BoundaryConditions.BC{i}.Attributes.var);
            var_ids = [var_ids; var_id];
            bcs{i} = strtrim(xdoc.Model.BoundaryConditions.BC{i}.Text);
        end
    end
    [dummy ordered_ids] = sort(var_ids);
    data('bcs') = bcs(ordered_ids);
    % Initial conditions
    nb_ics = length(xdoc.Model.InitialConditions.IC);
    var_ids = [];
    ics = {};
    if nb_ics == 1
        var_id = str2num(xdoc.Model.InitialConditions.IC.Attributes.var);
        var_ids = [var_ids; var_id];
        ics{1} = strtrim(xdoc.Model.InitialConditions.IC.Text);
    else
        for i=1:nb_ics
            var_id = str2num(xdoc.Model.InitialConditions.IC{i}.Attributes.var);
            var_ids = [var_ids; var_id];
            ics{i} = strtrim(xdoc.Model.InitialConditions.IC{i}.Text);
        end
    end
    [dummy ordered_ids] = sort(var_ids);
    data('ics') = ics(ordered_ids);
    % Preprocessing commands
    nb_ppcs = length(xdoc.Model.PreProcessingCode.command);
    var_ids = [];
    ppcs = {};
    if nb_ppcs == 1
        var_id = str2num(xdoc.Model.PreProcessingCode.command.Attributes.index);
        var_ids = [var_ids; var_id];
        ppcs{1} = strtrim(xdoc.Model.PreProcessingCode.command.Text);
    else
        for i=1:nb_ppcs
            var_id = str2num(xdoc.Model.PreProcessingCode.command{i}.Attributes.index);
            var_ids = [var_ids; var_id];
            ppcs{i} = strtrim(xdoc.Model.PreProcessingCode.command{i}.Text);
        end
    end
    [dummy ordered_ids] = sort(var_ids);
    data('ppcs') = ppcs(ordered_ids);
    
    error = 0; % Parsing went OK
    return
end