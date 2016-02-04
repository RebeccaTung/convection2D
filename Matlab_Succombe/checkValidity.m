function error = checkValidity(data, stype, dim, nb_vars, bifurc_param_name)
    % Check validity of input parameters
    % Still TODO:
    % * check all keys per category
    % * check types, value ranges...
    
    error = 1; % Return 1 if error, 0 if all parameters are valid
    
    % Check all categories are defined
    map = containers.Map(); % map between data attributes and xml keywords
    map('bcs') = 'BoundaryConditions';
    map('equations') = 'Equations';
    map('ics') = 'InitialConditions';
    map('numerical_params') = 'NumericalParameters';
    map('params') = 'EquationParameters';
    map('ppcs') = 'PreProcessingCode';
    expected_categories = ...
        {'bcs','equations','ics','numerical_params','params','ppcs'};
    categories = data.keys();
    for i=1:length(categories)
        category = categories{i};
        if isempty(strmatch(category, expected_categories,'exact'))
            fprintf('Missing input category "%s"\n',map(category))
            return
        end
    end
    
    % Check variable names used don't interfere with each other
    params = data('params');
    param_names = params.keys();
    for i=1:length(param_names)
        param_name = param_names{i};
        indices = 1:length(param_names);
        indices(i) = [];
        other_param_names = param_names(indices);
        for j=1:length(other_param_names)
            other_param_name = other_param_names{j};
            if ~isempty(strfind(other_param_name, param_name))
                fprintf('Parameter name "%s" is part of other parameter name "%s"\n',...
                    param_name, other_param_name);
                return
            end
        end
    end
    
    % Check variables names used don't interfere with reserved keywords
    reserved_keywords = {'D1', 'D2', 'L'};
    params = data('params');
    param_names = params.keys();
    for i=1:length(param_names)
        param_name = param_names{i};
        if ~isempty(strmatch(param_name, reserved_keywords,'exact'))
            fprintf(['Sorry, "%s" is a reserved keyword so you cannot ',...
                'use it as parameter name\n'],param_name)
            fprintf('Reserved keywords are: "%s"\n',...
                strjoin(reserved_keywords, '", "'))
            return
        end
    end
    
    % All checks passed, all good then
    error = 0;
    return
end

