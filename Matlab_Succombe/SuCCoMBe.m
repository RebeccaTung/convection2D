function SuCCoMBe(input_filename, type)
    % Spectral Coupled multiphysiCs and Mechanics Bifurcation mEthod: pseudo arclength continuation
    % Thomas POULET, Manolis VEVEAKIS
    % CSIRO 
    % Contributions from
    % Sotiris Alevizos (NTUA, Greece), Aubert LeBrozec (EP, France)
    % 2013
    
    % input_filename (string) - name of XML file with all parameters
    % type (string) - 'steady-state' or 'transient', simulation type
    % Example:  SuCCoMBe('input_files/input1D.xml', 'ss')
    
    global nb_vars dim
    
    % Check simulation type requested
    stype = getStandardType(type);
    if strcmp(stype,'invalid')
        fprintf('Wrong input type "%s"\nBye\n',type)
        return
    end
    % Parse input file
    [data, dim, nb_vars, bifurc_param_name, error] = parseInputFile(input_filename); % Parse input file for parameters
    if error
        fprintf('Error parsing input "%s"\nBye\n',input_filename)
        return
    end
    % Check validity of input parameters
    validation_error = checkValidity(data, stype, dim, nb_vars, bifurc_param_name);
    if validation_error
        fprintf('Modify your input and try again. Bye\n')
        return
    end
    % Run appropriate code
    if strcmp(stype,'ss') % steady-state
        SuCCoMBe_S(data, bifurc_param_name);
    elseif strcmp(stype,'t') % transient
        SuCCoMBe_T(data);
    else
        fprintf('Wrong input type "%s". Bye\n',type)
        return
    end
end

function standard_type = getStandardType(type)
    % Get standard type "ss" (for steady state) or "t" (for transient)
    if strcmp(lower(type),'steady-state') || strcmp(lower(type),'ss') ||...
            strcmp(lower(type),'steady state')|| strcmp(lower(type),'steady')
        standard_type = 'ss';
    elseif strcmp(lower(type),'transient') || strcmp(lower(type),'t')
        standard_type = 't';
    else
        standard_type = 'invalid';
    end
end

function SuCCoMBe_T(data)
    % Run transient simulation
    % Note: could not find better way than global variables for ode23s
    global N nb_vars bcs equations ppcs dim indices_interior boundary_mask
    global D1 D2 L
    
    if dim > 1
        TODO % not implemented yet...
    end
    num_params = data('numerical_params');
    restartfilename = num_params('restartFileNameT'); % '' => start from beginning (no restart)
    if restartfilename
        TODO % not implemented yet...
        N_new = N;
        [N_old,T0,Ti,P0,Pi,Le,n,Gr,T,P,Tb,Tc,alpha,delta,mu,...
            Ar,Kc,deltaE,phi0,depth_km] = ...
            readParametersFromFile(restartfilename,start_time);
        N = N_new;
        if N ~= N_old
            % interpolate to get initial solution on proper base
            x_old = cos([0:N_old]*pi/N_old)'; % Chebyshev collocation points
            ini_T = bcinterp(w, x_old, T, x);
            ini_P = bcinterp(w, x_old, P, x);
        else
            ini_T = T;
            ini_P = P;
        end
    else
        % Running new simulation
        
        % Get numerical parameters
        N = str2num(num_params('N')); % Number of collocation points in [-1,1], must be even.
        savefilename = num_params('saveFileNameT');
        start_time = str2num(num_params('startTime'));
        final_time = str2num(num_params('finalTime'));
        solver_rel_tol = str2num(num_params('solverRelTol'));
        i_param_plot = str2num(num_params('indexVarPlot')); % index of parameter to plot
        
        equations = {}; % equations as strings ready to be evaluated
        bcs = data('bcs');
        eqs = data('equations');
        ppcs = data('ppcs'); % post processing commands
        params = data('params');
        param_names = params.keys;
        param_values = params.values;
        for i = 1:nb_vars
            equations{i} = eqs(i);
            % Replace names of parameters including bifurcation parameter
            for j =1:length(params)
                equations{i} = strrep(equations{i}, param_names{j}, param_values{j});
            end
        end
        for i = 1:length(ppcs)
            % Replace names of parameters including bifurcation parameter
            for j =1:length(params)
                ppcs{i} = strrep(ppcs{i}, param_names{j}, param_values{j});
            end
        end
        
        if start_time ~= 0
            start_time = 0;
            fprintf('Warning: start time taken as "0" (ignoring your input "%s")\n',...
                num_params('startTime'));
        end
    end
    
    % Run simulation
    tic;
    fprintf('Starting transient simulation on %s...\n', datestr(now));
    [export_dir,~,~] = fileparts(savefilename);
    if exist(export_dir) ~= 7
        mkdir(export_dir)
        disp(sprintf('\nCreated export directory "%s"', export_dir));
    end
    % Get index of points defining the boundaries of the domain
    [indices_interior, boundary_mask] = getBoundaries(dim, N);
    % Get Chebyshev points and weight factors
    [z,w] = getChebyshevPointsAndWeightFactors(N);
    % Differentiation matrix
    [D1,D2] = bcmatrix(w,z);
    L = 0; % just to be defined if dim==1
    if dim == 2
        I = eye(N+1);
        L = kron(I,D2) + kron(D2,I); % 2D Laplacian
    end
    % Solver options
    solv_opt = odeset('RelTol', solver_rel_tol);
    % Initial conditions
    ics = data('ics');
    initial_condition = zeros(1,nb_vars*(N-1)^dim); % Not including boundary values (Dirichlet BC)
    for i = 1:nb_vars
        initial_condition((i-1)*(N-1)^dim+1:i*(N-1)^dim) = eval(ics{i});
    end
    [times, total_sol] = ode23s(@syspde, [start_time,final_time], initial_condition, solv_opt);
    disp('Finished calculation')
    
    % Save to file
    u = {}; % Individual solutions
    [nb_lines, ~] = size(total_sol);
    for i = 1:nb_vars
        tmp = eval(bcs{i}).*ones(nb_lines,(N+1)^dim).*repmat(boundary_mask,nb_lines,1);
        tmp(:,indices_interior) = total_sol(:,(i-1)*(N-1)^dim+1:i*(N-1)^dim);
        u{i} = tmp;
    end
    save(savefilename, 'data', 'dim', 'nb_vars', 'times', 'u');
    disp(sprintf('Time elapsed: %.2f seconds',toc))
    
    % Plot results
    plotTransientResults(dim, N, times, u{i_param_plot}, i_param_plot);
end

function plotTransientResults(dim, N, times, values, i_param_plot)
    % Plot results from transient simulation
    k_centre = getSpatialCentreIndex(dim, N);
    figure('name', 'SuCCoMBe (transient)');
    plot(times, values(:,k_centre));
    title(sprintf('Time evolution of parameter %d', i_param_plot));
    xlabel('Normalised time');
    ylabel(sprintf('Parameter %d', i_param_plot));
end

function k_centre = getSpatialCentreIndex(dim, N)
    % Get index of the centre of the spatial domain for a line/column
    % vector of all spatial coordinates
    if dim == 1
        k_centre = N/2+1; % Spatial index of interval centre
    elseif dim == 2
        %k_centre = (N/2-1)*(N-1)+N/2; % Spatial index of interval centre
        k_centre = (N+2)*(N/2)+1; % Spatial index of interval centre
    else
        TODO % dim > 2 not implemented yet...
    end
end

function [N_old,T0,Ti,P0,Pi,Le,n,Gr,start_time,T,P,Tb,Tc,alpha,delta,mu,Ar,Kc,deltaE,phi0,depth_km] ...
    = readInitialParamsFromFile(initial_file,k)
    % Load parameters from restart file for transient simulation
    TODO
	disp(sprintf('\nLoading %s', initial_file));
    load(initial_file, 'times','eps_dot_centre','T_centre',...
        'P_centre','N','Ps','Ts','eps_dot','Ti','T0','Pi','P0','Le','n',...
        'Gr','Tb','Tc','alpha','delta','mu','Ar','Kc','deltaE','phi0','depth_km');
    N_old = N;
    TODO %s find k so that start_time = times(k);
    disp(sprintf('Restarting from time %f (k=%d)', start_time, k));
    T = Ts(k,:);
    P = Ps(k,:);
end

function new_values = syspde(t, old_values)
    % Evaluate PDE system, in format to use in ode23s
    global N nb_vars bcs equations ppcs dim indices_interior boundary_mask
    global D1 D2 L
    
    new_values = Func(old_values',0,... % value of lambda is not used, so 0 will do
        N,nb_vars,bcs,equations,ppcs,dim,indices_interior,boundary_mask,...
        D1,D2,L);
end

function [indices_interior, boundary_mask] = getBoundaries(dim, N)
    % Helper function to get boundary mask and indices of interior domain
    if dim == 1
        boundary_indices = [1 N+1]';
    elseif dim == 2
        boundary_indices = [1:N+1]';
        for i = 1:N-1
            boundary_indices = cat(1, boundary_indices, [i*(N+1)+1 (i+1)*(N+1)]');
        end
        boundary_indices = cat(1, boundary_indices, [(N+1)^2-N:(N+1)^2]');
    else
        TODO
    end
    boundary_mask = zeros(1,(N+1)^dim);
    boundary_mask(boundary_indices) = 1;
    interior_mask = 1 - boundary_mask;
    indices_interior = find(interior_mask == 1);
end

function SuCCoMBe_S(data, bifurc_param_name)
    % Run steady state analysis
    
    global dim nb_vars
    
    num_params = data('numerical_params');
    bcs = data('bcs');
    eqs = data('equations');
    params = data('params');
    remove(params, bifurc_param_name);
    ppcs = data('ppcs'); % post processing commands
    
    equations = {}; % equations as strings ready to be evaluated
    param_names = params.keys;
    param_values = params.values;
    for i =1:nb_vars
        equations{i} = eqs(i);
        % Replace names of parameters (Ar, alpha,...)
        for j =1:length(params)
            equations{i} = strrep(equations{i}, param_names{j}, param_values{j});
        end
        % Replace bifurcation parameter
        equations{i} = strrep(equations{i}, bifurc_param_name, 'lambda');
    end
    for i = 1:length(ppcs)
        % Replace names of parameters including bifurcation parameter
        for j =1:length(params)
            ppcs{i} = strrep(ppcs{i}, param_names{j}, param_values{j});
        end
    end
    
    %% NUMERICAL PARAMETERS %%
    N = str2num(num_params('N')); % Number of collocation points in [-1,1], must be even.
    ds_target = str2num(num_params('ds')); % target arc-length step
    IG = str2num(num_params('LambdaInitialGuess')); % initial guess for lambda
    notification_freq = str2num(num_params('notificationFrequency')); % percentage (0-100) of total smax
    smax = str2num(num_params('sMax')); % maximum length
    restartfilename = num_params('restartFileNameSS'); % '' => start from beginning (no restart)
    savefilename = num_params('saveFileNameSS');
    i_param_plot = str2num(num_params('indexVarPlot')); % index of parameter to plot
    %% END OF PARAMETERS %%
    
    tic;
    warning ('off','all');
    fprintf('Starting steady state analysis on %s...\n', datestr(now));
    disp('(All warnings disabled)');
    if mod(N,2)==1
        fprintf('Error: N=%d is not even. End of code.\n', N);
        return
    end
    
    % Get index of points defining the boundaries of the domain
    [indices_interior, boundary_mask] = getBoundaries(dim, N);
    
    % Get Chebyshev points and weight factors
    [z,w] = getChebyshevPointsAndWeightFactors(N);
    % Differentiation matrix
    [D1,D2] = bcmatrix(w,z);
    L = 0; % just to be defined if dim==1
    if dim == 2
        I = eye(N+1);
        L = kron(I,D2) + kron(D2,I); % 2D Laplacian
    end
    % Initial solution [P0,T0]
    G = zeros(1,nb_vars*(N-1)^dim); % Not including boundary values (Dirichlet BC)
    for i = 1:nb_vars
        G((i-1)*(N-1)^dim+1:i*(N-1)^dim) = eval(bcs{i});
    end
    
    if restartfilename
        % Starting from previous save file
        fprintf('Restarting from "%s"\n', restartfilename);
        disp('(Assuming no parameter has changed except "smax")');
        load(restartfilename,'u','eigenv','bif_param','err','s','v','s_values');
        [nb_lines ~] = size(u{1});
        total_sol = zeros(nb_lines, nb_vars*(N-1)^dim+1);
        for i = 1:nb_vars
            total_sol(:,(i-1)*(N-1)^dim+1:i*(N-1)^dim) = u{i}(:,indices_interior);
        end
        total_sol(:,end) = bif_param;
        solz = total_sol(end-1,:);
        solf = total_sol(end,:);
        s_percent = min(100, s*100/smax);
        fprintf('%3.1f %% completed (at restart).\n', s_percent);
        next_notification = ceil((s_percent+0.1)/notification_freq)*notification_freq;
    else
        % Starting from the beginning
        % Initialisation
        s = 0;
        v = 1; % iteration index
        s_values = [];
        total_sol = []; % solutions, will be size nb_lines x (nb_vars*(N+1)+1)
        eigenv = []; % eigenvalues, will be size nb_lines x (nb_vars*(N-1))
        err = []; % errors, will be size nb_lines x 1

        %% INITIAL APPROXIMATION %%
        lambda0 = IG;
        f = @(x)Func(x,lambda0,N,nb_vars,bcs,equations,ppcs,dim,...
            indices_interior,boundary_mask,D1,D2,L);
        options = optimset('Display','off','Algorithm','levenberg-marquardt','TolFun',1e-6,'MaxIter',1000);
        [solz,F,exitflag,~,jacobian] = fsolve(f,G,options);
        if (exitflag <= 0)
            fprintf('Error at 1st approximation, exitflag = %3.5f. End of code.\n', exitflag);
            return
        end
        J = schur(jacobian);
        d = ordeig(J);
        if dim == 1
            E = jacobian\(-F);
            er = solz.\E';
            err = [err; max(abs(er))];
        end
        eigenv = [eigenv;d'];
        solz = [solz lambda0];
        total_sol = [total_sol;solz];

        %% SECOND APPROXIMATION %%
        lambda1 = 2*IG;
        f = @(x)Func(x,lambda1,N,nb_vars,bcs,equations,ppcs,dim,...
            indices_interior,boundary_mask,D1,D2,L);
        [solf,F,exitflag,~,jacobian] = fsolve(f,solz(1:nb_vars*(N-1)^dim),options);
        if (exitflag <= 0)
            fprintf('Error at 2nd approximation, exitflag = %3.5f. End of code.\n', exitflag);
            return
        end
        J = schur(jacobian);
        v = v+1;
        d = ordeig(J); 
        if dim == 1
            E = jacobian\(-F);
            er = solf.\E';
            err = [err; max(abs(er))];
        end
        eigenv = [eigenv;d'];
        solf = [solf lambda1];
        total_sol = [total_sol;solf];
        
        next_notification = notification_freq;
    end
    
    %% CODE INITIATION %%
    while s < smax
        v = v+1;
        test = 0; % number of recalculations in case of errors
        ds = ds_target;
        s1 = s;
        s_values = [s_values s];
        options = optimset('Display','off','Algorithm','levenberg-marquardt','TolFun',1e-4);

        f = @(x)Funct(x,ds,solz,solf,N,nb_vars,dim,bcs,equations,ppcs,...
            indices_interior,boundary_mask,D1,D2,L);
        [sol1,F,exitflag,~,jacobian] = fsolve(f,2*solf-solz,options);

        % Recalculation in case of error and end of the while loop %
        if (exitflag <= 0)
            fprintf('Error at length %3.5f, exitflag = %3.5f. Reduce step.\n', s, exitflag);    
            while (test < 5) && (exitflag <= 0)
                test = test + 1;
                f = @(x)Funct(x,ds_target/(10^test),solz,solf,N,nb_vars,...
                    dim,bcs,equations,ppcs,indices_interior,boundary_mask,D1,D2,L);
                [sol1,F,exitflag,~,jacobian]  = fsolve(f,2*solf-solz,options);
            end
            ds = ds_target/(10^test);
        end

        if (exitflag <= 0)
            fprintf('Error at length %3.4f , end of code .\n', s);
            return
        end

        %% SAVING SOLUTION %%
        total_sol = [total_sol;sol1];
        % Jacobian and eigenvalues
        J = jacobian(1:nb_vars*(N-1)^dim,1:nb_vars*(N-1)^dim);
        J = schur(J);
        d = ordeig(J);
        eigenv = [eigenv;d'];
        if dim == 1
            E = jacobian\(-F);
            er = sol1.\E';
            err = [err; max(abs(er))];
        end
        
        % reassigning the 2 approximations
        solz = solf;
        solf = sol1;

        s_percent = min(100, (s1 + ds)*100/smax);
        if s_percent >= next_notification
            fprintf('%3.1f %% completed.\n', s_percent);
            next_notification = next_notification + notification_freq;
        end

        s = s+ds;

    end

    %% SAVE TO FILE
    nb_lines = v;
    u = {}; % Individual solutions
    for i = 1:nb_vars
        tmp = eval(bcs{i}).*ones(nb_lines,(N+1)^dim).*repmat(boundary_mask,nb_lines,1);
        tmp(:,indices_interior) = total_sol(:,(i-1)*(N-1)^dim+1:i*(N-1)^dim);
        u{i} = tmp;
    end
    bif_param = total_sol(:,nb_vars*(N-1)^dim+1);

    save(savefilename,'u','eigenv','bif_param','z','err','s','v','s_values','N');
    fprintf('File "%s" saved\n',savefilename);

    %% FIND CRITICAL POINTS AND PLOT %%
    [saddlepoints, hopfs] = plotResult(bif_param, u{i_param_plot}, N, ...
        eigenv(:,(i_param_plot-1)*(N-1)^dim+1:i_param_plot*(N-1)^dim), ...
        s_values, dim, 0);
    printPointsCoords(saddlepoints, hopfs);
    fprintf('Time elapsed: %.2f seconds\n',toc)
    
    return
end

function [saddlepoints, hopfs] = plotResult(Gr, T, N, eigenv, s_values, dim, verbose)
    % Plot "S-curve" with colours indicating stability type
    % Returns list of Hopf and saddle points
    % First, determine types to plot with colours
    k_centre = getSpatialCentreIndex(dim, N);
    saddlepoints = [];
    hopfs = [];
    if 0
        % Plot only "S-curve"
        figure('name', 'SuCCoMBe (steady state)');plot(Gr,T(:,k_centre));
    else
        % Plot with colours
        figure('name', 'SuCCoMBe (steady state)');
        clf; hold on;xlabel('Gr');ylabel('T_{centre}');
        title('Steady state')
        [nb_lines ~] = size(T);
        v = 3;
        types = [-1;-1];
        while v < nb_lines+1
            if any(imag(eigenv(v,:)))
                real_vals = sort(real(eigenv(v,:)));
                if real_vals(end) < 0
                    % unstable limit cycle
                    printVerbose(sprintf('s=%f\tUnstable limit cycle\n',s_values(v-2)), verbose);
                    types = [types; 1];
                elseif real_vals(end) > 0
                    % stable limit cycle
                    printVerbose(sprintf('s=%f\tStable limit cycle\n',s_values(v-2)), verbose);
                    types = [types; 2];
                else
                    % Hopf point
                    printVerbose(sprintf('s=%f\tHopf bifurcation\n',s_values(v-2)), verbose);
                    types = [types; 3];
                    hopfs = [hopfs; Gr(v) T(v,k_centre)];
                end
            else
                real_vals = sort(eigenv(v,:));
                if real_vals(end) < 0
                    % stable
                    printVerbose(sprintf('s=%f\tStable\n',s_values(v-2)), verbose);
                    types = [types; 4];
                elseif real_vals(end) > 0
                    % unstable
                    printVerbose(sprintf('s=%f\tUnstable\n',s_values(v-2)), verbose);
                    types = [types; 5];
                else
                    % Saddle point
                    printVerbose(sprintf('s=%f\tSaddle point\n',s_values(v-2)), verbose);
                    types = [types; 6];
                    saddlepoints = [saddlepoints; Gr(v) T(v,k_centre)];
                end
            end
            v = v + 1;
        end

        % Plot with colours
        colour = {'red','red','red','black','black','black'};
        linestyle = {'-','--','x','-','--','o'};
        extra_x_start = []; % no previous segment to connect from
        extra_y_start = []; % no previous segment to connect from
        i_start = 3;
        type = types(i_start);
        i_end = i_start + 1;
        while i_end < nb_lines+1
            if types(i_end) == types(i_start)
                i_end = i_end + 1;
                continue
            else
                % Identify Hopf and saddle points
                if isequal(sort([types(i_end) types(i_start)]),[4  5])
                    % Saddle point somewhere in between
                    plot_separator = 1;
                    x_separator = (Gr(i_end-1)+Gr(i_end))/2;
                    y_separator = (T(i_end-1,k_centre)+T(i_end,k_centre))/2;
                    plot(x_separator, y_separator, ...
                        'Color','black', 'Marker','o', ...
                        'MarkerFaceColor','black', 'MarkerSize',3);
                    saddlepoints = [saddlepoints; x_separator y_separator];
                elseif any(types(i_start)==[4 5 6]) && any(types(i_end)==[1 2 3])
                    % Hopf somewhere in between
                    plot_separator = 1;
                    x_separator = (Gr(i_end-1)+Gr(i_end))/2;
                    y_separator = (T(i_end-1,k_centre)+T(i_end,k_centre))/2;
                    plot(x_separator, y_separator, ...
                        'Color','red', 'Marker','x', 'MarkerSize',8);
                    hopfs = [hopfs; x_separator y_separator];
                else
                    plot_separator = 0;
                end
                % Plot segment
                printVerbose(sprintf('Plotting s from %d to %d with style %d\n',...
                    i_start,i_end, types(i_start)), verbose);
                if plot_separator
                    x_array = [Gr(i_start:i_end-1);x_separator];
                    y_array = [T(i_start:i_end-1,k_centre);y_separator];
                else
                    x_array = Gr(i_start:i_end);
                    y_array = T(i_start:i_end,k_centre);
                end
                x_array = [extra_x_start;x_array];
                y_array = [extra_y_start;y_array];
                if plot_separator
                    extra_x_start = x_separator;
                    extra_y_start = y_separator;
                else
                    extra_x_start = [];
                    extra_y_start = [];
                end
                plot(x_array,y_array,'Color',char(colour(types(i_start))),...
                    'LineStyle',char(linestyle(types(i_start))));
                % Start next segment
                i_start = i_end;
                i_end = i_end + 1;
                continue
            end
        end
        if i_start < nb_lines
            printVerbose(sprintf('Plotting s from %d to %d with style %d\n',...
                i_start,nb_lines, types(i_start)), verbose);
            x_array = [extra_x_start;Gr(i_start:nb_lines)];
            y_array = [extra_y_start;T(i_start:nb_lines,k_centre)];
            plot(x_array,y_array,'Color',char(colour(types(i_start))),...
                 'LineStyle',char(linestyle(types(i_start))));
        end
    end
    return
end

function printPointsCoords(saddlepoints, hopfs)
    % Print message to show coordinates of Hopf and saddle points
    % Saddle points
    [nb_sp ~] = size(saddlepoints);
    tmp1 = '';
    tmp2 = '';
    if nb_sp == 0
        tmp1 = 'No';
    elseif nb_sp == 1
        tmp1 = '1';
        tmp2 = sprintf(':%s',char(10));
        tmp2 = [tmp2 sprintf('\tGr=%s\tT=%s',saddlepoints(1,1),saddlepoints(1,2))];
    else
        tmp1 = sprintf('%d',nb_sp);
        tmp2 = sprintf('s:%s',char(10));
        for i = 1:nb_sp
            tmp2 = [tmp2 sprintf('\tGr=%s\tT=%s',saddlepoints(i,1),saddlepoints(i,2))];
            if i < nb_sp
                tmp2 = [tmp2 char(10)];
            end
        end
    end
    fprintf('%s saddle point%s\n',tmp1,tmp2);
    % Hopf
    [nb_hp ~] = size(hopfs);
    tmp1 = '';
    tmp2 = '';
    if nb_hp == 0
        tmp1 = 'No';
    elseif nb_hp == 1
        tmp1 = '1';
        tmp2 = sprintf(':%s',char(10));
        tmp2 = [tmp2 sprintf('\tGr=%s\tT=%s',hopfs(1,1),hopfs(1,2))];
    else
        tmp1 = sprintf('%d',nb_hp);
        tmp2 = sprintf(':%s',char(10));
        for i = 1:nb_hp
            tmp2 = [tmp2 sprintf('\tGr=%s\tT=%s%s',hopfs(i,1),hopfs(i,2),char(10))];
        end
    end
    fprintf('%s Hopf point%s\n',tmp1,tmp2);
end

function [x,w] = getChebyshevPointsAndWeightFactors(N)
    x = cos((0:N)*pi/N)'; % Chebyshev collocation points
    % Weight factors (which will stay constant all the time)
    w = (-1).^(0:N)';
    w(1) = w(1)/2.;
    w(end)= w(end)/2;
end

function [D1,D2] = bcmatrix(w,x)
    % BCMATRIX constructs the first and second order differentiation matrices, D1
    % and D2, corresponding to the barycentric weights w and grid points x. Note
    % that w and x must be column vectors of the same size.
    N = length(x)-1; ii = (1:N+2:(N+1)^2)';
    Dw = repmat(w',N+1,1) ./ repmat(w,1,N+1) - eye(N+1);
    Dx = repmat(x ,1,N+1) - repmat(x',N+1,1) + eye(N+1);
    D1 = Dw ./ Dx;
    D1(ii) = 0; D1(ii) = - sum(D1,2);
    D2 = 2*D1 .* (repmat(D1(ii),1,N+1) - 1./Dx);
    D2(ii) = 0; D2(ii) = - sum(D2,2);
end

function [] = printVerbose(my_string, verbose)
    % print message only if verbose is on
    if verbose
        fprintf(my_string)
    end
end

function y = Funct(x,ds,solz,solf,N,nb_vars,dim,bcs,equations,ppcs,...
    indices_interior,boundary_mask,D1,D2,L)
    % Arclength continuation function
    y = Func(x,x(nb_vars*(N-1)^dim+1),...
        N,nb_vars,bcs,equations,ppcs,dim,indices_interior,boundary_mask,...
        D1,D2,L);
    y(nb_vars*(N-1)^dim+1) = (sum((x - solf).*conj(solf-solz))/ds)-ds;
    return
end

function y = Func(x,lambda,...
    N,nb_vars,bcs,equations,ppcs,dim,indices_interior,boundary_mask,...
    D1,D2,L)
    % Function evaluation of PDE system to return time derivative of
    % unknowns.
    
    % Break vector of variables back in user variables
    u = {};
    for i = 1:nb_vars
        tmp = eval(bcs{i}).*ones(1,(N+1)^dim).*boundary_mask;
        tmp(:,indices_interior) = x(:,(i-1)*(N-1)^dim+1:i*(N-1)^dim);
        u{i} = tmp';
    end
    % Apply preprocessing commands
    for i = 1:length(ppcs)
        eval(char(ppcs{i}));
    end
    % Apply equations 
    y = [];
    for i = 1:nb_vars
        tmp = eval(char(equations{i}));
        y = [y; tmp(indices_interior)];
    end
    return
end
