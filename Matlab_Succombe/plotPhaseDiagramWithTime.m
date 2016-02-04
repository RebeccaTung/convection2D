% To create figs 2a, 2b, 3b and a movie for supplementary material

function plotPhaseDiagramWithTime()
    mat_file = 'exports/export_Napoli_smooth.mat';
    % Load variables saved in export file
    load(mat_file, 'times','dim','u');
    
    % Load manually other values
    N = 8;
    Ti = 0;
    T0 = 0;
    Pi = 0;
    P0 = 0;
    Le = 4;
    n = 3;
    
    Gr = 1.3e-7;
    Tc = 360;
    alpha = 0.5;
    delta = 1e-3; 
    mu = 1e-4;
    Ar = 40;
    Kc = 1e10;
    deltaE = 80e3;
    phi0 = 0.03;
    
    R = 8.3144621; % J.mol^-1.K^-1
    M_cao = 0.15; % kg/mol
    Mcaco3 = 0.35; % kg/mol
    Mco2 = 0.018; % kg/mol
    rho_cao = 3.3e3; % kg/m3
    rho_caco3 = 2e3; % kg/m3
    rho_co2 = 0.9e3; % kg/m3
   
    %%%%%% New params needed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Kappa_m = 1e-6; % Thermal diffusion coefficient of the mixture
    half_width = 0.05;%0.75; % half width (m) of the shear zone
    background_vel = 0.01 / 3.14e7; % 1cm/year
    gammadotb = background_vel/(2*half_width);
    gamma_dot_zero = 1e8; % pre-exponential factor in flow law
    %slope = 8.2e10; % smooth=6.17e9, sharp=8.2e10
    slope = 6.17e9; % smooth=6.17e9, sharp=8.2e10
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	Ps = u{1};
    Ts = u{2};
	eps_dot = (1-Ps).^n.*exp(-(1-alpha)*Ar./(1+delta*Ts));
    i_c = N/2+1;%N+1;%N/2+1; % index of fault centre point
    P_centre = Ps(:,i_c);
    T_centre = Ts(:,i_c);
    eps_dot_centre = eps_dot(:,i_c);
    
    day = 1/365.25;
    delta_t_yr = 4*day; % time step for movie (in years)
    
    year = 3.15569e7;
    times_yr = (half_width^2)/Kappa_m*times/year;
    x = cos([0:N]*pi/N)';
    export_dir = 'pics_DRT';
    if exist(export_dir) ~= 7
        mkdir(export_dir)
        disp(sprintf('\nCreated export directory "%s"', export_dir));
    end

    % Chemistry. See paper part I, appendix equation A.4

    w_rel = (rho_caco3/rho_cao)*(M_cao/Mcaco3)*Kc*exp(-deltaE./(R*Tc*(1+delta.*Ts)));
    s = w_rel./(1+w_rel);
    delta_phi_chem = (1-phi0)./(1+(rho_co2/rho_cao)*(M_cao/Mco2).*(1./s));
    phi = phi0 + delta_phi_chem;

    % Chebychev points and weight factors (for plotting on fine interpolation)
    x = cos([0:N]*pi/N)'; % Chebyshev collocation points
    w = (-1).^[0:N]';
    w(1) = w(1)/2.;
    w(end)= w(end)/2;

    ref_eps_dot = eps_dot(1,1); % to normalise strain rate wrt boundary rate
    
    % Plot displacement with time
    N_times = length(times);
    start_k = 1; % Starting time index from which we track displacement
    end_k = N_times; % Ending time index
    i_middle = N+1; % spatial index locating the centre of the shear zone
    %V = zeros(N_times, N+1);
    V = 2*half_width*gamma_dot_zero*ones(N_times, N+1);
    for i=2:N+1
        %V(:,i) = V(:,i-1) + eps_dot(:,i-1)*(x(i-1)-x(i));
        V(:,i) = V(:,i-1) + 2*half_width*gamma_dot_zero*...
            (max(0,eps_dot(:,i-1))+max(0,eps_dot(:,i)))*(x(i-1)-x(i))/2;
    end
    V = V - V(:,i_middle)*ones(1,N+1); % impose V=0 at centre of shear zone

    % Integrate Velocity V in time to get displacement
    u = zeros(N_times, N+1);
    for k=start_k+1:end_k
        %u(k,:) = u(k-1,:) + V(k-1,:)*(times(k)-times(k-1));
        u(k,:) = u(k-1,:) + half_width^2/Kappa_m*((V(k-1,:)+V(k,:))*(times(k)-times(k-1))/2);
    end
    format long;
    year = 3.15569e7;
    mm = 0.001;
    times_yr = (half_width^2)/Kappa_m*times/year;
    disp_mm = 2*u(:,1)/mm;
    
    values = slope*times_yr+disp_mm; min_val=min(values); max_val=max(values);
    fig5 = figure(55);
    set(fig5, 'Position',[500,360,500,300])%[left bottom width height])
    ax0 = subplot(311);% dummy axes for the box and background color
    set (ax0, 'Box', 'on', 'Color', 'white', 'XTick', [], 'YTick', []);
    ax1 = axes ('Position', get (ax0, 'Position'));
    set (ax1, 'Box','off','Color','none','YAxisLocation','left','XTickLabel',[],'YTickLabel',[]);
    hold(ax1, 'on');
    plot(times_yr, values*1e7/gamma_dot_zero, 'color', [0 0 0], 'LineStyle', '-', 'LineWidth',2);
    xlabel('Time');% (yr)
    ylabel('Displacement');% (mm)
    grid on;

    ax0 = subplot(312);% dummy axes for the box and background color
    set (ax0, 'Box', 'on', 'Color', 'white', 'XTick', [], 'YTick', []);
    ax1 = axes ('Position', get (ax0, 'Position'));
    set (ax1, 'Box','off','Color','none','YAxisLocation','left','XTickLabel',[],'YTickLabel',[]);
    hold(ax1, 'on');
    plot(times_yr, phi(:,i_c), 'color', [0 0 0], 'LineStyle', '-', 'LineWidth',2);
    xlabel('Time');% (yr)
    ylabel({'CO2 release';'(centre)'});% (mm) 
    grid on;

    ax0 = subplot(313);% dummy axes for the box and background color
    set (ax0, 'Box', 'on', 'Color', 'white', 'XTick', [], 'YTick', []);
    ax1 = axes ('Position', get (ax0, 'Position'));
    set (ax1, 'Box','off','Color','none','YAxisLocation','left','XTickLabel',[],'YTickLabel',[]);
    hold(ax1, 'on');
    cumulated_phi = sum(phi,2); % cumulated over whole fault (in space)
    plot(times_yr, cumulated_phi, 'color', [0 0 0], 'LineStyle', '-', 'LineWidth',2);
    xlabel('Time');% (yr)
    ylabel({'CO2 release';'(cumulated)'});% (mm) 
    grid on;
return
    
    
    
    
    
    
    
    
    

    % Plot phase diagram (figure 2a of the paper)
    fig1 = figure(1);
    set(fig1, 'NumberTitle','off','Name', 'fig. 2a for the paper');
    set(fig1, 'Position',[30,0,800,300])
    ax0=subplot(211);% dummy axes for the box and background color
    set(ax0,'position',[0.13 0.2 0.78 0.66])%[left bottom width height])
    set (ax0, 'Box', 'on', 'Color', 'white', 'XTick', [], 'YTick', []);
    % first axes for left y-axis (phi, s)
    ax1 = axes ('Position', get (ax0, 'Position'));
    set (ax1, 'Box', 'off', 'Color', 'none', 'YAxisLocation', 'left');
    % second axes for right y-axis assuming common x-axis controlled by ax1
    ax2 = axes ('Position', get (ax0, 'Position'));
    set (ax2, 'Box', 'off', 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right');
    hold(ax1, 'on');
    hold(ax2, 'on');
    %xlim(ax1, [1,7.5])
    ylim(ax1, [0,1])
    set(ax1,'Ytick',0:0.5:1)
    set(ax2,'Ytick',0:1e-7:2e-7)
    %xlim(ax2, [1,7.5]);
    ylim(ax2, [0,4e-7])
    %set(ax1,'Xtick',0:10:49)
    %set(ax1,'XtickLabel',0:1:10)
    %set(ax2,'Ytick',0:5.e-5:2.e-4)
    %set(ax2,'YtickLabel',0:5.e-8:1.e-7)
    max_eps_dot_centre = max(eps_dot_centre);
    set(ax1, 'FontUnits', 'points', 'FontSize', 14);
    set(ax2, 'FontUnits', 'points', 'FontSize', 14);
    p3 = plot(ax1, times_yr, 1*s(:,i_c), 'color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth',2);
    p2 = plot(ax1, times_yr, 1*phi(:,i_c), 'color', [0.2 0.2 0.2], 'LineStyle', ':', 'LineWidth',3);
    p1 = plot(ax2, times_yr, max(0,eps_dot_centre), 'k', 'LineWidth',3);
    xlabel(ax1, 'Time (yr)');
    ylabel(ax1, {'Porosity','oxide content'});
    ylabel(ax2, '\textsf{Strain rate ($\dot{\epsilon}/\dot{\epsilon_0}$)}','interpreter','latex');
    leg2a = legend([p1,p2,p3], 'strain rate', 'porosity', 'oxide content');
    set(leg2a, 'Location' ,'NorthWest');
    leg2atext = findobj(leg2a,'type','text');
    set(leg2atext,'FontSize',10)
    save_filename = sprintf('%s/fig2a.png',export_dir);
    saveas(gcf, save_filename, 'png');
    
    % Plot phase diagram (figure 2b of the paper)
    fig2 = figure(2);
    s1=subplot(211);
    set(s1,'position',[0.13 0.22 0.78 0.70])%[left bottom width height])
    set(fig2, 'NumberTitle','off','Name', 'fig. 2b for the paper');
    set(fig2, 'Position',[30,400,800,300])
    start_index = length(T_centre)/3;
    plot((1+delta*T_centre(start_index:end))*Tc - 273.15, P_centre(start_index:end)*1.81/2.71+0.9/2.71, 'k', 'LineWidth',2);
    set(gca, 'FontUnits', 'points', 'FontSize', 16);
    xlabel('T_{core} (^{\circ}{\rm C})'); 
    ylabel('P^f_{core} / P^{rock}_{boundary}');
    %ylabel('P^f_{core} / \rho g h');
    %xlim([0.9;2.45]);
    ylim([0.3;1.7]);
    %set(gca,'Xtick',1:0.5:3)
    %set(gca,'XtickLabel',1:0.1:2)
    %set(gca,'Ytick',0:0.5:3)
    %set(gca,'YtickLabel',0:0.5:3)
    save_filename = sprintf('%s/fig2b.png',export_dir);
    saveas(gcf, save_filename, 'png');

    % Plot Batman diagram (figure 3b of the paper)
    fig3 = figure(3);
    s3=subplot(211);
    x_fine =[-1:0.001:1]'; % points to plot fine curves in between collocation points
    % Find index of max strain rate (where centre is still the highest, before Batman!)
    shift_start_index = floor(length(eps_dot_centre)/4);
    [max_eps_dot_centre, k_peak] = max(eps_dot_centre(shift_start_index:end));% index of strain rate peak
    k_peak = k_peak + shift_start_index - 1;
    eps_dot_fine = bcinterp(w, x, eps_dot(k_peak,:), x_fine)';
    while eps_dot_fine(1) < eps_dot_fine(2)
        k_peak = k_peak-7;
        eps_dot_fine = bcinterp(w, x, eps_dot(k_peak,:), x_fine)';
        disp(sprintf('  Moving k_peak to %d', k_peak));
    end
    %set(s3,'position',[0 0 1 1])%[left bottom width height])
    set(fig3, 'NumberTitle','off','Name', 'fig. 3b for the paper');
    ax10=subplot(211);
    set(ax10,'position',[0.17 0.08 0.72 0.85])%[left bottom width height])
    set (ax10, 'Box', 'on', 'Color', 'white', 'XTick', [], 'YTick', []);
    % second axes for right y-axis assuming common x-axis controlled by ax1
    ax12 = axes ('Position', get (ax10, 'Position'));
    set (ax12, 'Box', 'off', 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right');
    % first axes for left y-axis (phi, s)
    ax11 = axes ('Position', get (ax10, 'Position'));
    set (ax11, 'Box', 'off', 'Color', 'none', 'YAxisLocation', 'left');
    hold(ax11, 'on');
    hold(ax12, 'on');
    phi_peak_fine = bcinterp(w, x, phi(k_peak,:), x_fine)';
    eps_dot_fine = bcinterp(w, x, eps_dot(k_peak,:), x_fine)';
    full_x = [-x-1;x(end:-1:1)+1];
    full_x_fine = [-x_fine-1;x_fine(end:-1:1)+1];
    full_phi_peak = [phi(k_peak,:),phi(k_peak,end:-1:1)];
    full_phi_peak_fine = [phi_peak_fine,phi_peak_fine(end:-1:1)];
    full_s_peak = [s(k_peak,:),s(k_peak,end:-1:1)];
    full_eps_dot_peak = [eps_dot(k_peak,:),eps_dot(k_peak,end:-1:1)];
    full_eps_dot_peak_fine = [eps_dot_fine,eps_dot_fine(end:-1:1)];
    p3b1 = plot(ax12, full_x, full_eps_dot_peak, 'r', 'LineWidth',2);
    %p3b1b = plot(ax12, full_x_fine, full_eps_dot_peak_fine, 'k', 'LineWidth',1);
    p3b2 = plot(ax11, full_x, full_phi_peak, 'b', 'LineWidth',2);
    p3b3 = plot(ax11, full_x, full_s_peak, 'g', 'LineWidth',2);
    set(ax11, 'FontUnits', 'points', 'FontSize', 14);
    set(ax12, 'FontUnits', 'points', 'FontSize', 14);
    xlabel(ax11, 'Distance across Process Zone');
    ylabel(ax11, {'Porosity','oxide content'});
    ylabel(ax12, '\textsf{Strain rate $\dot{\gamma}$ (1/s)}','interpreter','latex');
    ylim(ax11, [0,1]);
    %ylim(ax2, [0,4e-7]);
    set(ax11,'Xtick',-2:2:2)
    set(ax11,'XtickLabel',[])
    %set(ax11,'Ytick',0:0.1:1)
    %set(ax11,'YtickLabel',0:0.1:1)
    %set(ax12,'Ytick',0:5e-8:2e-7)
    leg3b = legend([p3b1,p3b2,p3b3], 'strain rate', 'porosity', 'oxide content');
    set(leg3b, 'Location' ,'NorthWest');
    leg3btext = findobj(leg3b,'type','text');
    set(leg3btext,'FontSize',12)
    save_filename = sprintf('%s/fig3b.png',export_dir);
    saveas(gcf, save_filename, 'png');

    % Figure 3 for Geology GSA paper
    fig3gsa = figure(31);
    s3=subplot(211);
    %set(s3,'position',[0 0 1 1])%[left bottom width height])
    set(fig3gsa, 'NumberTitle','off','Name', 'fig. 3 for Geology GSA paper');
    ax30 = subplot(211);
    set(ax30,'position',[0.17 0.12 0.72 0.76])%[left bottom width height])
    set (ax30, 'Box', 'on', 'Color', 'white', 'XTick', [], 'YTick', []);
    %hold(ax30, 'on');
    % second axes for right y-axis assuming common x-axis controlled by ax1
    ax32 = axes ('Position', get (ax30, 'Position'));
    set (ax32, 'Box', 'off', 'Color', 'none', 'YAxisLocation', 'right');
    set(ax32, 'FontUnits', 'points', 'FontSize', 14, 'XAxisLocation', 'top');
    % first axes for left y-axis (phi, s)
    ax31 = axes ('Position', get (ax30, 'Position'));
    set (ax31, 'Box', 'off', 'Color', 'none', 'YAxisLocation', 'left');
    hold(ax31, 'on');
    hold(ax32, 'on');
    p3gsa1 = plot(ax32, full_eps_dot_peak, full_x, 'k', 'LineWidth',2);
    p3gsa2 = plot(ax31, 100*full_phi_peak, full_x, 'k--', 'LineWidth',2);
    p3gsa3 = plot(ax31, 100*full_s_peak, full_x, 'Color', 0.6*[1 1 1], 'LineWidth',2);
    set(ax31, 'FontUnits', 'points', 'FontSize', 14);
    xlabel(ax31, 'Porosity, oxide content (%)');
    ylabel(ax31, '(Half) Fault Core');
    xlabel(ax32, '\textsf{Strain rate ($\dot{\epsilon}/\dot{\epsilon_0}$)}','interpreter','latex');
    ylim(ax31, [-2,0.15]);
    ylim(ax32, [-2,0.15]);
%     set(ax32,'Xtick',0:0.5e-7:2e-7)
%     set(ax32,'XtickLabel',0:0.5e-7:2e-7)
    set(ax31,'Ytick',0)
    set(ax31,'YtickLabel',[])
    set(ax32,'Ytick',0)
    set(ax32,'YtickLabel',[])
    set(ax32,'Xtick',[0,1e-7,2e-7])
    set(ax32,'XtickLabel',[0,1e-7,2e-7])
    leg3gsa = legend([p3gsa1,p3gsa2,p3gsa3], 'strain rate', 'porosity', 'oxide content');
    set(leg3gsa, 'Location' ,'SouthEast');
    leg3gsatext = findobj(leg3gsa,'type','text');
    set(leg3gsatext,'FontSize',12)
    save_filename = sprintf('%s/fig3_GSA.png',export_dir);
    saveas(gcf, save_filename, 'png');
    disp(sprintf('Saved figure %s', save_filename));

    % Plot time evolution (supplementary material for paper)
    fig3 = figure(4);clf;
    set(fig3, 'NumberTitle','off','Name', 'To create movie for supplementary material');
    markersize =5;
    i_pic = 1;
    prev_pic_time_yr = -1e99;
    start_index = floor(length(times)/3); % Thomas cutting the initial loop (well, 1/3 of the total)
    for k=start_index:1:length(times)
      if mod(k-1,40)== 0 || k==length(times)
      
      
      %if k>9624 && k<13000 && mod(k-1,10)== 0
      %if times_yr(k)>prev_pic_time_yr + delta_t_yr
        prev_pic_time_yr = times_yr(k);
        clf;
        subplot(3,1,1); hold on
        % title
        nb_yr = floor(times_yr(k));
        nb_days = (times_yr(k) - nb_yr)*365.25;
        title(sprintf('Time=%d years, %.2f days',nb_yr,nb_days));
        %title(sprintf('Time=%d years, %.2f days (time index=%d)',nb_yr,nb_days, k));
        plot(times_yr, log10(eps_dot_centre/ref_eps_dot)); xlabel('Time (yr)'); 
        ylabel('$log_{10}(\dot{\gamma}/\dot{\gamma_0})$','interpreter','latex');
        ylim([-8,5]);
        plot(times_yr(k), log10(eps_dot_centre(k)/ref_eps_dot), 'ro','MarkerSize',markersize,'MarkerFaceColor','red') ;
        subplot(3,1,2); hold on
        plot(times_yr, s(:,i_c), 'g'); 
        plot(times_yr, phi(:,i_c)); 
        xlabel('Time (yr)'); ylabel('Porosity/s');
        plot(times_yr(k), s(k,i_c), 'ro','MarkerSize',markersize,'MarkerFaceColor','red');
        subplot(3,1,3); hold on
        plot((1+delta*T_centre(start_index:end))*Tc - 273.15, P_centre(start_index:end)*1.81/2.71+0.9/2.71);
        xlabel('T_{core} (C)'); ylabel('P^f_{core} / \rho g h');
        plot((1+delta*T_centre(k))*Tc - 273.15, P_centre(k)*1.81/2.71+0.9/2.71, 'ro','MarkerSize',markersize,'MarkerFaceColor','red');
        
        save_filename = sprintf('%s/frame%05d.png',export_dir,i_pic);
        saveas(gcf, save_filename, 'png');
        i_pic = i_pic + 1;
        return
        %pause
      end
    end  
end

function ff = bcinterp(w,x,f,xx)
    % BCINTERP evaluates the barycentric interpolation formula to find the values
    % ff of a rational interpolant with barycentric weights w, grid points x, and
    % grid data f, at the points xx. Note that w, x, and f must be column vectors
    % of the same size.
    [mask,index] = ismember(xx,x);
    invmask = (mask==0);
    xxx = xx(invmask);
    ff = zeros(length(xx),1);
    numer = zeros(length(xxx),1);
    denom = zeros(length(xxx),1);
    for i=1:length(x)
        temp = w(i)./(xxx-x(i));
        numer = numer + temp*f(i);
        denom = denom + temp;
    end
    ff(invmask) = numer./denom;
    ff(mask) = f(index(mask));
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