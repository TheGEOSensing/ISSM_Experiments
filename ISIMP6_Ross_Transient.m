

% Load the model we have derived from the inversion
cd '/Users/rishi/Desktop/ISSM_Project_Ross/Ross_Ice_Shelf/'
md = loadmodel('/Users/rishi/Desktop/ISSM_Project_Ross/Ross_Ice_Shelf/Ross_Grounded.mat');


%% USE THE ISIMP6 forcing for the transient simulation
% ------------------------------------------------------------------
%     - model_name (string): name of the climate model and scenario
%     - suppported models:
%             2.6 scenario             8.5 scenario
%             ---------------------------------------------
%             cnrm-cm6-1_ssp126        cnrm-cm6-1_ssp585
%             ipsl-cm5a-mr_rcp2.6      ipsl-cm5a-mr_rcp8.5
%             noresm1-m_rcp2.6         noresm1-m_rcp8.5
% ------------------------------------------------------------------
md.smb = interpISMIP6AntarcticaSMB (md,'noresm1-m_rcp8.5');
md.basalforcings = interpISMIP6AntarcticaOcn (md,'noresm1-m_rcp8.5');


%% Indicate the components of transient to activate
md.transient.ismasstransport = 1;
md.transient.isstressbalance = 1;
md.transient.isgroundingline = 1;
md.transient.ismovingfront = 0;
md.transient.isthermal = 0;

%Disable inverse method
md.inversion.iscontrol=0;

%Initialize fields for transient and add boundary conditions
md.initialization.vx = md.results.StressbalanceSolution.Vx;
md.initialization.vy = md.results.StressbalanceSolution.Vy;
md.initialization.vel = md.results.StressbalanceSolution.Vel;
md.masstransport.spcthickness = NaN*ones(md.mesh.numberofvertices,1);

% set the grounded ice melting rate
md.basalforcings.groundedice_melting_rate = ones(md.mesh.numberofvertices,1);
md.basalforcings.geothermalflux = zeros(md.mesh.numberofvertices,1);
md.basalforcings.islocal = 1;


%% Request additional outputs
md.transient.requested_outputs={'default','IceVolume',...
    'IceVolumeAboveFloatation','TotalSmb','TotalFloatingBmb',...
    'AverageButtressing'};

%Specify time steps and length of simulation (years)
md.timestepping.start_time = 2014;
md.timestepping.time_step = 0.5;
md.timestepping.final_time = 2050;

% Solve transient solution and change messages provided
md.cluster = generic('name',oshostname,'np',4);
md.verbose = verbose('solution',false);
md=solve(md,'Transient');


%% Plot results
plotmodel(md, 'data', md.results.TransientSolution(1).Vel,'title#1', 'Velocity 2014 (m/yr)',...
  'data', md.results.TransientSolution(end).Vel,'title#2', 'Velocity 2050 years (m/yr)',...
  'caxis#1',([0 3500]),'caxis#2',([0 3500]));


%%

% buttressing and backstress for the initial time
for i = 1:size(md.results.TransientSolution,2)
    
    % load the velocity datasets
    vx = md.results.TransientSolution(i).Vx;
    vy = md.results.TransientSolution(i).Vy;
    [Kn_i, Backstress_Furst_i] = ...
        Furst_buttressing(md, vx, vy, false);

    Kn(:,i) = Kn_i; % Calculate Buttressing ratio for all time
    Backstress_Furst(:,i) = Backstress_Furst_i; % Backstress
end

%% Plot them and show the difference between them
plotmodel(md,'data',Backstress_Furst(:,1),'title', 'Backstress t=0 years','caxis#1',[0 3e5], ...
              'data',Backstress_Furst(:,end),'title', 'Backstress t=20 years','caxis#2',[0 3e5]),...
colormap(brewermap(50,'Spectral'))

%%

plotmodel(md,'data',(Backstress_Furst(:,1) - Backstress_Furst(:,end)),'title','Difference','caxis#1', [-1e5 1e5]);
colormap(brewermap(50,'RdBu'))

%%

% time-series of transient solution:
Time_period = [md.results.TransientSolution.time];

% the parameters for time series analysis
Backstress_TS = [md.results.TransientSolution.AverageButtressing];
VAF_TS = [md.results.TransientSolution.IceVolume];
Thickness_TS = nanmean([md.results.TransientSolution.Thickness],1);
Velocity_TS = nanmean([md.results.TransientSolution.Vel],1);

% plot the time series analysis
% Average backstress
subplot(4,1,1)
plot(Time_period,Backstress_TS, 'LineWidth', 1.5, 'color',"#80B3FF")
ylabel('GL Backstress (%)', 'FontSize', 12)
ax1 = gca;

% Thickness
subplot(4,1,2)
plot(Time_period,Thickness_TS, 'LineWidth', 1.5, 'color',"blue")
ylabel('Thickness (meters)', 'FontSize', 12)
ax2 = gca;

% Volume above floatation
subplot(4,1,3)
plot(Time_period,VAF_TS, 'LineWidth', 1.5, 'color',"blue")
ylabel('Volume Above Floatation (m^3)', 'FontSize', 12)
ax3 = gca;

subplot(4,1,4)
plot(Time_period,Velocity_TS, 'LineWidth', 1.5, 'color',"blue")
ylabel('Velocity (m/yr)', 'FontSize', 12)
ax4 = gca;


% Apply settings to all subplots
axs = [ax1, ax2, ax3, ax4]; % Store all axes handles

for i = 1:length(axs)
    axs(i).FontSize = 14;
    axs(i).FontWeight = 'bold';
    axs(i).LineWidth = 1.5; % Optional: Adjust grid color
    box(axs(i), 'on'); 
end
