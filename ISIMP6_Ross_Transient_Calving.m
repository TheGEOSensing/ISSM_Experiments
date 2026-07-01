% Load the model we have derived from the inversion
cd '/Users/rishi/Desktop/ISSM_Project_Ross/Ross_Ice_Shelf/'
md = loadmodel('/Users/rishi/Desktop/ISSM_Project_Ross/Ross_Ice_Shelf/Ross_Grounded.mat');

%% ISMIP6 forcing - unchanged
md.smb = interpISMIP6AntarcticaSMB(md,'noresm1-m_rcp8.5');
md.basalforcings = interpISMIP6AntarcticaOcn(md,'noresm1-m_rcp8.5');


%% Transient components
md.transient.ismasstransport = 1;
md.transient.isstressbalance = 1;
md.transient.isgroundingline = 1;
md.transient.ismovingfront   = 1;
md.transient.isthermal       = 0;
md.inversion.iscontrol       = 0;

% Implement Von Mises calving law 
md.calving = calvingvonmises();
md.calving.stress_threshold_floatingice = 100e3; % 150 kPa - adjust to calibrate
md.calving.stress_threshold_groundedice = 1e6;   % 1 MPa - effectively no calving on grounded ice

% Level set setup (required for moving front) 
md.mask.ice_levelset        = reinitializelevelset(md, md.mask.ice_levelset);
md.levelset.stabilization   = 5;     % SUPG - most stable option
md.levelset.reinit_frequency = 10;   % reinitialize every 10 steps
md.levelset.migration_max   = 5e3;   % max front migration per year [m] - prevents blow-up

% Fix level set at domain boundary (front cannot cross model boundary)
md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
pos = find(md.mesh.vertexonboundary);
md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);

% Frontal melt (required field when ismovingfront=1) 
md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);



%% Initialize fields - unchanged
md.initialization.vx  = md.results.StressbalanceSolution.Vx;
md.initialization.vy  = md.results.StressbalanceSolution.Vy;
md.initialization.vel = md.results.StressbalanceSolution.Vel;
md.masstransport.spcthickness = NaN * ones(md.mesh.numberofvertices,1);

md.basalforcings.groundedice_melting_rate = ones(md.mesh.numberofvertices,1);
md.basalforcings.geothermalflux = zeros(md.mesh.numberofvertices,1);
md.basalforcings.islocal = 1;

%% Requested outputs - ADD calving outputs
md.transient.requested_outputs = {'default','IceVolume',...
    'IceVolumeAboveFloatation','TotalSmb','TotalFloatingBmb',...
    'AverageButtressing','CalvingCalvingrate'};


%% --- Switch to adaptive timestepping ---
% Fixed 0.5 yr step can crash when front migrates. Adaptive is safer.
md.timestepping = timesteppingadaptive();
md.timestepping.start_time    = 2014;
md.timestepping.final_time    = 2100;
md.timestepping.time_step_min = 0.01;
md.timestepping.time_step_max = 0.5;  % same upper bound as your original


%% Solve - unchanged
md.cluster = generic('name',oshostname,'np',4);
md.verbose = verbose('solution',false);
md = solve(md,'Transient');

%% Plot results
plotmodel(md, 'data', md.results.TransientSolution(1).Vel,'title#1', 'Velocity 2014 (m/yr)',...
  'data', md.results.TransientSolution(end-20).Vel,'title#2', 'Velocity 2050 years (m/yr)',...
  'caxis#1',([0 4000]),'caxis#2',([0 4000]));


%% Time-series of transient solution:
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
plot(Time_period,VAF_TS, 'LineWidth', 1.5, 'color',"#4DBEEE")
ylabel('Volume Above Floatation (m^3)', 'FontSize', 12)
ax3 = gca;

subplot(4,1,4)
plot(Time_period,Velocity_TS, 'LineWidth', 1.5, 'color',"#0072BD")
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



%% buttressing and backstress for the initial time
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
plotmodel(md,'data',Backstress_Furst(:,1),'title', 'Backstress: 2014','caxis#1',[0 3e5], ...
              'data',Backstress_Furst(:,end),'title', 'Backstress: 2100','caxis#2',[0 3e5]),...
colormap(brewermap(50,'Spectral'))

%%

plotmodel(md,'data',(Backstress_Furst(:,1) - Backstress_Furst(:,end)),'title','Difference','caxis#1', [-2e5 2e5]);
colormap(brewermap(50,'RdBu'))
