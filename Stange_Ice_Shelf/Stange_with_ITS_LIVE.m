
%% Mesh parameters
domain = '/Users/rishi/Desktop/ISSM/examples/Pig/Ross_Ice_shelf.exp';
hinit=1000;	% element size for the initial mesh
hmax=10000;		% maximum element size of the final mesh
hmin=500;		% minimum element size of the final mesh
gradation=1.7;	% maximum size ratio between two neighboring elements
err= 8;			% maximum error between interpolated and control field

% Generate an initial uniform mesh (resolution = hinit m)
md=bamg(model,'domain',domain,'hmax',hinit);
plotmodel(md,'data','mesh')

% Genreate Anisotropic Mesh based on the velocity of ice
% Get observed velocity field on mesh nodes
ncdata='/Users/rishi/Desktop/ISSM/examples/Data/ITS_LIVE_velocity_120m_RGI19A_0000_v02.nc';

x1		= ncread(ncdata,'x');
y1		= ncread(ncdata,'y'); y1 = flipud(y1);
velx	= ncread(ncdata,'vx');
vely	= ncread(ncdata,'vy');
vx		= InterpFromGridToMesh(x1,y1,velx',md.mesh.x,md.mesh.y,0);
vy		= InterpFromGridToMesh(x1,y1,vely',md.mesh.x,md.mesh.y,0);
vel     = sqrt(vx.^2+vy.^2);

%refine mesh using surface velocities as metric
md=bamg(md,'hmax',hmax,'hmin',hmin,'field',vel,'err',err);
plotmodel(md,'data','mesh')


%% 

% Step-2.1 : set the ice mask (using a level-set)
% -----------------------------------------
% read the msak from the Bedmachine dataset
mask = interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'mask','linear',...
'/Users/rishi/Desktop/ISSM/examples/Data/BedMachineAntarctica-v3.nc');


% initialize the mesh of ice mask == same as the number of vertices
md.mask.ice_levelset = ones(md.mesh.numberofvertices,1);
md.mask.ice_levelset(mask>1.5) = -1; 
% the ice mask should be >1.5 in bedmachine
md.mask.ice_levelset = reinitializelevelset(md,md.mask.ice_levelset);


%% Step-2.2: interpolate the geometry from BedMachine (bed, surface, base, thickness, grounding line)

% Surface and Bed information from BedMachine
md.geometry.surface = interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'surface','linear',...
'/Users/rishi/Desktop/ISSM/examples/Data/BedMachineAntarctica-v3.nc');

md.geometry.bed = interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'bed','linear',...
'/Users/rishi/Desktop/ISSM/examples/Data/BedMachineAntarctica-v3.nc');
% ice-shelf base same as the bed
md.geometry.base = md.geometry.bed; 


% calculation for using hydrostatic equilibrium for the ice shelf
rho_water = 1028.9;
di = md.materials.rho_ice/rho_water;
pos = find(-di/(1-di)*md.geometry.surface > md.geometry.bed);
md.geometry.base(pos) = di/(di-1)*md.geometry.surface(pos);
md.geometry.thickness = md.geometry.surface - md.geometry.base;
md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed/di;


%% Step-2.3 : interpolate observed velocities to prepare the inversion
% -----------------------------------------------------------------
% Observed velocity
% interpolate observed velocity from the velocity netcdf file
[md.inversion.vx_obs, md.inversion.vy_obs] = interpMouginotAnt2019 (md.mesh.x,md.mesh.y,...
'/Users/rishi/Desktop/ISSM/examples/Data/antarctic_ice_vel_phase_map_v01.nc');

md.inversion.vx_obs(isnan(md.inversion.vx_obs)) = 0; % make all the nan values to be zero
md.inversion.vy_obs(isnan(md.inversion.vy_obs)) = 0;
md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2); % compute velocity from vx and vy



%% Step-2.4 : initialize friction parameters and rheology parameters
% -----------------------------------------------------------------

% initialize Sliding parameters (uniform 200 by default)
md.friction.coefficient = 200 * ones(md.mesh.numberofvertices, 1);
md.friction.p = ones(md.mesh.numberofelements, 1);
md.friction.q = ones(md.mesh.numberofelements, 1);

% Initialize Rheology parameters, assuming ice is at -10°C
md.materials.rheology_B = cuffey(273.15-10) * ones(md.mesh.numberofvertices, 1);
md.materials.rheology_n = 3 * ones(md.mesh.numberofelements, 1);

% Other parameters (just ignore)
md.stressbalance.referential = NaN(md.mesh.numberofvertices,6);
md.stressbalance.loadingforce = zeros(md.mesh.numberofvertices,3);



%% Step-2.5 : Set Boundary conditions
% ----------------------------------

% Set Boundary conditions: constrain inflow velocities only
md.stressbalance.spcvx = NaN(md.mesh.numberofvertices,1);
md.stressbalance.spcvy = NaN(md.mesh.numberofvertices,1);
md.stressbalance.spcvz = NaN(md.mesh.numberofvertices,1);
pos = find(md.mesh.vertexonboundary & md.mask.ocean_levelset>0);
md.stressbalance.spcvx(pos) = md.inversion.vx_obs(pos);
md.stressbalance.spcvy(pos) = md.inversion.vy_obs(pos);



%% Step-2.6 : Set flowequation as SSA
% ---------------------------------- 
md = setflowequation(md,'SSA','all');

% Model name
md.miscellaneous.name = 'Stange_Inversion';
%save ./Models/Stange_Parameterized md;

% Plot some parameters
plotmodel(md,'data','BC', ... % boundary conditions
    'data',md.inversion.vel_obs,'log#2',10,'caxis#2',[0.1 4000],... % observed velocity 
    'data',md.mask.ice_levelset,'caxis#3',[-1 1], ... % masked ice and no ice
    'data',md.mask.ocean_levelset,'caxis#4',[-1 1]) % masked ice and ocean


%%  Step 3: Infer ice shelf rheology

% Control general
md.inversion= m1qn3inversion(md.inversion);
md.inversion.iscontrol = 1; % controlled experiments
md.inversion.maxsteps = 80; % maximum steps
md.inversion.maxiter = 80;  % maximum interation
md.inversion.dxmin = 0.1; 
md.inversion.gttol = 1.0e-6;

% Cost functions
md.inversion.cost_functions = [101 103 502];
md.inversion.cost_functions_coefficients= ones(md.mesh.numberofvertices,3);
md.inversion.cost_functions_coefficients(:,1) = 1;
md.inversion.cost_functions_coefficients(:,2) = 0.01;
md.inversion.cost_functions_coefficients(:,3) = 1e-18;
pos=find(md.inversion.vx_obs==0 | md.mask.ice_levelset>0);
md.inversion.cost_functions_coefficients(pos,1:2)=0;

% Controls
md.inversion.control_parameters={'MaterialsRheologyBbar'};
md.inversion.min_parameters = md.materials.rheology_B;
md.inversion.max_parameters = md.materials.rheology_B;
pos = find(md.mask.ocean_levelset<0);
md.inversion.min_parameters(pos) = cuffey(273);
md.inversion.max_parameters(pos) = cuffey(200);


% Solve
md.verbose=verbose('control',true);
md.cluster=generic('name',oshostname,'np',2);
mds = extract(md,md.mask.ocean_levelset<0); % only for the ice shelf
mds = solve(mds,'Stressbalance');
%save ./Models/Stange_Inverted mds; % save the model at specific location

%%

% Plot the results of the inversion

plotmodel(mds,...
'data',mds.results.StressbalanceSolution.Vel,'title','Modeled velocity','colormap','turbo',...
'data',mds.inversion.vel_obs,'title','Observed velocity','colormap','turbo',...
'data',mds.results.StressbalanceSolution.MaterialsRheologyBbar,'title','Infered rheology',...
'data',mds.results.StressbalanceSolution.Vel-mds.inversion.vel_obs,'title','Velocity difference',...
'caxis#1',[0 4000],'caxis#2',[0 4000],'caxis#4',[-100 100])