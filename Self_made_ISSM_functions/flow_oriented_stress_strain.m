function md = flow_oriented_stress_strain(md, vx, vy)
%--------------------------------------------------------------------------
% FLOW_ORIENTED_STRESS_STRAIN: Computes flow-aligned strain and stress
% components: longitudinal, transverse, and shear (strainrate + deviatoric stress).
%
% Inputs:
%   md      - ISSM model structure (with strainrate and deviatoric stress)
%   vx, vy  - Nodal velocity components (from solution or transient step)
%
% Output:
%   md      - Updated model structure with new results:
%             md.results.strainrate.elon, .etrans, .eshear
%             md.results.deviatoricstress.longitudinal, .transverse, .shear
%--------------------------------------------------------------------------
% Example Case
%--------------------------------------------------------------------------
% vx = md.results.TransientSolution(1).Vx;
% vy = md.results.TransientSolution(1).Vy;
% 
% md = flow_oriented_stress_strain(md, vx, vy);
%--------------------------------------------------------------------------

% Update deviatoric stress from velocities
md = mechanicalproperties(md, vx, vy);

% Number of elements
n_el = md.mesh.numberofelements;

% Element-averaged velocities
vx_el = mean(vx(md.mesh.elements), 2);
vy_el = mean(vy(md.mesh.elements), 2);

% Allocate outputs
elon   = NaN(n_el, 1); etrans = NaN(n_el, 1); eshear = NaN(n_el, 1);
slon   = NaN(n_el, 1); strans = NaN(n_el, 1); sshear = NaN(n_el, 1);

% Loop through elements
for el = 1:n_el
    % --- Strain rate tensor ---
    e_xx = md.results.strainrate.xx(el);
    e_yy = md.results.strainrate.yy(el);
    e_xy = md.results.strainrate.xy(el);
    Rs = [e_xx, e_xy;
          e_xy, e_yy];

    % --- Stress tensor ---
    sig_xx = md.results.deviatoricstress.xx(el);
    sig_yy = md.results.deviatoricstress.yy(el);
    sig_xy = md.results.deviatoricstress.xy(el);
    S = [sig_xx, sig_xy;
         sig_xy, sig_yy];

    % --- Flow direction ---
    vec = [vx_el(el); vy_el(el)];
    if norm(vec) == 0
        continue;
    end
    n = vec / norm(vec);           % unit vector along flow
    n_perp = [-n(2); n(1)];         % orthogonal unit vector

    % --- Flow-oriented strain rate components ---
    elon(el)   = n'      * Rs * n;
    etrans(el) = n_perp' * Rs * n_perp;
    eshear(el) = n'      * Rs * n_perp;

    % --- Flow-oriented stress components ---
    slon(el)   = n'      * S * n;
    strans(el) = n_perp' * S * n_perp;
    sshear(el) = n'      * S * n_perp;
end

% Store in structure
md.results.strainrate.elon   = elon;
md.results.strainrate.etrans = etrans;
md.results.strainrate.eshear = eshear;

md.results.deviatoricstress.longitudinal   = slon;
md.results.deviatoricstress.transverse = strans;
md.results.deviatoricstress.shear = sshear;

end