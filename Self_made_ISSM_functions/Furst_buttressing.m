function [Kn, Backstress_Furst, Integrated_Backstress] = Furst_buttressing(md, vx, vy, use_flow_direction)
%--------------------------------------------------------------------------
% FURST_BUTTRESSING:  Calculate Furst Buttressing Number and Backstress
% from given velocities (steady or transient)
%
% Inputs:
%   md                  - ISSM model structure
%   vx, vy              - Velocity fields (nodal)
%   use_flow_direction  - true = flow-aligned, false = principal stress
%
% Outputs:
%   Kn                  - Furst buttressing number (0 = no buttressing)
%   Backstress_Furst    - Estimated backstress in Pa
%   Integrated_Backstress - Integrated backstress over the shelf (scalar)
%--------------------------------------------------------------------------

% Update stress from provided velocities
md = mechanicalproperties(md, vx, vy);
dev = md.results.deviatoricstress;

% Constants
g = md.constants.g;
rho_i = md.materials.rho_ice;
rho_w = md.materials.rho_water;

% Mesh info
n_el = md.mesh.numberofelements;
H = mean(md.geometry.thickness(md.mesh.elements), 2);
vx_el = mean(vx(md.mesh.elements), 2);
vy_el = mean(vy(md.mesh.elements), 2);

% Allocate
Kn = NaN(n_el, 1);
Backstress_Furst = NaN(n_el, 1);

% Loop through elements
for el = 1:n_el
    tau_xx = dev.xx(el);
    tau_yy = dev.yy(el);
    tau_xy = dev.xy(el);
    R = [2*tau_xx + tau_yy, tau_xy; tau_xy, 2*tau_yy + tau_xx];

    if use_flow_direction
        vec = [vx_el(el); vy_el(el)];
    else
        vec = dev.principalaxis2(el, :)';
    end

    if norm(vec) == 0
        continue;
    end
    n = vec / norm(vec);

    % calculate Buttressing ratio and backstress
    N = n' * R * n;
    N0 = 0.5 * g * rho_i * (1 - rho_i / rho_w) * H(el);
    Kn(el) = 1 - N / N0;
    Backstress_Furst(el) = N0 - N;
end
Backstress_Furst(Backstress_Furst <= 0) = NaN;

% calculate integrated backstress as a scalar property
Integrated_Backstress = sum(Backstress_Furst,'omitnan');
end