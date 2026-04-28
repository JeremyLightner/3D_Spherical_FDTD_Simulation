%% 3D Spherical FDTD with MTLE Lightning Source
% Near-equator over-ocean ground-wave propagation
% f_max = 500 kHz, ds = lambda/20
%
% Coordinate conventions:
%   theta = polar angle from north pole (0 at N pole, pi/2 at equator)
%   phi   = azimuth measured WESTWARD from prime meridian (phi > 0 = west)
%   r     = radial distance from Earth center
%
% Yee staggering (integer index = node, half-index = cell center):
%   E_r  (i+1/2, j,     k    )  -> size (Nr-1, Nthe,   Nphi  )
%   E_th (i,     j+1/2, k    )  -> size (Nr,   Nthe-1, Nphi  )
%   E_ph (i,     j,     k+1/2)  -> size (Nr,   Nthe,   Nphi-1)
%   H_r  (i,     j+1/2, k+1/2)  -> size (Nr,   Nthe-1, Nphi-1)
%   H_th (i+1/2, j,     k+1/2)  -> size (Nr-1, Nthe,   Nphi-1)
%   H_ph (i+1/2, j+1/2, k    )  -> size (Nr-1, Nthe-1, Nphi  )
%
% Array index (i,j,k) <-> MATLAB dim (1,2,3) <-> coord (r, theta, phi)

clear; clc;

%% Physical constants
c    = 299792458;
mu0  = 4*pi*1e-7;
eps0 = 1/(c^2*mu0);
R_Earth = 6371e3;

%% Problem parameters (hardcoded for this run -- no prompts)
f_max       = 500e3;       % Hz
start_point = [0.0, 0.0];  % lat, lon (equator, prime meridian)
dist_max    = 10e3;       % m (westward propagation range)
width_max   = 1e3;         % m (N/S width)
height_max  = 8e3;         % m (top of domain above surface -- extended to fit channel)
below_ground = 2;          % m (thin seawater layer below surface; skin depth ~0.36 m)
n_max  = 1;                % maximum refractive index of environment

% Grid resolution: lambda/20 baseline
lam_min = c / (n_max*f_max);
ds      = lam_min / 20;    % ~30 m

fprintf('Grid resolution: ds = %.2f m (lambda/20 at f_max=%.0f kHz)\n', ds, f_max/1e3);

%% Spherical patch bounds
lat0 = start_point(1);
lon0 = start_point(2);

r_inner = R_Earth - below_ground;
r_outer = R_Earth + height_max;

% N/S width, centered on start lat
dlat_half_deg = rad2deg((width_max/2) / R_Earth);
lat_min = lat0 - dlat_half_deg;
lat_max = lat0 + dlat_half_deg;

% Westward range
lat0_rad = deg2rad(lat0);
if abs(cos(lat0_rad)) < 1e-12
    error('Start latitude too close to pole.');
end
dlon_range_deg = rad2deg(dist_max / (R_Earth * cos(lat0_rad)));

lon_east = lon0;
lon_west = lon0 - dlon_range_deg;

theta_min = deg2rad(90 - lat_max);
theta_max = deg2rad(90 - lat_min);
phi_start = deg2rad(-lon_east);
phi_end   = deg2rad(-lon_west);
while phi_end < phi_start
    phi_end = phi_end + 2*pi;
end

fprintf('Patch: theta=[%.6f, %.6f] rad, phi=[%.6f, %.6f] rad\n', ...
        theta_min, theta_max, phi_start, phi_end);
fprintf('       r=[%.1f, %.1f] km\n', r_inner/1e3, r_outer/1e3);

%% Grid step sizes
% Use the LARGEST radius for the angular steps so ds_physical <= ds everywhere
% (smaller r -> smaller physical arc length; we want ds to be the ceiling)
del_r   = ds;
del_the = ds / r_outer;
del_phi = ds / (r_outer * max(sin(theta_min), sin(theta_max)));

Nr   = ceil((r_outer - r_inner) / del_r)   + 1;   % number of r-nodes
Nthe = ceil((theta_max - theta_min) / del_the) + 1;
Nphi = ceil((phi_end   - phi_start ) / del_phi) + 1;

fprintf('Grid: Nr=%d, Nthe=%d, Nphi=%d  (total cells ~ %.2f M)\n', ...
        Nr, Nthe, Nphi, (Nr*Nthe*Nphi)/1e6);

% Coordinate arrays
r_node   = r_inner   + (0:Nr-1)   * del_r;       % length Nr
the_node = theta_min + (0:Nthe-1) * del_the;     % length Nthe
phi_node = phi_start + (0:Nphi-1) * del_phi;     % length Nphi

r_half   = r_node(1:end-1)   + del_r/2;          % length Nr-1
the_half = the_node(1:end-1) + del_the/2;        % length Nthe-1
phi_half = phi_node(1:end-1) + del_phi/2;        % length Nphi-1

sin_the_node = sin(the_node);
sin_the_half = sin(the_half);

%% CFL time step
% Smallest physical cell length in each direction
dl_r   = del_r;
dl_the = r_inner * del_the;                                    % smallest at r_inner
dl_phi = r_inner * min(sin_the_node) * del_phi;                % smallest at r_inner, min sin(theta)

cfl_safety = 0.95;
del_t = cfl_safety / (c * sqrt(1/dl_r^2 + 1/dl_the^2 + 1/dl_phi^2));

fprintf('Time step: dt = %.3e s (CFL safety factor %.2f)\n', del_t, cfl_safety);

% Choose total simulation time: let pulse cross the domain
t_transit = dist_max / c;
num_steps = ceil(3 * t_transit / del_t);  % 3x transit time
fprintf('Num time steps: %d (%.2f ms simulated)\n', num_steps, num_steps*del_t*1e3);

%% Material properties
% Seawater: eps_r ~ 81, sigma ~ 4 S/m
% Air:      eps_r = 1,  sigma = 0
eps_r_sea   = 81;
sigma_sea   = 4.0;
eps_r_air   = 1.0;
sigma_air   = 0.0;
mu_r_all    = 1.0;  % nonmagnetic

% Material is defined at cell centers (i+1/2, j+1/2, k+1/2), size (Nr-1, Nthe-1, Nphi-1)
% Radial index determines sea vs air: r_cell_center = r_half(i)
r_cell_center = r_half;  % (Nr-1,)
is_sea_r = r_cell_center < R_Earth;   % logical vector along r

eps_cell_r   = eps0 * (eps_r_air  + (eps_r_sea - eps_r_air) * is_sea_r);    % length Nr-1
sigma_cell_r = sigma_air + (sigma_sea - sigma_air) * is_sea_r;              % length Nr-1
mu_cell_r    = mu0 * mu_r_all * ones(1, Nr-1);

% Broadcast to full 3D cell-center arrays (uniform in theta, phi)
eps_cell   = repmat(eps_cell_r(:),   [1, Nthe-1, Nphi-1]);
sigma_cell = repmat(sigma_cell_r(:), [1, Nthe-1, Nphi-1]);
mu_cell    = repmat(mu_cell_r(:),    [1, Nthe-1, Nphi-1]);

%% Averaging material to Yee field locations
% E-field components live on cell EDGES.
% Standard 2-point averaging: average the material from the two cells sharing that edge.
%
% For our uniform-in-(theta,phi) material, this reduces to:
%   eps at E_r edge      (i+1/2, j, k):    = eps_cell at (i+1/2) -- edge lies inside one r-column, no averaging needed
%   eps at E_th edge     (i,     j+1/2,k): average of r-cells (i-1/2) and (i+1/2)
%   eps at E_phi edge    (i, j, k+1/2):    average of r-cells (i-1/2) and (i+1/2)
%
% For full 3D inhomogeneous materials you would do proper 4-cell averaging on edges.
% Keeping it simple here given layered material.

% eps for E_r: size (Nr-1, Nthe, Nphi)
eps_Er   = repmat(eps_cell_r(:),   [1, Nthe, Nphi]);
sigma_Er = repmat(sigma_cell_r(:), [1, Nthe, Nphi]);

% eps for E_th: size (Nr, Nthe-1, Nphi)
% Need averaging at r-nodes (interior). At r boundaries, just use nearest cell.
eps_cell_r_node = zeros(1, Nr);
eps_cell_r_node(1)     = eps_cell_r(1);
eps_cell_r_node(end)   = eps_cell_r(end);
eps_cell_r_node(2:end-1) = 0.5*(eps_cell_r(1:end-1) + eps_cell_r(2:end));
sigma_cell_r_node = zeros(1, Nr);
sigma_cell_r_node(1)     = sigma_cell_r(1);
sigma_cell_r_node(end)   = sigma_cell_r(end);
sigma_cell_r_node(2:end-1) = 0.5*(sigma_cell_r(1:end-1) + sigma_cell_r(2:end));

eps_Eth   = repmat(eps_cell_r_node(:),   [1, Nthe-1, Nphi]);
sigma_Eth = repmat(sigma_cell_r_node(:), [1, Nthe-1, Nphi]);

eps_Ephi   = repmat(eps_cell_r_node(:),   [1, Nthe, Nphi-1]);
sigma_Ephi = repmat(sigma_cell_r_node(:), [1, Nthe, Nphi-1]);

% mu for H components (all nonmagnetic -> just mu0)
mu_Hr   = mu0 * ones(Nr,   Nthe-1, Nphi-1);
mu_Hth  = mu0 * ones(Nr-1, Nthe,   Nphi-1);
mu_Hphi = mu0 * ones(Nr-1, Nthe-1, Nphi);

%% Precomputed update coefficients
% E updates: E^{n+1} = CEE * E^n + CEH * curlH
CErE   = (1 - sigma_Er   .* del_t ./ (2*eps_Er  )) ./ (1 + sigma_Er   .* del_t ./ (2*eps_Er  ));
CErH   = (del_t ./ eps_Er  ) ./ (1 + sigma_Er   .* del_t ./ (2*eps_Er  ));
CEthE  = (1 - sigma_Eth  .* del_t ./ (2*eps_Eth )) ./ (1 + sigma_Eth  .* del_t ./ (2*eps_Eth ));
CEthH  = (del_t ./ eps_Eth ) ./ (1 + sigma_Eth  .* del_t ./ (2*eps_Eth ));
CEphiE = (1 - sigma_Ephi .* del_t ./ (2*eps_Ephi)) ./ (1 + sigma_Ephi .* del_t ./ (2*eps_Ephi));
CEphiH = (del_t ./ eps_Ephi) ./ (1 + sigma_Ephi .* del_t ./ (2*eps_Ephi));

% H updates: H^{n+1/2} = H^{n-1/2} - CHE * curlE
CHrE   = del_t ./ mu_Hr;
CHthE  = del_t ./ mu_Hth;
CHphiE = del_t ./ mu_Hphi;

clear eps_Er sigma_Er eps_Eth sigma_Eth eps_Ephi sigma_Ephi ...
      mu_Hr mu_Hth mu_Hphi eps_cell sigma_cell mu_cell

%% Initialize fields
E_r   = zeros(Nr-1, Nthe,   Nphi  );
E_th  = zeros(Nr,   Nthe-1, Nphi  );
E_phi = zeros(Nr,   Nthe,   Nphi-1);

H_r   = zeros(Nr,   Nthe-1, Nphi-1);
H_th  = zeros(Nr-1, Nthe,   Nphi-1);
H_phi = zeros(Nr-1, Nthe-1, Nphi  );

%% Broadcast coordinate arrays for vectorized updates
r_node_1  = reshape(r_node,     [Nr   1 1]);
r_half_1  = reshape(r_half,     [Nr-1 1 1]);
sin_the_node_2 = reshape(sin_the_node, [1 Nthe   1]);
sin_the_half_2 = reshape(sin_the_half, [1 Nthe-1 1]);

%% MTLE source setup
% Channel is at the EASTERN edge of the domain: phi index k = 1
% Channel is a vertical line of cells (along r) at (k=1, some j)
% User accepted 1 km N/S truncation -- place channel at j = round(Nthe/2)
% Channel occupies r from R_Earth up to r_outer (full height above ground)

src_k = 1;                        % phi index (eastern edge)
src_j = round(Nthe/2);            % theta index (N/S center)
src_i_ground = find(r_half >= R_Earth, 1, 'first');  % first air cell
src_i_top    = Nr - 1;                               % top of domain

% Channel heights at each E_r cell center along the source line
% E_r at index (i, j, k) sits at r_half(i). We inject J_r there.
src_i_range = src_i_ground:src_i_top;
channel_z = r_half(src_i_range) - R_Earth;   % height above ocean surface

% Cross-sectional area perpendicular to r at this E_r location:
% dA = r^2 sin(theta) d_theta d_phi, evaluated at (r_half(i), the_node(j))
src_area = (r_half(src_i_range).^2) * sin_the_node(src_j) * del_the * del_phi;

% ------------------------------------------------------------------
% MTLE CURRENT FUNCTION HOOK
% Define i(0, t) and MTLE parameters here. This is the *channel-base* current.
% Replace the body of this function with your own i(0, t).
% ------------------------------------------------------------------
v_return = 1.5e8;   % return stroke speed [m/s]  (typical: 1-2e8)
lambda_decay = 2e3; % current decay constant [m] (typical: 1-2 km)

i_0 = @(t) channel_base_current(t);   % user-defined; see function at bottom

% MTLE full current: i(z, t) = i(0, t - z/v) * exp(-z/lambda)  for t >= z/v, else 0
mtle_current = @(z, t) (t >= z./v_return) .* i_0(t - z./v_return) .* exp(-z./lambda_decay);

%% Time stepping
fprintf('\nStarting time stepping...\n');
tic_total = tic;

% Probe setup: record E_r vs time at a few ranges
% Probes at 10, 25, 50, 75, 99 km from source (western direction)
probe_ranges_km = [10, 25, 50, 75, 99];
probe_k_indices = zeros(size(probe_ranges_km));
for p = 1:length(probe_ranges_km)
    % phi coordinate corresponding to range_km west of the source
    delta_phi = probe_ranges_km(p)*1e3 / (R_Earth * cos(lat0_rad));
    target_phi = phi_start + delta_phi;
    [~, probe_k_indices(p)] = min(abs(phi_node - target_phi));
end
probe_i = src_i_ground;            % just above ground
probe_j = round(Nthe/2);           % center of N/S
probe_record = zeros(num_steps, length(probe_ranges_km));

report_every = max(1, round(num_steps/50));

for n = 0:num_steps-1
    t_now = n * del_t;
    
    %% ---- H field update (uses E at time n) ----
    % H_r at (i, j+1/2, k+1/2):
    %   curlE_r = (1/(r sin(th))) * [ d(E_phi sin(th))/dth - dE_th/dphi ]
    Ephi_sin = E_phi .* sin_the_node_2;                  % (Nr, Nthe, Nphi-1)
    dEphiSin_dth = diff(Ephi_sin, 1, 2) / del_the;       % (Nr, Nthe-1, Nphi-1)
    dEth_dphi    = diff(E_th,     1, 3) / del_phi;       % (Nr, Nthe-1, Nphi-1)
    curlE_r = (1 ./ (r_node_1 .* sin_the_half_2)) .* (dEphiSin_dth - dEth_dphi);
    H_r = H_r - CHrE .* curlE_r;
    
    % H_th at (i+1/2, j, k+1/2):
    %   curlE_th = (1/r) * [ (1/sin(th)) dE_r/dphi - d(r E_phi)/dr ]
    dEr_dphi = diff(E_r, 1, 3) / del_phi;                % (Nr-1, Nthe, Nphi-1)
    rEphi    = r_node_1 .* E_phi;                        % (Nr, Nthe, Nphi-1)
    drEphi_dr = diff(rEphi, 1, 1) / del_r;               % (Nr-1, Nthe, Nphi-1)
    curlE_th = (1 ./ r_half_1) .* ((1 ./ sin_the_node_2) .* dEr_dphi - drEphi_dr);
    H_th = H_th - CHthE .* curlE_th;
    
    % H_phi at (i+1/2, j+1/2, k):
    %   curlE_phi = (1/r) * [ d(r E_th)/dr - dE_r/dth ]
    rEth   = r_node_1 .* E_th;                           % (Nr, Nthe-1, Nphi)
    drEth_dr = diff(rEth, 1, 1) / del_r;                 % (Nr-1, Nthe-1, Nphi)
    dEr_dth  = diff(E_r,  1, 2) / del_the;               % (Nr-1, Nthe-1, Nphi)
    curlE_phi = (1 ./ r_half_1) .* (drEth_dr - dEr_dth);
    H_phi = H_phi - CHphiE .* curlE_phi;
    
    %% ---- E field update (uses H at time n+1/2) ----
    % E_r at (i+1/2, j, k) interior: j=2..Nthe-1, k=2..Nphi-1
    %   curlH_r = (1/(r sin(th))) * [ d(H_phi sin(th))/dth - dH_th/dphi ]
    Hphi_sin = H_phi .* sin_the_half_2;                  % (Nr-1, Nthe-1, Nphi)
    dHphiSin_dth = diff(Hphi_sin, 1, 2) / del_the;       % (Nr-1, Nthe-2, Nphi)
    dHth_dphi    = diff(H_th,     1, 3) / del_phi;       % (Nr-1, Nthe, Nphi-2)
    
    % Align to interior region (j=2..Nthe-1, k=2..Nphi-1)
    curlH_r_int = (1 ./ (r_half_1 .* reshape(sin_the_node(2:Nthe-1), [1 Nthe-2 1]))) ...
                  .* ( dHphiSin_dth(:, :, 2:Nphi-1) - dHth_dphi(:, 2:Nthe-1, :) );
    
    E_r(:, 2:Nthe-1, 2:Nphi-1) = CErE(:, 2:Nthe-1, 2:Nphi-1) .* E_r(:, 2:Nthe-1, 2:Nphi-1) ...
                              + CErH(:, 2:Nthe-1, 2:Nphi-1) .* curlH_r_int;
    
    % E_th at (i, j+1/2, k) interior: i=2..Nr-1, k=2..Nphi-1
    %   curlH_th = (1/r) * [ (1/sin(th)) dH_r/dphi - d(r H_phi)/dr ]
    dHr_dphi = diff(H_r, 1, 3) / del_phi;                % (Nr, Nthe-1, Nphi-2)
    rHphi    = r_half_1 .* H_phi;                        % (Nr-1, Nthe-1, Nphi)
    drHphi_dr = diff(rHphi, 1, 1) / del_r;               % (Nr-2, Nthe-1, Nphi)
    
    curlH_th_int = (1 ./ reshape(r_node(2:Nr-1), [Nr-2 1 1])) ...
                   .* ( (1 ./ sin_the_half_2) .* dHr_dphi(2:Nr-1, :, :) ...
                        - drHphi_dr(:, :, 2:Nphi-1) );
    
    E_th(2:Nr-1, :, 2:Nphi-1) = CEthE(2:Nr-1, :, 2:Nphi-1) .* E_th(2:Nr-1, :, 2:Nphi-1) ...
                              + CEthH(2:Nr-1, :, 2:Nphi-1) .* curlH_th_int;
    
    % E_phi at (i, j, k+1/2) interior: i=2..Nr-1, j=2..Nthe-1
    %   curlH_phi = (1/r) * [ d(r H_th)/dr - dH_r/dth ]
    rHth   = r_half_1 .* H_th;                           % (Nr-1, Nthe, Nphi-1)
    drHth_dr = diff(rHth, 1, 1) / del_r;                 % (Nr-2, Nthe, Nphi-1)
    dHr_dth  = diff(H_r,  1, 2) / del_the;               % (Nr, Nthe-2, Nphi-1)
    
    curlH_phi_int = (1 ./ reshape(r_node(2:Nr-1), [Nr-2 1 1])) ...
                    .* ( drHth_dr(:, 2:Nthe-1, :) - dHr_dth(2:Nr-1, :, :) );
    
    E_phi(2:Nr-1, 2:Nthe-1, :) = CEphiE(2:Nr-1, 2:Nthe-1, :) .* E_phi(2:Nr-1, 2:Nthe-1, :) ...
                               + CEphiH(2:Nr-1, 2:Nthe-1, :) .* curlH_phi_int;
    
    %% ---- Source injection: MTLE distributed current ----
    % J_r = i(z, t) / cross_sectional_area, injected at E_r line (src_j, src_k)
    % E_r update with current: E_r -= (dt/eps) * J_r    (already included CErE,CErH factors,
    % but since we apply AFTER the normal update, just subtract dt*J/eps directly)
    I_channel = mtle_current(channel_z(:), t_now);       % column vector
    J_r_src = I_channel ./ src_area(:);                  % A/m^2
    
    % Use same averaging approach as CErH for consistency
    eps_src_line = eps_cell_r(src_i_range(:));
    sigma_src_line = sigma_cell_r(src_i_range(:));
    src_injection = -(del_t ./ (eps_src_line(:) + sigma_src_line(:).*del_t/2)) .* J_r_src;
    
    E_r(src_i_range, src_j, src_k) = E_r(src_i_range, src_j, src_k) + reshape(src_injection, [], 1, 1);
    
    %% ---- Boundary conditions: PEC (placeholder) ----
    % Tangential E components on outer boundary faces stay at zero (they were never updated).
    % This is REFLECTIVE -- replace with Mur ABC or PML for production use.
    % The interior update above already leaves the boundary untouched.
    
    %% ---- Record probes ----
    for p = 1:length(probe_ranges_km)
        probe_record(n+1, p) = E_r(probe_i, probe_j, probe_k_indices(p));
    end
    
    %% ---- Progress ----
    if mod(n, report_every) == 0
        elapsed = toc(tic_total);
        frac = (n+1)/num_steps;
        eta = elapsed * (1/frac - 1);
        fprintf('Step %d/%d (%.1f%%)  elapsed %.1f s  ETA %.1f s\n', ...
                n+1, num_steps, 100*frac, elapsed, eta);
    end
end

fprintf('\nTotal runtime: %.1f s\n', toc(tic_total));

%% Plot probe traces
figure;
t_axis = (0:num_steps-1) * del_t * 1e3;   % ms
for p = 1:length(probe_ranges_km)
    subplot(length(probe_ranges_km), 1, p);
    plot(t_axis, probe_record(:, p));
    ylabel(sprintf('E_r at %d km', probe_ranges_km(p)));
    grid on;
    if p == length(probe_ranges_km)
        xlabel('Time (ms)');
    end
end
sgtitle('E_r probes vs time at various ranges from source');

%% ------------------------------------------------------------------
%% USER-DEFINED CHANNEL-BASE CURRENT  i(0, t)
%% ------------------------------------------------------------------
function i = channel_base_current(t)
    % Write your own here. Example placeholder: double Heidler for first stroke.
    % Replace freely.
    I01 = 10e3;  tau11 = 0.1e-6; tau12 = 5e-6;  n1 = 2;
    I02 = 6e3;   tau21 = 2.0e-6; tau22 = 50e-6; n2 = 2;
    eta1 = exp(-(tau11/tau12)*(n1*tau12/tau11)^(1/n1));
    eta2 = exp(-(tau21/tau22)*(n2*tau21/tau22)^(1/n2));
    
    t_safe = max(t, 0);
    term1 = (I01/eta1) .* ( (t_safe./tau11).^n1 ./ (1 + (t_safe./tau11).^n1) ) .* exp(-t_safe./tau12);
    term2 = (I02/eta2) .* ( (t_safe./tau21).^n2 ./ (1 + (t_safe./tau21).^n2) ) .* exp(-t_safe./tau22);
    i = term1 + term2;
    i(t < 0) = 0;
end