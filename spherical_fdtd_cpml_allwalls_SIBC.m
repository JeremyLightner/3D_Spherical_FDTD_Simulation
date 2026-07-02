function spherical_fdtd_cpml_allwalls_SIBC()
%% 3D Spherical FDTD with MTLE Lightning Source -- Multi-wall CPML
%
% IMPORTANT: This file is a FUNCTION, not a script (required for nested
% functions to share workspace). Run by calling spherical_fdtd_cpml_allwalls.
%
% CPML on west, east, north, south, and top walls. Bottom is PEC (ocean
% stand-in). Architecture: PER-COORDINATE PROFILES + PML-SIZED psi arrays.
%
%   - Each coordinate (r, theta, phi) has profile arrays (kappa, sigma,
%     alpha -> b, a) defined along its length, with PML ramps at the ends
%     that have nonzero thickness.
%   - Every curl derivative is stretched along its OWN differentiation
%     direction using that coordinate's profile.
%   - Corners and edges are handled AUTOMATICALLY: a corner cell is in the
%     PML index set of multiple coordinates, so each of its derivatives
%     gets the appropriate stretching with no special-case code.
%
% ENABLING / DISABLING WALLS:
%   Set the per-wall thickness (cpml_ct_*) to 0 to disable that wall, or
%   to a positive integer (e.g. 10) to enable it. No other code changes.
%   ** For first testing, enable ONLY north+south (set top=east=0), verify,
%      then enable east, then top. See notes at end of file. **
%
% Coordinate conventions:
%   theta = 0 at north pole, pi/2 at equator
%   phi   = 0 at prime meridian, increases WESTWARD

clc;

%% ============================================================
%% Physical constants
%% ============================================================
c    = 299792458;
mu0  = 4*pi*1e-7;
eps0 = 1/(c^2*mu0);
R_Earth = 6371e3;

%% ============================================================
%% Problem parameters
%% ============================================================
f_max        = 500e3;
start_point  = [0.0, 0.0];
dist_max     = 40e3;
width_max    = 7e3;
height_max   = 5e3;
below_ground = 0; % need to be careful to ensure that this doesnt put the surface of the earth in the wrong spot

lam_min = c / f_max;
ds      = lam_min / 20;
fprintf('Grid resolution: ds = %.2f m (lambda/20 at f_max=%.0f kHz)\n', ds, f_max/1e3);

%% ============================================================
%% Geometry
%% ============================================================
% Assume Spherical Earth (Close enough)
lat0 = start_point(1);
lon0 = start_point(2);
r_inner = R_Earth - below_ground;
r_outer = R_Earth + height_max;

dlat_half_deg = rad2deg((width_max/2) / R_Earth); % Half arc length converted
lat_min = lat0 - dlat_half_deg; % Centered on start point
lat_max = lat0 + dlat_half_deg;

lat0_rad = deg2rad(lat0);
if abs(cos(lat0_rad)) < 1e-12
    error('Start latitude too close to pole.');
end
dlon_range_deg = rad2deg(dist_max / (R_Earth * cos(lat0_rad))); % Modified arc length to account for latitude
lon_east = lon0; % Gridding established from eastern side but source can be anywhere
lon_west = lon0 - dlon_range_deg;

theta_min = deg2rad(90 - lat_max);
theta_max = deg2rad(90 - lat_min);
phi_start = deg2rad(-lon_east);
phi_end   = deg2rad(-lon_west);
while phi_end < phi_start
    phi_end = phi_end + 2*pi;
end

%% ============================================================
%% Grid step sizes and CFL
%% ============================================================
del_r   = ds;
del_the = ds / r_outer;
del_phi = ds / (r_outer * max(sin(theta_min), sin(theta_max)));

Nr   = ceil((r_outer - r_inner) / del_r)   + 1;
Nthe = ceil((theta_max - theta_min) / del_the) + 1;
Nphi = ceil((phi_end - phi_start) / del_phi)   + 1;
fprintf('Grid: Nr=%d, Nthe=%d, Nphi=%d (total cells ~ %.2f M)\n', ...
        Nr, Nthe, Nphi, (Nr*Nthe*Nphi)/1e6);

r_node   = r_inner   + (0:Nr-1)   * del_r;
the_node = theta_min + (0:Nthe-1) * del_the;
phi_node = phi_start + (0:Nphi-1) * del_phi;
r_half   = r_node(1:end-1)   + del_r/2;
the_half = the_node(1:end-1) + del_the/2;
phi_half = phi_node(1:end-1) + del_phi/2;
sin_the_node = sin(the_node);
sin_the_half = sin(the_half);

dl_r   = del_r;
dl_the = r_inner * del_the;
dl_phi = r_inner * min(sin_the_node) * del_phi;
cfl_safety = 0.95;
del_t = cfl_safety / (c * sqrt(1/dl_r^2 + 1/dl_the^2 + 1/dl_phi^2));
fprintf('Time step: dt = %.3e s (CFL safety factor %.2f)\n', del_t, cfl_safety);

t_transit = dist_max / c;
num_steps = ceil(3 * t_transit / del_t);
fprintf('Num time steps: %d (%.2f ms simulated)\n', num_steps, num_steps*del_t*1e3);

%% ============================================================
%% Broadcast helpers
%% ============================================================
r_node_1       = reshape(r_node,        [Nr   1 1]);
r_half_1       = reshape(r_half,        [Nr-1 1 1]);
sin_the_node_2 = reshape(sin_the_node,  [1 Nthe   1]);
sin_the_half_2 = reshape(sin_the_half,  [1 Nthe-1 1]);

%% ============================================================
%% Material properties (seawater below R_Earth, air above)
%% ============================================================
eps_r_sea = 81;  sigma_sea = 4.0;
eps_r_air = 1.0; sigma_air = 0.0;
mu_r_all  = 1.0;

is_sea_r = r_half < R_Earth;
eps_cell_r   = eps0 * (eps_r_air + (eps_r_sea - eps_r_air) * is_sea_r);
sigma_cell_r = sigma_air + (sigma_sea - sigma_air) * is_sea_r;

eps_cell_r_node = zeros(1, Nr);
eps_cell_r_node(1)       = eps_cell_r(1);
eps_cell_r_node(end)     = eps_cell_r(end);
eps_cell_r_node(2:end-1) = 0.5*(eps_cell_r(1:end-1) + eps_cell_r(2:end));
sigma_cell_r_node = zeros(1, Nr);
sigma_cell_r_node(1)       = sigma_cell_r(1);
sigma_cell_r_node(end)     = sigma_cell_r(end);
sigma_cell_r_node(2:end-1) = 0.5*(sigma_cell_r(1:end-1) + sigma_cell_r(2:end));

eps_Er    = repmat(eps_cell_r(:),        [1, Nthe, Nphi]);
sigma_Er  = repmat(sigma_cell_r(:),      [1, Nthe, Nphi]);
eps_Eth   = repmat(eps_cell_r_node(:),   [1, Nthe-1, Nphi]);
sigma_Eth = repmat(sigma_cell_r_node(:), [1, Nthe-1, Nphi]);
eps_Ephi  = repmat(eps_cell_r_node(:),   [1, Nthe, Nphi-1]);
sigma_Ephi= repmat(sigma_cell_r_node(:), [1, Nthe, Nphi-1]);

CErE   = (1 - sigma_Er   .* del_t ./ (2*eps_Er  )) ./ (1 + sigma_Er   .* del_t ./ (2*eps_Er  ));
CErH   = (del_t ./ eps_Er  ) ./ (1 + sigma_Er   .* del_t ./ (2*eps_Er  ));
CEthE  = (1 - sigma_Eth  .* del_t ./ (2*eps_Eth )) ./ (1 + sigma_Eth  .* del_t ./ (2*eps_Eth ));
CEthH  = (del_t ./ eps_Eth ) ./ (1 + sigma_Eth  .* del_t ./ (2*eps_Eth ));
CEphiE = (1 - sigma_Ephi .* del_t ./ (2*eps_Ephi)) ./ (1 + sigma_Ephi .* del_t ./ (2*eps_Ephi));
CEphiH = (del_t ./ eps_Ephi) ./ (1 + sigma_Ephi .* del_t ./ (2*eps_Ephi));
CHrE   = del_t / (mu0 * mu_r_all);
CHthE  = del_t / (mu0 * mu_r_all);
CHphiE = del_t / (mu0 * mu_r_all);

clear eps_Er sigma_Er eps_Eth sigma_Eth eps_Ephi sigma_Ephi

%% ============================================================
%% Field allocations
%% ============================================================
E_r   = zeros(Nr-1, Nthe,   Nphi  );
E_th  = zeros(Nr,   Nthe-1, Nphi  );
E_phi = zeros(Nr,   Nthe,   Nphi-1);
H_r   = zeros(Nr,   Nthe-1, Nphi-1);
H_th  = zeros(Nr-1, Nthe,   Nphi-1);
H_phi = zeros(Nr-1, Nthe-1, Nphi  );

%% ============================================================
%% CPML setup -- per-coordinate profiles, all walls
%% ============================================================
% Profile parameters (tuned: alpha_max scaled for low frequency, see notes)
cpml_m         = 3;
cpml_kappa_max = 2.07;  % 2.07 comes from Perry Optimization (should be max 4)
cpml_alpha_max = 2*pi*f_max*eps0;     % ~ omega*eps0 at f_max; tune per your band (2.7816e-5)

% Per-wall thicknesses (cells). Set to 0 to DISABLE that wall.
% phi: east = low-index end, west = high-index end
% theta: north = low-index end, south = high-index end
% r: top = high-index end (bottom is PEC, always 0)
cpml_ct_west  = 10;
cpml_ct_east  = 10;
cpml_ct_north = 10;
cpml_ct_south = 10;
cpml_ct_top   = 10;
cpml_ct_bottom = 0;   % PEC ocean stand-in -- keep 0

% sigma_max per coordinate (Roden-Gedney optimal, using representative cell size)
% Representative physical cell sizes (at surface / equator):
dphys_r   = del_r;
dphys_the = R_Earth * del_the;
dphys_phi = R_Earth * sin(pi/2) * del_phi;
% sigma_max_r   = (cpml_m+1) / (150*pi*dphys_r);
% sigma_max_the = (cpml_m+1) / (150*pi*dphys_the);
% sigma_max_phi = (cpml_m+1) / (150*pi*dphys_phi);
sigma_max_r   = - 0.25 * 5 * (cpml_m+1) * log(10E-6) / (2*120*pi*dphys_r);
sigma_max_the = - 0.25 * (cpml_m+1) * log(10E-6) / (2*120*pi*dphys_the);
sigma_max_phi = - 0.25 * (cpml_m+1) * log(10E-6) / (2*120*pi*dphys_phi);

% disp(sigma_max_r)
% disp(sigma_max_the)
% disp(sigma_max_phi)

% Build per-coordinate profiles (both staggerings, ramps at enabled ends)
[b_r_int,   a_r_int,   kap_r_int,   b_r_half,   a_r_half,   kap_r_half] = ...
    build_coord_profiles(Nr,   cpml_ct_bottom, cpml_ct_top, ...
                         cpml_m, cpml_kappa_max, cpml_alpha_max, sigma_max_r,   del_t, eps0);
[b_the_int, a_the_int, kap_the_int, b_the_half, a_the_half, kap_the_half] = ...
    build_coord_profiles(Nthe, cpml_ct_north,  cpml_ct_south, ...
                         cpml_m, cpml_kappa_max, cpml_alpha_max, sigma_max_the, del_t, eps0);
[b_phi_int, a_phi_int, kap_phi_int, b_phi_half, a_phi_half, kap_phi_half] = ...
    build_coord_profiles(Nphi, cpml_ct_east,   cpml_ct_west, ...
                         cpml_m, cpml_kappa_max, cpml_alpha_max, sigma_max_phi, del_t, eps0);

% ---- PML index sets ----
% "half" sets: indices into a diff array that is at half-positions along the
%   differenced dim (i.e. diff of a node-staggered field). Used in H updates.
% "node-diff" sets: indices into a diff array at interior-node positions
%   (i.e. diff of a half-staggered field). Used in E updates. The profile is
%   sampled at node position (diff index + 1).

% phi
idx_phi_half = [1:cpml_ct_east, (Nphi-cpml_ct_west):(Nphi-1)];
idx_phi_half = idx_phi_half(idx_phi_half>=1 & idx_phi_half<=Nphi-1);
idx_phi_ndf  = [1:(cpml_ct_east-1), (Nphi-cpml_ct_west):(Nphi-2)];
idx_phi_ndf  = idx_phi_ndf(idx_phi_ndf>=1 & idx_phi_ndf<=Nphi-2);
pos_phi_ndf  = idx_phi_ndf + 1;   % node positions for profile sampling

% theta
idx_the_half = [1:cpml_ct_north, (Nthe-cpml_ct_south):(Nthe-1)];
idx_the_half = idx_the_half(idx_the_half>=1 & idx_the_half<=Nthe-1);
idx_the_ndf  = [1:(cpml_ct_north-1), (Nthe-cpml_ct_south):(Nthe-2)];
idx_the_ndf  = idx_the_ndf(idx_the_ndf>=1 & idx_the_ndf<=Nthe-2);
pos_the_ndf  = idx_the_ndf + 1;

% r (top only)
idx_r_half = (Nr-cpml_ct_top):(Nr-1);
idx_r_half = idx_r_half(idx_r_half>=1 & idx_r_half<=Nr-1);
idx_r_ndf  = (Nr-cpml_ct_top):(Nr-2);
idx_r_ndf  = idx_r_ndf(idx_r_ndf>=1 & idx_r_ndf<=Nr-2);
pos_r_ndf  = idx_r_ndf + 1;

% ---- Reshaped coefficient vectors (broadcast along the differenced dim) ----
% phi -> dim 3 -> [1 1 L]
bph_h = reshape(b_phi_half(idx_phi_half), [1 1 numel(idx_phi_half)]);
aph_h = reshape(a_phi_half(idx_phi_half), [1 1 numel(idx_phi_half)]);
kph_h = reshape(kap_phi_half(idx_phi_half), [1 1 numel(idx_phi_half)]);
bph_n = reshape(b_phi_int(pos_phi_ndf), [1 1 numel(idx_phi_ndf)]);
aph_n = reshape(a_phi_int(pos_phi_ndf), [1 1 numel(idx_phi_ndf)]);
kph_n = reshape(kap_phi_int(pos_phi_ndf), [1 1 numel(idx_phi_ndf)]);

% theta -> dim 2 -> [1 L 1]
bth_h = reshape(b_the_half(idx_the_half), [1 numel(idx_the_half) 1]);
ath_h = reshape(a_the_half(idx_the_half), [1 numel(idx_the_half) 1]);
kth_h = reshape(kap_the_half(idx_the_half), [1 numel(idx_the_half) 1]);
bth_n = reshape(b_the_int(pos_the_ndf), [1 numel(idx_the_ndf) 1]);
ath_n = reshape(a_the_int(pos_the_ndf), [1 numel(idx_the_ndf) 1]);
kth_n = reshape(kap_the_int(pos_the_ndf), [1 numel(idx_the_ndf) 1]);

% r -> dim 1 -> [L 1 1]
br_h = reshape(b_r_half(idx_r_half), [numel(idx_r_half) 1 1]);
ar_h = reshape(a_r_half(idx_r_half), [numel(idx_r_half) 1 1]);
kr_h = reshape(kap_r_half(idx_r_half), [numel(idx_r_half) 1 1]);
br_n = reshape(b_r_int(pos_r_ndf), [numel(idx_r_ndf) 1 1]);
ar_n = reshape(a_r_int(pos_r_ndf), [numel(idx_r_ndf) 1 1]);
kr_n = reshape(kap_r_int(pos_r_ndf), [numel(idx_r_ndf) 1 1]);

% ---- psi allocations (PML-sized along the differenced dimension) ----
Lph_h = numel(idx_phi_half);  Lph_n = numel(idx_phi_ndf);
Lth_h = numel(idx_the_half);  Lth_n = numel(idx_the_ndf);
Lr_h  = numel(idx_r_half);    Lr_n  = numel(idx_r_ndf);

% H-update derivatives (of E)
psi_Eth_phi    = zeros(Nr,   Nthe-1, Lph_h);   % dE_th/dphi      at H_r stagger
psi_Ephisin_th = zeros(Nr,   Lth_h,  Nphi-1);  % d(E_phi sinth)/dth at H_r
psi_Er_phi     = zeros(Nr-1, Nthe,   Lph_h);   % dE_r/dphi       at H_th
psi_rEphi_r    = zeros(Lr_h, Nthe,   Nphi-1);  % d(r E_phi)/dr   at H_th
psi_rEth_r     = zeros(Lr_h, Nthe-1, Nphi  );  % d(r E_th)/dr    at H_phi
psi_Er_th      = zeros(Nr-1, Lth_h,  Nphi  );  % dE_r/dth        at H_phi

% E-update derivatives (of H)
psi_Hphisin_th = zeros(Nr-1, Lth_n,  Nphi  );  % d(H_phi sinth)/dth at E_r
psi_Hth_phi    = zeros(Nr-1, Nthe,   Lph_n);   % dH_th/dphi      at E_r
psi_Hr_phi     = zeros(Nr,   Nthe-1, Lph_n);   % dH_r/dphi       at E_th
psi_rHphi_r    = zeros(Lr_n, Nthe-1, Nphi  );  % d(r H_phi)/dr   at E_th
psi_rHth_r     = zeros(Lr_n, Nthe,   Nphi-1);  % d(r H_th)/dr    at E_phi
psi_Hr_th      = zeros(Nr,   Lth_n,  Nphi-1);  % dH_r/dth        at E_phi

fprintf('CPML enabled walls: W=%d E=%d N=%d S=%d Top=%d (cells)\n', ...
        cpml_ct_west, cpml_ct_east, cpml_ct_north, cpml_ct_south, cpml_ct_top);
%% ============================================================
%% TEST SOURCE: charge-neutral soft source (differentiated Gaussian)
%% ============================================================
% Use this INSTEAD of the MTLE source for validating propagation/CPML/SIBC.
% Zero-time-integral waveform -> no charge accumulation, no static field.
test_source_enabled = false;     % set false to use the MTLE source instead

% Center the spectrum in your band of interest:
test_f0   = 200e3;              % spectral peak [Hz] (in your 100-300 kHz band)
test_tau  = 1/(2*pi*test_f0);  % Gaussian width giving that peak
test_t0   = 6*test_tau;        % delay so the pulse starts ~from zero (6 sigma)
test_amp  = 1.0;               % amplitude scaling (arbitrary for a test)

% Where to inject (a single point, a few cells off the boundaries/source col):
test_i = round(Nr/2);          % mid-height (well away from top/bottom)
test_j = round(Nthe/2);        % mid-theta
test_k = round(Nphi/2);        % mid-phi (well inside, away from CPML walls)
if test_source_enabled
    fprintf('TEST soft source: f0=%.0f kHz, tau=%.3e s, at (i,j,k)=(%d,%d,%d)\n', ...
        test_f0/1e3, test_tau, test_i, test_j, test_k);
end

%% ============================================================
%% Surface Impedance Boundary Condition (SIBC) -- contour-path / Beggs form
%% ============================================================
% Replaces the PEC ocean bottom (i=1) with a Leontovich impedance closure,
% implemented by substituting the SIBC into the discrete Faraday update for
% the surface tangential H (H_th, H_phi at the first half-plane r_half(1)).
% This makes the surface-H update IMPLICIT in H (H appears on both sides),
% which is solved per-cell by a scalar division -- no lag, no extrapolation.
%
% Derivation summary (H_phi shown; H_th identical by symmetry):
%   Faraday (phi): dH_phi/dt = -(1/mu0)(1/r_h1)[ (r2 E_th(2) - r1 E_th(1))/dr
%                                                 - dE_r/dth ]
%   The lower theta-edge E_th(1) sits ON the surface (r=R_e); substitute the
%   monochromatic SIBC, time-centered at n (your correction):
%     E_th(1)^n = -R_s[ (H_phi^{n+1/2}+H_phi^{n-1/2})/2
%                       + (1/w0)(H_phi^{n+1/2}-H_phi^{n-1/2})/dt ]
%   Collecting H_phi^{n+1/2} gives the implicit update below. The +Z_s for
%   E_phi (in the H_th contour) combines with the opposite curl sign to yield
%   the SAME coefficients, so one (delta,c1,c2) triple serves both components.
%
%   H^{n+1/2}(1) = ( H_explicit(1) + c1 * H^{n-1/2}(1) ) / c2
%     delta = dt*R_s*r_node(1) / (mu0 * r_half(1) * dr)
%     c1    = delta*( 1/(w0*dt) - 1/2 )
%     c2    = 1 + delta/2 + delta/(w0*dt)
%
% PEC limit: sigma->inf => R_s->0 => delta->0 => c1->0, c2->1 => H unchanged.
%
% << CONV UPGRADE >>: for the dispersive variant, the instantaneous part of the
% Z_s convolution is absorbed into c2 and the history part replaces c1*H^{n-1/2}
% with a sum of recursive accumulators Psi_k (Psi_k <- Q_k Psi_k + S_k H). The
% skeleton "explicit update, then (.+history)/c2 on the i=1 rows" is unchanged.

sibc_enabled = false;       % set true to enable SIBC (false = PEC bottom)
sibc_mode    = 'conv';      % 'mono' (single-freq) or 'conv' (dispersive)   <== NEW
sibc_f0      = 200e3;       % monochromatic design frequency [Hz] ('mono' only)
sibc_sigma   = 4.0;         % seawater conductivity [S/m] (match GRWAVE run)

% ---- geometric coupling (shared by both modes) ----   <== NEW (shared)
% g = dt * r_node(1) / (mu0 * r_half(1) * dr); coefficient on E_surf in the
% Faraday correction  H(1) += g * E_surf.
sibc_g = del_t * r_node(1) / (mu0 * r_half(1) * del_r);

% ---- monochromatic coefficients ('mono' mode) ----   (unchanged math)
sibc_w0 = 2*pi*sibc_f0;
sibc_Rs = sqrt(sibc_w0 * mu0 / (2*sibc_sigma));          % surface resistance [ohm]
sibc_delta = sibc_g * sibc_Rs;                           % = dt*Rs*r1/(mu0*rh1*dr)
sibc_c1 = sibc_delta * ( 1/(sibc_w0*del_t) - 0.5 );
sibc_c2 = 1 + sibc_delta/2 + sibc_delta/(sibc_w0*del_t);

% ================= CONVOLUTIONAL SIBC ('conv' mode) =====================  <== NEW BLOCK
% Dispersive surface impedance via recursive convolution. Z_s(s) is fit by a
% rational (pole/residue) function -- from Vector Fitting (vectfit3.m) of the
% good-conductor Z_s(w)=sqrt(j w mu0/sigma) over your band (100-300 kHz + low
% tail). Each pole gives one recursive accumulator (like a CPML psi):
%   Psi_k^n = Q_k Psi_k^{n-1} + S_k H^{n-1/2}      (causal -> EXPLICIT)
%   E_surf^n = -[ Z_inf H^n + sum_k Psi_k^n ]
%   Q_k = exp(p_k dt),  S_k = c_k (exp(p_k dt) - 1)/p_k
% Poles p_k must be in the LEFT half-plane (Re(p_k)<0) for stability/passivity.
%
%   <<< REPLACE THESE PLACEHOLDER POLES WITH YOUR vectfit3 OUTPUT >>>
% Column vectors, conjugate pairs included so the accumulator sum is real. The
% placeholder is a crude 1-pole stand-in ONLY so the code runs; it is NOT an
% accurate fit of sqrt(jw) and WILL give wrong attenuation.
sibc_poles = [ -2*pi*200e3 ];            % p_k [rad/s] (Np x 1, complex)  PLACEHOLDER
sibc_res   = [  sibc_Rs * 2*pi*200e3 ];  % c_k [ohm*rad/s]                PLACEHOLDER
sibc_Zinf  = 0;                          % constant term; 0 for sqrt(s)
sibc_poles_are_placeholder = true;       % set false once a real fit is inserted

sibc_Np = numel(sibc_poles);
sibc_Q  = exp(sibc_poles(:) * del_t);                                       % (Np x 1)
sibc_S  = sibc_res(:) .* (exp(sibc_poles(:)*del_t) - 1) ./ sibc_poles(:);   % (Np x 1)
sibc_Qb = reshape(sibc_Q, [sibc_Np 1 1]);   % broadcast shapes for accumulators
sibc_Sb = reshape(sibc_S, [sibc_Np 1 1]);
sibc_c2conv = 1 + sibc_g * sibc_Zinf / 2;   % implicit denom; =1 when Zinf=0

% Pole accumulators on the surface sheet (Np x tangential face).  <== NEW STATE
Psi_phi = zeros(sibc_Np, Nthe-1, Nphi  );   % drives E_th  (paired with H_phi)
Psi_th  = zeros(sibc_Np, Nthe,   Nphi-1);   % drives E_phi (paired with H_th )
% ========================================================================

% Previous-step surface H (= H^{n-1/2}); rest initial condition -> zero.
Hth_prev1  = zeros(1, Nthe,   Nphi-1);   % matches H_th(1,:,:)
Hphi_prev1 = zeros(1, Nthe-1, Nphi  );   % matches H_phi(1,:,:)

% Optional: store the diagnostic surface tangential E (NOT fed back into E_th/
% E_phi(1,:,:), which must stay zero or update_H would double-count).
sibc_store_Esurf = false;
if sibc_store_Esurf
    Eth_surf_diag  = zeros(1, Nthe-1, Nphi  );
    Ephi_surf_diag = zeros(1, Nthe,   Nphi-1);
end

if sibc_enabled
    fprintf('SIBC (contour) mode=%s  sigma=%.1f S/m  g=%.3e\n', sibc_mode, sibc_sigma, sibc_g);
    if strcmp(sibc_mode, 'mono')
        fprintf('   mono: f0=%.0f kHz  R_s=%.4f ohm  c1=%.3e c2=%.5f\n', ...
                sibc_f0/1e3, sibc_Rs, sibc_c1, sibc_c2);
    else
        fprintf('   conv: Np=%d poles  Zinf=%.4g  c2conv=%.5f\n', sibc_Np, sibc_Zinf, sibc_c2conv);
        if any(real(sibc_poles(:)) >= 0)
            warning('SIBC conv: pole(s) NOT in left half-plane -> UNSTABLE.');
        end
        if sibc_poles_are_placeholder
            warning(['SIBC conv: PLACEHOLDER poles -- attenuation WILL be wrong. ' ...
                     'Replace sibc_poles/sibc_res with a vectfit3 fit of Z_s(w).']);
        end
    end
end

%% ============================================================
%% MTLE source setup
%% ============================================================
% If east CPML is enabled, the source must sit INSIDE the domain, not in the
% PML. Inset it a few cells past the east PML inner edge.

src_j        = round(Nthe/2);
src_i_ground = find(r_half >= R_Earth, 1, 'first');
% src_i_top    = Nr - 1;
% src_i_top    = floor(Nr/2);
src_i_top    = round(1000/del_r) + src_i_ground; % 1000 is determined by the height of a -CG
src_i_range  = src_i_ground:src_i_top;
src_inset    = round((2*((src_i_top*del_r)^2)/lam_min) / del_r) + 1; % Set by the Fraunhofer 
disp(src_inset)
    % Distance to capture near field effects given a dr and minimum wavelength
if cpml_ct_east > 0
    src_k = cpml_ct_east + src_inset;
else
    src_k = src_inset;
end

fprintf('Source at phi-index k=%d (inset %d past east PML)\n', src_k, src_inset);

channel_z = r_half(src_i_range) - R_Earth;
src_area  = (r_half(src_i_range).^2) * sin_the_node(src_j) * del_the * del_phi;

v_return     = (7/12)*c;
lambda_decay = 1e3;
i_0 = @(t) channel_base_current(t);
mtle_current = @(z, t) (t >= z./v_return) .* i_0(t - z./v_return) .* exp(-z./lambda_decay);


%% ============================================================
%% Probe setup
%% ============================================================
probe_ranges_km = [1, 5, 10, 20, 38];
probe_k_indices = zeros(size(probe_ranges_km));
for p = 1:length(probe_ranges_km)
    delta_phi = probe_ranges_km(p)*1e3 / (R_Earth * cos(lat0_rad));
    [~, probe_k_indices(p)] = min(abs(phi_node - (phi_start + delta_phi)));
end
probe_i = src_i_ground;
probe_j = round(Nthe/2);
probe_record = zeros(num_steps, length(probe_ranges_km));

%% ============================================================
%% Visualization setup
%% ============================================================
viz_every          = 75; % Normally set to 75, reduced to 5 for testing
viz_noise_floor    = 1e-6;
viz_save_pngs      = false;
viz_png_dir        = './fdtd_frames';
viz_cross_range_km = dist_max/2E3; % Halfway through the distance
delta_phi_viz = viz_cross_range_km*1e3 / (R_Earth * cos(lat0_rad));
[~, viz_k_cross] = min(abs(phi_node - (phi_start + delta_phi_viz)));
viz_x = (phi_node - phi_start) * R_Earth * cos(lat0_rad) / 1e3;
viz_y = (the_node - the_node(probe_j)) * R_Earth / 1e3;
viz_z = (r_half - R_Earth) / 1e3;
[VZ, VY, VX] = ndgrid(viz_z, viz_y, viz_x);
slice_x = viz_x(viz_k_cross);
slice_y = 0;
slice_z = viz_z(1);
viz_fig = [];
viz_handles = struct();
if viz_save_pngs && ~exist(viz_png_dir, 'dir'); mkdir(viz_png_dir); end

report_every = max(1, round(num_steps/50));

%% ============================================================
%% Time loop
%% ============================================================
fprintf('\nStarting time stepping...\n');
tic_total = tic;
for n = 0:num_steps-1
    t_now = n * del_t;
    update_H();
    if sibc_enabled; apply_sibc_bottom(); end
    update_E();
    if test_source_enabled
        apply_test_source(t_now);
    else
        apply_source(t_now);    % the MTLE source
    end
    record_probes(n);
    if mod(n, viz_every) == 0
        update_visualization(n, t_now);
        drawnow limitrate;
    end
    if mod(n, report_every) == 0
        report_progress(n, tic_total);
    end
end
fprintf('\nTotal runtime: %.1f s\n', toc(tic_total));

%% ============================================================
%% Post-loop plotting
%% ============================================================
figure;
t_axis_ms = (0:num_steps-1) * del_t * 1e3;
for p = 1:length(probe_ranges_km)
    subplot(length(probe_ranges_km), 1, p);
    plot(t_axis_ms, 20*log10(probe_record(:, p))); grid on;
    % ylim([-40, max(20*log10(probe_record(:, p)))])
    ylabel(sprintf('E_r at %d km', probe_ranges_km(p)));
    if p == length(probe_ranges_km); xlabel('Time (ms)'); end
end
sgtitle('E_r probes vs time at various ranges from source');
%% Probe Filtered
fs = 1/del_t;                        % sampling frequency
fLo = 100e3; fHi = 300e3;                          % band-pass cutoff (below your 100-300 kHz band)
% [b, a] = butter(4, [fLo fHi], 'bandpass');   % 6th-order Butterworth bandpass
[z,p,k] = butter(6, [fLo fHi]/(fs/2), 'bandpass');
[sos,g] = zp2sos(z,p,k);
F
probe_filtered = zeros(size(probe_record));
for p = 1:length(probe_ranges_km)
    probe_filtered(:, p) = filtfilt(sos, g, probe_record(:, p));
end

% Compare raw vs filtered at one probe
figure;
t_axis_us = (0:num_steps-1) * del_t * 1e6;
for p = 1:length(probe_ranges_km)
    subplot(2,length(probe_ranges_km),p); 
    plot(t_axis_us, probe_record(:,p));
    title(sprintf('Raw E_r (static + radiation) at %d km', probe_ranges_km(p))); 
    ylabel('E_r'); grid on;
    subplot(2,length(probe_ranges_km),p+length(probe_ranges_km)); 
    plot(t_axis_us, probe_filtered(:,p));
    title('High-pass filtered (radiation only)'); 
    xlabel('Time (\mus)'); ylabel('E_r'); grid on;
end
sgtitle('E_r probes vs time at various ranges from source filtered');

%% Probe Spectral Content (PSD or Spectrogram)
figure;
t_axis_us = (0:num_steps-1) * del_t * 1e6;
for p = 1:length(probe_ranges_km)
    subplot(1,length(probe_ranges_km),p)
    spectrogram(probe_filtered(:,p));
end
sgtitle('E_r probes vs time at various ranges from source filtered');
validSaveName = matlab.lang.makeValidName("probe_data" +  string(datetime("now")) +  ".mat");
save(validSaveName, 'probe_record', 'del_t', 'probe_ranges_km', ...
    'sibc_enabled','cpml_ct_top','cpml_ct_bottom','cpml_ct_east', ...
    'cpml_ct_west','cpml_ct_north','cpml_ct_south')
%% Current Waveform
t_vec = (0:num_steps-1) * del_t;
i0_vec = arrayfun(@(tt) channel_base_current(tt), t_vec);

figure; plot(t_vec*1e6, i0_vec/1e3);
xlabel('Time (\mus)'); ylabel('i(0,t) [kA]');
title('Channel-base current waveform');

figure;
Fs = 1 / del_t;
N = length(i0_vec);
xdft = fft(i0_vec);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(i0_vec):Fs/2;
plot(freq/1E3,pow2db(psdx),'Linewidth',2)
xline(500, LineWidth=2) %places a line at 500
grid on
title("Periodogram of MTLE Current Waveform")
xlabel("Frequency (kHz)")
ylabel("Power/Frequency (dB/Hz)")

%% ============================================================
%% NESTED FUNCTIONS
%% ============================================================

function update_H()
    % --- H_r: curlE_r = (1/(r sinth)) [ d(E_phi sinth)/dth - dE_th/dphi ] ---
    Ephi_sin = E_phi .* sin_the_node_2;
    dEphiSin_dth = diff(Ephi_sin, 1, 2) / del_the;     % (Nr,Nthe-1,Nphi-1) theta-half
    dEth_dphi    = diff(E_th,     1, 3) / del_phi;     % (Nr,Nthe-1,Nphi-1) phi-half

    % theta-stretch on dEphiSin_dth (dim 2, half)
    [dEphiSin_dth, psi_Ephisin_th] = cpml_stretch(dEphiSin_dth, psi_Ephisin_th, 2, idx_the_half, bth_h, ath_h, kth_h);
    % phi-stretch on dEth_dphi (dim 3, half)
    [dEth_dphi, psi_Eth_phi] = cpml_stretch(dEth_dphi, psi_Eth_phi, 3, idx_phi_half, bph_h, aph_h, kph_h);

    curlE_r = (1 ./ (r_node_1 .* sin_the_half_2)) .* (dEphiSin_dth - dEth_dphi);
    H_r = H_r - CHrE * curlE_r;

    % --- H_th: curlE_th = (1/r) [ (1/sinth) dE_r/dphi - d(r E_phi)/dr ] ---
    dEr_dphi  = diff(E_r, 1, 3) / del_phi;             % (Nr-1,Nthe,Nphi-1) phi-half
    rEphi     = r_node_1 .* E_phi;
    drEphi_dr = diff(rEphi, 1, 1) / del_r;             % (Nr-1,Nthe,Nphi-1) r-half

    [dEr_dphi, psi_Er_phi]   = cpml_stretch(dEr_dphi, psi_Er_phi, 3, idx_phi_half, bph_h, aph_h, kph_h);
    [drEphi_dr, psi_rEphi_r] = cpml_stretch(drEphi_dr, psi_rEphi_r, 1, idx_r_half, br_h, ar_h, kr_h);

    curlE_th = (1 ./ r_half_1) .* ((1 ./ sin_the_node_2) .* dEr_dphi - drEphi_dr);
    H_th = H_th - CHthE * curlE_th;

    % --- H_phi: curlE_phi = (1/r) [ d(r E_th)/dr - dE_r/dth ] ---
    rEth     = r_node_1 .* E_th;
    drEth_dr = diff(rEth, 1, 1) / del_r;               % (Nr-1,Nthe-1,Nphi) r-half
    dEr_dth  = diff(E_r,  1, 2) / del_the;             % (Nr-1,Nthe-1,Nphi) theta-half

    [drEth_dr, psi_rEth_r] = cpml_stretch(drEth_dr, psi_rEth_r, 1, idx_r_half, br_h, ar_h, kr_h);
    [dEr_dth, psi_Er_th]   = cpml_stretch(dEr_dth, psi_Er_th, 2, idx_the_half, bth_h, ath_h, kth_h);

    curlE_phi = (1 ./ r_half_1) .* (drEth_dr - dEr_dth);
    H_phi = H_phi - CHphiE * curlE_phi;
end

function update_E()
    % --- E_r interior: curlH_r = (1/(r sinth)) [ d(H_phi sinth)/dth - dH_th/dphi ] ---
    Hphi_sin = H_phi .* sin_the_half_2;
    dHphiSin_dth = diff(Hphi_sin, 1, 2) / del_the;     % (Nr-1,Nthe-2,Nphi) theta-node(interior)
    dHth_dphi    = diff(H_th,     1, 3) / del_phi;     % (Nr-1,Nthe,Nphi-2)  phi-node(interior)

    % theta-stretch on dHphiSin_dth (dim 2, node). Its dim-2 length is Nthe-2,
    % indexed by interior nodes 2..Nthe-1 -> use idx_the_ndf directly.
    [dHphiSin_dth, psi_Hphisin_th] = cpml_stretch(dHphiSin_dth, psi_Hphisin_th, 2, idx_the_ndf, bth_n, ath_n, kth_n);
    % phi-stretch on dHth_dphi (dim 3, node)
    [dHth_dphi, psi_Hth_phi] = cpml_stretch(dHth_dphi, psi_Hth_phi, 3, idx_phi_ndf, bph_n, aph_n, kph_n);

    curlH_r_int = (1 ./ (r_half_1 .* reshape(sin_the_node(2:Nthe-1), [1 Nthe-2 1]))) ...
                  .* (dHphiSin_dth(:, :, 2:Nphi-1) - dHth_dphi(:, 2:Nthe-1, :));
    E_r(:, 2:Nthe-1, 2:Nphi-1) = CErE(:, 2:Nthe-1, 2:Nphi-1) .* E_r(:, 2:Nthe-1, 2:Nphi-1) ...
                              + CErH(:, 2:Nthe-1, 2:Nphi-1) .* curlH_r_int;

    % --- E_th interior: curlH_th = (1/r) [ (1/sinth) dH_r/dphi - d(r H_phi)/dr ] ---
    dHr_dphi  = diff(H_r, 1, 3) / del_phi;             % (Nr,Nthe-1,Nphi-2)  phi-node
    rHphi     = r_half_1 .* H_phi;
    drHphi_dr = diff(rHphi, 1, 1) / del_r;             % (Nr-2,Nthe-1,Nphi)  r-node

    [dHr_dphi, psi_Hr_phi]   = cpml_stretch(dHr_dphi, psi_Hr_phi, 3, idx_phi_ndf, bph_n, aph_n, kph_n);
    [drHphi_dr, psi_rHphi_r] = cpml_stretch(drHphi_dr, psi_rHphi_r, 1, idx_r_ndf, br_n, ar_n, kr_n);

    curlH_th_int = (1 ./ reshape(r_node(2:Nr-1), [Nr-2 1 1])) ...
                   .* ((1 ./ sin_the_half_2) .* dHr_dphi(2:Nr-1, :, :) ...
                       - drHphi_dr(:, :, 2:Nphi-1));
    E_th(2:Nr-1, :, 2:Nphi-1) = CEthE(2:Nr-1, :, 2:Nphi-1) .* E_th(2:Nr-1, :, 2:Nphi-1) ...
                              + CEthH(2:Nr-1, :, 2:Nphi-1) .* curlH_th_int;

    % --- E_phi interior: curlH_phi = (1/r) [ d(r H_th)/dr - dH_r/dth ] ---
    rHth     = r_half_1 .* H_th;
    drHth_dr = diff(rHth, 1, 1) / del_r;               % (Nr-2,Nthe,Nphi-1) r-node
    dHr_dth  = diff(H_r,  1, 2) / del_the;             % (Nr,Nthe-2,Nphi-1) theta-node

    [drHth_dr, psi_rHth_r] = cpml_stretch(drHth_dr, psi_rHth_r, 1, idx_r_ndf, br_n, ar_n, kr_n);
    [dHr_dth, psi_Hr_th]   = cpml_stretch(dHr_dth, psi_Hr_th, 2, idx_the_ndf, bth_n, ath_n, kth_n);

    curlH_phi_int = (1 ./ reshape(r_node(2:Nr-1), [Nr-2 1 1])) ...
                    .* (drHth_dr(:, 2:Nthe-1, :) - dHr_dth(2:Nr-1, :, :));
    E_phi(2:Nr-1, 2:Nthe-1, :) = CEphiE(2:Nr-1, 2:Nthe-1, :) .* E_phi(2:Nr-1, 2:Nthe-1, :) ...
                               + CEphiH(2:Nr-1, 2:Nthe-1, :) .* curlH_phi_int;
end

function apply_sibc_bottom()
    % Contour-path SIBC at the ocean surface (bottom row i=1).
    % Called AFTER update_H(): on entry, H_th(1,:,:) and H_phi(1,:,:) hold the
    % "explicit" (PEC-style) update, because the surface E rows E_th(1,:,:) and
    % E_phi(1,:,:) are held at zero (never written). We correct those rows in
    % place. Hth_prev1/Hphi_prev1 hold surface H^{n-1/2} (zero at rest).
    %
    % Both tangential H components share the same coefficients: the radii/dr at
    % i=1 are identical, and the sign difference between the E_th and E_phi SIBC
    % relations is cancelled by the opposite curl sign (see header derivation).

    Hphi_expl = H_phi(1, :, :);     % H_explicit^{n+1/2}, (1,Nthe-1,Nphi)
    Hth_expl  = H_th(1, :, :);      % H_explicit^{n+1/2}, (1,Nthe,Nphi-1)

    if strcmp(sibc_mode, 'mono')
        % -------- MONOCHROMATIC (single-freq, implicit) --------
        %   H^{n+1/2}(1) = ( H_expl(1) + c1 * H^{n-1/2}(1) ) / c2
        Hphi_new = (Hphi_expl + sibc_c1 * Hphi_prev1) / sibc_c2;
        Hth_new  = (Hth_expl  + sibc_c1 * Hth_prev1 ) / sibc_c2;

        if sibc_store_Esurf
            dHphi_dt = (Hphi_new - Hphi_prev1) / del_t;
            dHth_dt  = (Hth_new  - Hth_prev1 ) / del_t;
            Eth_surf_diag  = -sibc_Rs * ( 0.5*(Hphi_new + Hphi_prev1) + dHphi_dt/sibc_w0 );
            Ephi_surf_diag =  sibc_Rs * ( 0.5*(Hth_new  + Hth_prev1 ) + dHth_dt /sibc_w0 );
        end

    else
        % -------- CONVOLUTIONAL (dispersive, causal -> explicit) --------   <== NEW PATH
        % 1) Advance pole accumulators to step n using surface H^{n-1/2} (=prev):
        Psi_phi = sibc_Qb .* Psi_phi + sibc_Sb .* Hphi_prev1;   % (Np,Nthe-1,Nphi)
        Psi_th  = sibc_Qb .* Psi_th  + sibc_Sb .* Hth_prev1;    % (Np,Nthe,Nphi-1)

        % 2) Surface E from convolution history (real part of the pole sum):
        sumPsi_phi = real(sum(Psi_phi, 1));   % (1,Nthe-1,Nphi)
        sumPsi_th  = real(sum(Psi_th , 1));   % (1,Nthe,Nphi-1)

        % 3) Explicit surface-H correction (c2conv=1 unless Zinf~=0):
        Hphi_new = (Hphi_expl - sibc_g*(sibc_Zinf/2)*Hphi_prev1 - sibc_g*sumPsi_phi) / sibc_c2conv;
        Hth_new  = (Hth_expl  - sibc_g*(sibc_Zinf/2)*Hth_prev1  - sibc_g*sumPsi_th ) / sibc_c2conv;

        if sibc_store_Esurf
            Eth_surf_diag  = -( sibc_Zinf*0.5*(Hphi_new+Hphi_prev1) + sumPsi_phi );
            Ephi_surf_diag =  ( sibc_Zinf*0.5*(Hth_new +Hth_prev1 ) + sumPsi_th  );
        end
    end

    H_phi(1, :, :) = Hphi_new;
    H_th(1, :, :)  = Hth_new;

    % Store corrected H^{n+1/2} -> becomes H^{n-1/2} next step (both modes).
    Hphi_prev1 = Hphi_new;
    Hth_prev1  = Hth_new;
end

function apply_source(t_now)
    I_channel = mtle_current(channel_z(:), t_now);
    J_r_src   = I_channel ./ src_area(:);
    eps_src_line   = eps_cell_r(src_i_range(:));
    sigma_src_line = sigma_cell_r(src_i_range(:));
    src_injection  = -(del_t ./ (eps_src_line(:) + sigma_src_line(:).*del_t/2)) .* J_r_src;
    E_r(src_i_range, src_j, src_k) = E_r(src_i_range, src_j, src_k) + reshape(src_injection, [], 1, 1);
end

function apply_test_source(t_now)
    % Soft additive source into E_r at a single interior point.
    % Additive (+=) so it does not act as a hard scatterer.
    src_val = diff_gaussian(t_now, test_t0, test_tau, test_amp);
    E_r(test_i, test_j, test_k) = E_r(test_i, test_j, test_k) + src_val;
end

function record_probes(n)
    for p = 1:length(probe_ranges_km)
        probe_record(n+1, p) = E_r(probe_i, probe_j, probe_k_indices(p));
    end
end

function update_visualization(step_n, t_now)
    log_data = log10(max(abs(E_r), viz_noise_floor));
    cmin = log10(viz_noise_floor);
    cmax = max(log_data(:)); if cmax - cmin < 0.5; cmax = cmin + 0.5; end
    log_data_p = permute(log_data, [2 3 1]);
    VX_p = permute(VX, [2 3 1]); VY_p = permute(VY, [2 3 1]); VZ_p = permute(VZ, [2 3 1]);
    if isempty(viz_fig) || ~isvalid(viz_fig)
        viz_fig = figure('Name','E_r field (log scale)','Position',[100 100 1200 600],'Color','w');
        viz_handles.ax = axes('Parent', viz_fig);
        viz_handles.first_call = true;
    else
        figure(viz_fig);
        if isfield(viz_handles,'slice_handle')
            delete(viz_handles.slice_handle(isvalid(viz_handles.slice_handle)));
        end
    end
    viz_handles.slice_handle = slice(viz_handles.ax, VX_p, VY_p, VZ_p, log_data_p, slice_x, slice_y, slice_z);
    for kk = 1:length(viz_handles.slice_handle)
        set(viz_handles.slice_handle(kk), 'EdgeColor','none','FaceColor','interp');
    end
    if isfield(viz_handles,'first_call') && viz_handles.first_call
        colormap(viz_handles.ax, hot); hcb = colorbar(viz_handles.ax);
        ylabel(hcb,'log_{10}|E_r| (V/m)');
        xlabel(viz_handles.ax,'Range west (km)'); ylabel(viz_handles.ax,'N/S offset (km)');
        zlabel(viz_handles.ax,'Altitude (km)'); view(viz_handles.ax,[-30 25]); grid(viz_handles.ax,'on');
        viz_handles.first_call = false;
        viz_handles.xlim = [min(VX_p(:)) max(VX_p(:))];
        viz_handles.ylim = [min(VY_p(:)) max(VY_p(:))];
        viz_handles.zlim = [min(VZ_p(:)) max(VZ_p(:))];
    end
    xlim(viz_handles.ax, viz_handles.xlim); ylim(viz_handles.ax, viz_handles.ylim); zlim(viz_handles.ax, viz_handles.zlim);
    caxis(viz_handles.ax, [cmin cmax]);
    title(viz_handles.ax, sprintf('Step %d  t = %.2f \\mus  max|E_r|=%.2e V/m', step_n, t_now*1e6, max(abs(E_r(:)))));
    if viz_save_pngs
        exportgraphics(viz_fig, fullfile(viz_png_dir, sprintf('frame_%06d.png', step_n)), 'Resolution', 100);
    end
end

function report_progress(n_step, tic_handle)
    elapsed = toc(tic_handle); frac = (n_step+1)/num_steps; eta = elapsed*(1/frac - 1);
    fprintf('Step %d/%d (%.1f%%)  elapsed %.1f s  ETA %.1f s  max|E_r|=%.2e\n', ...
            n_step+1, num_steps, 100*frac, elapsed, eta, max(abs(E_r(:))));
end

end  % end main function

%% ============================================================
%% LOCAL (non-nested) HELPER FUNCTIONS
%% ============================================================

function [deriv, psi] = cpml_stretch(deriv, psi, dimn, idx, b_s, a_s, kap_s)
    % Apply CPML coordinate stretching to `deriv` along dimension `dimn`,
    % at index set `idx`. Coefficients b_s,a_s,kap_s are pre-reshaped to
    % broadcast along `dimn`. psi is updated and returned.
    % Reduces to a no-op where the profile is trivial (b=1,a=0,kappa=1).
    if isempty(idx); return; end
    if dimn == 1
        d = deriv(idx, :, :);
        psi = b_s .* psi + a_s .* d;
        deriv(idx, :, :) = d ./ kap_s + psi;
    elseif dimn == 2
        d = deriv(:, idx, :);
        psi = b_s .* psi + a_s .* d;
        deriv(:, idx, :) = d ./ kap_s + psi;
    else
        d = deriv(:, :, idx);
        psi = b_s .* psi + a_s .* d;
        deriv(:, :, idx) = d ./ kap_s + psi;
    end
end

function [b_int,a_int,kap_int, b_half,a_half,kap_half] = build_coord_profiles(...
            N, ct_low, ct_high, m, kappa_max, alpha_max, sigma_max, del_t, eps0)
    % Per-coordinate CPML profiles with ramps at low and/or high ends.
    % Integer positions 1..N, half positions 1..N-1 (half k between node k,k+1).
    % kappa=1, sigma=0, alpha=0 outside PML regions.
    kap_int=ones(1,N);   sig_int=zeros(1,N);   alf_int=zeros(1,N);
    kap_half=ones(1,N-1); sig_half=zeros(1,N-1); alf_half=zeros(1,N-1);

    if ct_high > 0
        for i = (N-ct_high+1):N
            d = (i-(N-ct_high))/ct_high;            % 1/ct .. 1 (wall)
            kap_int(i)=1+(kappa_max-1)*d^m; sig_int(i)=sigma_max*d^m; alf_int(i)=alpha_max*(1-d)^m;
        end
        for k = (N-ct_high):(N-1)
            d = (k+0.5-(N-ct_high))/ct_high; d=min(max(d,0),1);
            kap_half(k)=1+(kappa_max-1)*d^m; sig_half(k)=sigma_max*d^m; alf_half(k)=alpha_max*(1-d)^m;
        end
    end
    if ct_low > 0
        for i = 1:ct_low
            d = (ct_low-i+1)/ct_low;                % 1 (wall, i=1) .. 1/ct
            kap_int(i)=1+(kappa_max-1)*d^m; sig_int(i)=sigma_max*d^m; alf_int(i)=alpha_max*(1-d)^m;
        end
        for k = 1:ct_low
            d = (ct_low-k+0.5)/ct_low; d=min(max(d,0),1);
            kap_half(k)=1+(kappa_max-1)*d^m; sig_half(k)=sigma_max*d^m; alf_half(k)=alpha_max*(1-d)^m;
        end
    end

    b_int = exp(-(sig_int./kap_int + alf_int)*del_t/eps0);
    den = kap_int.*(sig_int+kap_int.*alf_int); a_int=zeros(1,N); mm=den>0;
    a_int(mm)=sig_int(mm)./den(mm).*(b_int(mm)-1);

    b_half = exp(-(sig_half./kap_half + alf_half)*del_t/eps0);
    den = kap_half.*(sig_half+kap_half.*alf_half); a_half=zeros(1,N-1); mm=den>0;
    a_half(mm)=sig_half(mm)./den(mm).*(b_half(mm)-1);
end

function i = channel_base_current(t)
    % User-defined MTLE channel-base current. Edit freely.
    % I01 = 10e3;  tau11 = 0.1e-6; tau12 = 5e-6;  n1 = 2;
    % I02 = 6e3;   tau21 = 2.0e-6; tau22 = 50e-6; n2 = 2;
    % eta1 = exp(-(tau11/tau12)*(n1*tau12/tau11)^(1/n1));
    % eta2 = exp(-(tau21/tau22)*(n2*tau21/tau22)^(1/n2));
    % t_safe = max(t, 0);
    % term1 = (I01/eta1) .* ((t_safe./tau11).^n1 ./ (1 + (t_safe./tau11).^n1)) .* exp(-t_safe./tau12);
    % term2 = (I02/eta2) .* ((t_safe./tau21).^n2 ./ (1 + (t_safe./tau21).^n2)) .* exp(-t_safe./tau22);
    % i = term1 + term2;
    I0 = 10E3;
    t_safe = max(t, 0);
    tauf = 50E-6;
    taur = 0.3E-6;
    term1 = exp(-t_safe./tauf);
    term2 = exp(-t_safe./taur);
    i = I0 .* (term1 - term2);
    i(t < 0) = 0;
end

function f = diff_gaussian(t, t0, tau, amp)
    % First derivative of a Gaussian: zero time-integral (charge-neutral),
    % spectrum peaked at 1/(2*pi*tau), no DC content.
    f = -amp * (t - t0) ./ (tau^2) .* exp(-((t - t0).^2) ./ (2*tau^2));
end