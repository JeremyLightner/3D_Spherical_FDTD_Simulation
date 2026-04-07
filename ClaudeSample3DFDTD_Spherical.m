%% =========================================================================
%  3D Spherical FDTD Solver — Lat / Lon / Altitude Grid
%  Yee Scheme on a spherical shell domain
% ==========================================================================
clear; clc; close all;

%% =========================================================================
%  CONSTANTS
% ==========================================================================
c0   = 299792458;               % Speed of light (m/s)
mu0  = 4*pi*1e-7;         % Permeability of free space (H/m)
eps0 = 1/(mu0*c0^2);      % Permittivity of free space (F/m)

%% =========================================================================
%  GRID PARAMETERS
% ==========================================================================
% Angular extents (degrees)
lat_min =  -1;  lat_max =  1;   d_lat = 0.01;   % deg
lon_min = -30;  lon_max = 0;   d_lon = 0.5;   % deg

% Radial extents (meters above Earth surface)
R_earth = 6.371e6;        % Earth radius (m)
alt_min = 0;              % Surface
alt_max = 1e3;          % 1 km (e.g. lower ionosphere)
d_alt   = 100;            % 100 m radial step

% Build coordinate vectors
lat_vec = lat_min : d_lat : lat_max;       % (nlat x 1)  degrees
lon_vec = lon_min : d_lon : lon_max;       % (nlon x 1)  degrees
alt_vec = alt_min : d_alt : alt_max;       % (nalt x 1)  meters above surface
r_vec   = R_earth + alt_vec;               % (nalt x 1)  geocentric radius (m)

nlat = numel(lat_vec);
nlon = numel(lon_vec);
nalt = numel(alt_vec);

fprintf('Grid size: %d (lat) x %d (lon) x %d (alt) = %.2f M cells\n', ...
        nlat, nlon, nalt, nlat*nlon*nalt/1e6);

%% =========================================================================
%  METRIC COEFFICIENTS  (arc-length spacings vary with lat & alt)
% ==========================================================================
lat_rad = deg2rad(lat_vec(:));             % (nlat x 1)

% Meridional arc-length step (same at all lats & alts)
ds_lat = d_alt;                            % placeholder; real: r * d_lat_rad
d_lat_rad = deg2rad(d_lat);
d_lon_rad = deg2rad(d_lon);

% Zonal arc-length step:  ds_lon(lat, alt) = r(alt) * cos(lat) * d_lon_rad
%   → (nlat x nalt) matrix
cos_lat  = cos(lat_rad);                   % (nlat x 1)
[COS_LAT, R_MAT] = meshgrid(r_vec, cos_lat);   % both (nlat x nalt)
ds_lon   = COS_LAT .* R_MAT * d_lon_rad;       % (nlat x nalt)  zonal spacing

% Meridional arc-length:  ds_lat(alt) = r(alt) * d_lat_rad
ds_theta = r_vec(:)' * d_lat_rad;              % (1 x nalt)

% Radial step (uniform)
ds_r = d_alt;                              % scalar (m)

fprintf('Minimum zonal spacing (at poles): %.1f m\n', min(ds_lon(:)));
fprintf('Maximum zonal spacing (equator, top): %.1f m\n', max(ds_lon(:)));

%% =========================================================================
%  TIME STEPPING  (CFL condition)
% ==========================================================================
ds_min  = min([min(ds_lon(:)), min(ds_theta(:)), ds_r]);
courant = 0.9;
dt      = courant * ds_min / (c0 * sqrt(3));   % 3-D CFL
t_end   = 1e-3;                                % simulation duration (s)
nt      = ceil(t_end / dt);
t_vec   = (0:nt-1) * dt;

fprintf('dt = %.4e s,  Nt = %d steps\n', dt, nt);

%% =========================================================================
%  FIELD ARRAYS  (Yee staggered grid)
%  Dimensions: (nlat, nlon, nalt)
%  E-field components at cell edges; H-field at cell faces
% ==========================================================================
% Electric fields
E_r   = zeros(nlat, nlon, nalt);   % radial
E_lat = zeros(nlat, nlon, nalt);   % meridional (theta)
E_lon = zeros(nlat, nlon, nalt);   % zonal      (phi)

% Magnetic fields
H_r   = zeros(nlat, nlon, nalt);
H_lat = zeros(nlat, nlon, nalt);
H_lon = zeros(nlat, nlon, nalt);

%% =========================================================================
%  MATERIAL PARAMETERS
% ==========================================================================
eps_r = ones(nlat, nlon, nalt);    % Relative permittivity
mu_r  = ones(nlat, nlon, nalt);    % Relative permeability
sigma = zeros(nlat, nlon, nalt);   % Electric conductivity (S/m)

% --- TODO: Load realistic material maps ---
% e.g., ionospheric electron density, land/sea surface conductivity
% sigma = load_conductivity_map(lat_vec, lon_vec, alt_vec);

%% =========================================================================
%  PRECOMPUTE UPDATE COEFFICIENTS
% ==========================================================================
eps_abs = eps_r * eps0;
mu_abs  = mu_r  * mu0;

% E-field coefficients (accounts for ohmic loss)
Ca = (1 - sigma.*dt./(2*eps_abs)) ./ (1 + sigma.*dt./(2*eps_abs));
Cb =  (dt./eps_abs)               ./ (1 + sigma.*dt./(2*eps_abs));

% H-field coefficients (lossless magnetic medium)
Da = ones(nlat, nlon, nalt);
Db = dt ./ mu_abs;

%% =========================================================================
%  ABSORBING BOUNDARY CONDITIONS  (PML)
% ==========================================================================
pml_thickness = 10;   % cells

% --- TODO: Build PML conductivity profiles ---
% Polynomial grading inside PML layers at lat/lon/alt boundaries
% sigma_pml = build_pml(nlat, nlon, nalt, pml_thickness, ds_min, eps0);

%% =========================================================================
%  SOURCE PARAMETERS
% ==========================================================================
f0     = 1e3;                   % Source frequency (Hz) — ELF example
t0     = 3.0 / f0;              % Pulse delay (3 cycles)
spread = 1.0 / f0;              % Gaussian width

% Source location (grid indices)
src_lat = round(nlat/2);
src_lon = round(nlon/2);
src_alt = 1;                    % Near surface

%% =========================================================================
%  OBSERVATION / OUTPUT
% ==========================================================================
obs_lat = src_lat + 10;
obs_lon = src_lon + 10;
obs_alt = 1;

E_r_ts  = zeros(1, nt);        % Time series at observer

%% =========================================================================
%  MAIN TIME-MARCH LOOP
% ==========================================================================
fprintf('Starting FDTD time march...\n');
tic;

for n = 1:nt

    %% --- 1. Update H fields (curl of E) ---
    % dE_lon/d_theta - dE_lat/dr  → H_r
    % dE_r/dr        - dE_lon/dr  → H_lat  (placeholder finite differences)
    % dE_lat/dr      - dE_r/d_lat → H_lon

    % Radial H component
    dElat_dlat = (circshift(E_lat, -1, 1) - E_lat) / d_lat_rad;  % approx
    dElon_dlon = (circshift(E_lon, -1, 2) - E_lon) / d_lon_rad;
    H_r = Da .* H_r - Db .* (dElat_dlat - dElon_dlon);           % TODO: add r*cos(lat) metric

    % Meridional H component
    dEr_dalt  = (circshift(E_r,  -1, 3) - E_r ) / ds_r;
    dElon_dr  = (circshift(E_lon,-1, 3) - E_lon) / ds_r;
    H_lat = Da .* H_lat - Db .* (dEr_dalt - dElon_dr);           % TODO: full curl

    % Zonal H component
    dElat_dr  = (circshift(E_lat,-1, 3) - E_lat) / ds_r;
    dEr_dlat  = (circshift(E_r,  -1, 1) - E_r ) / d_lat_rad;
    H_lon = Da .* H_lon - Db .* (dElat_dr - dEr_dlat);           % TODO: full curl

    %% --- 2. Apply H-field boundary conditions ---
    % TODO: PML damping, periodic lon wrap (already via circshift)

    %% --- 3. Update E fields (curl of H) ---
    dHlon_dalt = (H_lon - circshift(H_lon, 1, 3)) / ds_r;
    dHlat_dalt = (H_lat - circshift(H_lat, 1, 3)) / ds_r;
    dHr_dlat   = (H_r   - circshift(H_r,   1, 1)) / d_lat_rad;
    dHr_dlon   = (H_r   - circshift(H_r,   1, 2)) / d_lon_rad;
    dHlon_dlon = (H_lon - circshift(H_lon, 1, 2)) / d_lon_rad;
    dHlat_dlat = (H_lat - circshift(H_lat, 1, 1)) / d_lat_rad;

    E_r   = Ca .* E_r   + Cb .* (dHlon_dalt - dHlat_dalt);  % TODO: metric factors
    E_lat = Ca .* E_lat + Cb .* (dHr_dlat   - dHlon_dlon);
    E_lon = Ca .* E_lon + Cb .* (dHlat_dlat - dHr_dlon  );

    %% --- 4. Inject source ---
    src_val = gaussian_pulse(t_vec(n), t0, spread);
    E_r(src_lat, src_lon, src_alt) = E_r(src_lat, src_lon, src_alt) + src_val;

    %% --- 5. Apply E-field boundary conditions ---
    % TODO: PML damping on lat/alt edges; lon is naturally periodic

    %% --- 6. Record observables ---
    E_r_ts(n) = E_r(obs_lat, obs_lon, obs_alt);

    %% --- 7. (Optional) Snapshot output ---
    if mod(n, 100) == 0
        elapsed = toc;
        fprintf('  Step %5d / %5d  |  sim time = %.4e s  |  wall = %.1f s\n', ...
                n, nt, t_vec(n), elapsed);
    end

end

fprintf('FDTD complete. Total wall time: %.2f s\n', toc);

%% =========================================================================
%  POST-PROCESSING
% ==========================================================================
figure;
plot(t_vec*1e3, E_r_ts);
xlabel('Time (ms)'); ylabel('E_r (V/m)');
title('Time Series at Observer Point');
grid on;

%% =========================================================================
%  LOCAL FUNCTIONS
% ==========================================================================
function val = gaussian_pulse(t, t0, spread)
    val = exp(-((t - t0).^2) / (2 * spread^2));
end