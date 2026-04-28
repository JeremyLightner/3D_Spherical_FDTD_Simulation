%% Spherical FDTD
clear all; clc;

c = 299792458;
mu0 = (4*pi)*10^-7;
eps0 = (c*sqrt(mu0))^-2;
num_dims = 3;
R_Earth = 6371 * 1E3; % m


f_max = NaN;
dist_max = NaN;
width_max = NaN;
height_max = NaN;

while isnan(f_max)
    f_max = input('Enter the maximum relevant frequency (Hz): ', "s");
    f_max = str2double(f_max);
    if isnan(f_max)
        disp('Nice try, thats not a number...')
    end
end

% GPT generated prompt for lat lon
while true
    str = input('Enter the start lat,lon (e.g., 34.5,-120.2): ', "s");
    
    % Convert string → numeric vector
    vals = str2num(str); %#ok<ST2NM>
    
    % Validation
    if isempty(vals) || numel(vals) ~= 2 || any(isnan(vals))
        disp('Invalid input. Enter two numeric values like: 34.5,-120.2')
        continue
    end
    
    lat = vals(1);
    lon = vals(2);
    
    % Optional: enforce valid geographic bounds
    if lat < -90 || lat > 90 || lon < -180 || lon > 180
        disp('Lat must be [-90,90], Lon must be [-180,180]')
        continue
    end
    
    start_point = [lat, lon];
    break
end

% Look angle will always be westward (we want to look eastward for over
% ocean propagation)

while isnan(dist_max)
    dist_max = input('Enter the maximum distance in range (m): ', "s");
    dist_max = str2double(dist_max);
    if isnan(dist_max)
        disp('Nice try, thats not a number...')
    end
end

while isnan(width_max)
    width_max = input('Enter the maximum width of simulation (m): ', "s");
    width_max = str2double(width_max);
    if isnan(width_max)
        disp('Nice try, thats not a number...')
    end
end

while isnan(height_max)
    height_max = input('Enter the maximum height of simulation (m): ', "s");
    height_max = str2double(height_max);
    if isnan(height_max)
        disp('Nice try, thats not a number...')
    end
end

% convert dist_max width and height to a 3d volume on the surface of the
% earth, including a little bit of space below the ground
% The distance is the movement in azimuth in km, the width in km is the size along
% the latitude lines and the height is the height in km.

% Additional input: depth below surface
below_ground = NaN;
while isnan(below_ground)
    below_ground = input('Enter the depth below the surface (m): ', "s");
    below_ground = str2double(below_ground);
    if isnan(below_ground) || below_ground < 0
        below_ground = NaN;
        disp('Nice try, enter a nonnegative number...')
    end
end

% Build spherical patch
% Conventions:
%   latitude  = geographic latitude, + north
%   longitude = geographic longitude, + east
%
%   theta = polar angle measured from north pole
%         = 0 at north pole, pi/2 at equator, pi at south pole
%
%   phi   = azimuth measured WESTWARD from the prime meridian
%         = 0 at lon = 0
%         = positive when moving west
%
% Start point is the CENTER of the width and the START of the westward range.

lat0 = start_point(1);   % deg
lon0 = start_point(2);   % deg

% Radial bounds
r_inner = R_Earth - below_ground;
r_outer = R_Earth + height_max;

% Width is north/south centered on the start point
dlat_half_deg = rad2deg((width_max/2) / R_Earth);

lat_min = lat0 - dlat_half_deg;
lat_max = lat0 + dlat_half_deg;

% Keep within valid geographic limits
lat_min = max(lat_min, -90);
lat_max = min(lat_max,  90);

% Westward range along constant latitude line through lat0
lat0_rad = deg2rad(lat0);
if abs(cos(lat0_rad)) < 1e-12
    error('Start latitude is too close to a pole for westward longitude stepping.')
end

dlon_range_deg = rad2deg(dist_max / (R_Earth * cos(lat0_rad)));

% Geographic longitude decreases when moving west
lon_east = lon0;                  % starting edge
lon_west = lon0 - dlon_range_deg; % ending edge after moving west

% Convert geographic lat/lon to your spherical coordinates
% theta starts at north pole
theta_min = deg2rad(90 - lat_max);
theta_max = deg2rad(90 - lat_min);

% phi is positive westward from Greenwich:
%   phi = -lon_eastward
% phi_start = wrapTo2Pi_local(deg2rad(-lon_east));
% phi_end   = wrapTo2Pi_local(deg2rad(-lon_west));
% Unwrapped west-positive phi values
phi_start = deg2rad(-lon_east);
phi_end   = deg2rad(-lon_west);

% Force phi to increase westward over the small intended patch
while phi_end < phi_start
    phi_end = phi_end + 2*pi;
end

fprintf('\nPatch definition:\n');
fprintf('  Latitude bounds      : [%g, %g] deg\n', lat_min, lat_max);
fprintf('  Longitude bounds     : from %g deg to %g deg (westward)\n', lon_east, lon_west);
fprintf('  Theta bounds         : [%g, %g] rad\n', theta_min, theta_max);
fprintf('  Phi start/end        : [%g, %g] rad\n', phi_start, phi_end);
fprintf('  Radial bounds        : [%g, %g] km\n', r_inner, r_outer);

% Plot highlighted patch on Earth surface
plot_earth_volume(theta_min, theta_max, phi_start, phi_end, r_inner, r_outer, R_Earth);


%% Define Grid Structure

% This does not account for the smallest value of epsr: we will address
% this

lam_min    = c/f_max;
v_p        = lam_min*f_max;

%Grid Dispersion Error req ds << lambda
ds = lam_min/100; 

theta_temp = theta_min:0.0001:theta_max; %radians;

del_r      = ds;
del_the    = ds / (r_outer);
del_phi    = ds / (r_outer * max(sin(theta_temp)));

%This will be hardcoded once I know the space

Nr         = ceil((r_outer-r_inner)/del_r);
Nthe       = ceil(abs(theta_max-theta_min)/del_the);
Nphi       = ceil(abs(phi_start-phi_end)/del_phi);

E_r = rand(Nr-1, Nthe, Nphi);
E_the = rand(Nr, Nthe-1, Nphi);
E_phi = rand(Nr, Nthe, Nphi-1);

H_r = rand(Nr, Nthe-1, Nphi-1);
H_the = rand(Nr-1, Nthe, Nphi-1);
H_phi = rand(Nr-1, Nthe-1, Nphi);

%% Establish time step using CFL Criterion

%correct to spherical coord version
del_t = (1/c)/sqrt((1/(del_r)^2) + (1/(r_inner*del_the)^2) + (1/(r_inner*sin(theta_min)*del_phi)^2));

%% Create Spatial Property Matrices
% Material Properties are defined at the 
% center of cell volumes in spherical grid
eps_cell   = rand(Nr+2,Nthe+2,Nphi+2); % sweeps from 0:Nr+1
sigma_cell = rand(Nr+2,Nthe+2,Nphi+2); % sweeps from 0:Nr+1
% Use these to interpolate to the larger space
eps_r      = 0.25 * (eps_cell(2:Nr+2,2:Nthe+2,2:Nphi+2)...
                + eps_cell(2:Nr+2,1:Nthe+1,2:Nphi+2)...
                + eps_cell(2:Nr+2,2:Nthe+2,1:Nphi+1)...
                + eps_cell(2:Nr+2,1:Nthe+1,1:Nphi+1));

eps_the    = 0.25 * (eps_cell(2:Nr+2,2:Nthe+2,2:Nphi+2)...
                + eps_cell(1:Nr+1,2:Nthe+2,2:Nphi+2)...
                + eps_cell(2:Nr+2,2:Nthe+2,1:Nphi+1)...
                + eps_cell(1:Nr+1,2:Nthe+2,1:Nphi+1));

eps_phi    = 0.25 * (eps_cell(2:Nr+2,2:Nthe+2,2:Nphi+2)...
                + eps_cell(1:Nr+1,2:Nthe+2,2:Nphi+2)...
                + eps_cell(2:Nr+2,1:Nthe+1,2:Nphi+2)...
                + eps_cell(1:Nr+1,1:Nthe+1,2:Nphi+2));

sigma_r    = 0.25 * (sigma_cell(2:Nr+2,2:Nthe+2,2:Nphi+2)...
                 + sigma_cell(2:Nr+2,1:Nthe+1,2:Nphi+2)...
                 + sigma_cell(2:Nr+2,2:Nthe+2,1:Nphi+1)...
                 + sigma_cell(2:Nr+2,1:Nthe+1,1:Nphi+1));

sigma_the  = 0.25 * (sigma_cell(2:Nr+2,2:Nthe+2,2:Nphi+2)...
                 + sigma_cell(1:Nr+1,2:Nthe+2,2:Nphi+2)...
                 + sigma_cell(2:Nr+2,2:Nthe+2,1:Nphi+1)...
                 + sigma_cell(1:Nr+1,2:Nthe+2,1:Nphi+1));

sigma_phi  = 0.25 * (sigma_cell(2:Nr+2,2:Nthe+2,2:Nphi+2)...
                 + sigma_cell(1:Nr+1,2:Nthe+2,2:Nphi+2)...
                 + sigma_cell(2:Nr+2,1:Nthe+1,2:Nphi+2)...
                 + sigma_cell(1:Nr+1,1:Nthe+1,2:Nphi+2));

mu_cell    = rand(Nr+2,Nthe+2,Nphi+2); % sweeps from 0:Nr+1 (instead defined by 1:Nr+2)

mu_r       = 0.5 * (mu_cell(2:Nr+2,2:Nthe+2,2:Nphi+2))...
                + mu_cell(1:Nr+1,2:Nthe+2,2:Nphi+2);

mu_the       = 0.5 * (mu_cell(2:Nr+2,2:Nthe+2,2:Nphi+2))...
                + mu_cell(2:Nr+2,1:Nthe+1,2:Nphi+2);

mu_phi       = 0.5 * (mu_cell(2:Nr+2,2:Nthe+2,2:Nphi+2))...
                + mu_cell(2:Nr+2,2:Nthe+2,1:Nphi+1);

%Precompute Constants to save resources, takes more ram
CErE       = (1 - sigma_r*del_t./(2*eps_r)) ./...
                (1 + sigma_r*del_t./(2*eps_r));
CEtheE     = (1 - sigma_the*del_t./(2*eps_the)) ./...
                (1 + sigma_the*del_t./(2*eps_the));
CEphiE     = (1 - sigma_phi*del_t./(2*eps_phi)) ./...
                (1 + sigma_phi*del_t./(2*eps_phi));
CErH       = (del_t./ (eps_r + sigma_r*(del_t/2)));
CEtheH     = (del_t./ (eps_the + sigma_the*(del_t/2)));
CEphiH     = (del_t./ (eps_phi + sigma_phi*(del_t/2)));

%Need to precompute the H component multiplier CHr, CHthe, and CHphi

clear eps_r eps_the eps_phi sigma_r sigma_the sigma_phi ...
    eps_cell sigma_cell mu_cell


%% Create Fields

for t = 0:num_steps
    tic
    % Update H Field : UPDATE THESE TO BE VECTORIZED

    Hr = Hr + CHr*(  )

    % Update E Field: UPDATE THESE TO BE VECTORIZED

    % Inject Source
    
    % Probe the fields
    toc
end

%%
function plot_earth_volume(theta_min, theta_max, phi_start, phi_end, r_inner, r_outer, R_Earth)
    % Plots Earth plus the full 3D spherical patch volume.
    %
    % Conventions:
    %   theta = 0 at north pole
    %   phi   increases westward
    %
    % Cartesian mapping:
    %   x = r sin(theta) cos(phi)
    %   y = -r sin(theta) sin(phi)
    %   z = r cos(theta)

    figure;
    hold on;
    axis equal;
    grid on;
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Earth with Highlighted 3D Simulation Volume');

    % -----------------------------
    % Draw Earth surface
    % -----------------------------
    nSphere = 100;
    [Xs, Ys, Zs] = sphere(nSphere);
    surf(R_Earth*Xs, R_Earth*Ys, R_Earth*Zs, ...
        'FaceColor', [0.75 0.8 0.9], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.25);

    % Resolution for patch surfaces
    nTheta = 40;
    nPhi   = 50;
    nR     = 20;

    theta_vec = linspace(theta_min, theta_max, nTheta);
    phi_vec   = linspace(phi_start, phi_end, nPhi);
    r_vec     = linspace(r_inner, r_outer, nR);

    % =========================================================
    % 1) Bottom surface
    % =========================================================
    [Phi_b, Theta_b] = meshgrid(phi_vec, theta_vec);
    R_b = r_inner * ones(size(Phi_b));
    [Xb, Yb, Zb] = sph_west_to_cart(R_b, Theta_b, Phi_b);

    surf(Xb, Yb, Zb, ...
        'FaceColor', [1 0.2 0.2], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.35);

    % =========================================================
    % 2) Top surface
    % =========================================================
    [Phi_t, Theta_t] = meshgrid(phi_vec, theta_vec);
    R_t = r_outer * ones(size(Phi_t));
    [Xt, Yt, Zt] = sph_west_to_cart(R_t, Theta_t, Phi_t);

    surf(Xt, Yt, Zt, ...
        'FaceColor', [1 0.2 0.2], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.65);

    % =========================================================
    % 3) Theta = theta_min side wall
    % =========================================================
    [Phi_s1, R_s1] = meshgrid(phi_vec, r_vec);
    Theta_s1 = theta_min * ones(size(Phi_s1));
    [Xs1, Ys1, Zs1] = sph_west_to_cart(R_s1, Theta_s1, Phi_s1);

    surf(Xs1, Ys1, Zs1, ...
        'FaceColor', [0.85 0.1 0.1], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.45);

    % =========================================================
    % 4) Theta = theta_max side wall
    % =========================================================
    [Phi_s2, R_s2] = meshgrid(phi_vec, r_vec);
    Theta_s2 = theta_max * ones(size(Phi_s2));
    [Xs2, Ys2, Zs2] = sph_west_to_cart(R_s2, Theta_s2, Phi_s2);

    surf(Xs2, Ys2, Zs2, ...
        'FaceColor', [0.85 0.1 0.1], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.45);

    % =========================================================
    % 5) Phi = phi_start side wall
    % =========================================================
    [Theta_s3, R_s3] = meshgrid(theta_vec, r_vec);
    Phi_s3 = phi_start * ones(size(Theta_s3));
    [Xs3, Ys3, Zs3] = sph_west_to_cart(R_s3, Theta_s3, Phi_s3);

    surf(Xs3, Ys3, Zs3, ...
        'FaceColor', [0.85 0.1 0.1], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.45);

    % =========================================================
    % 6) Phi = phi_end side wall
    % =========================================================
    [Theta_s4, R_s4] = meshgrid(theta_vec, r_vec);
    Phi_s4 = phi_end * ones(size(Theta_s4));
    [Xs4, Ys4, Zs4] = sph_west_to_cart(R_s4, Theta_s4, Phi_s4);

    surf(Xs4, Ys4, Zs4, ...
        'FaceColor', [0.85 0.1 0.1], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.45);

    % -----------------------------
    % Draw outer edges for clarity
    % -----------------------------
    plot_edge(theta_min*ones(size(phi_vec)), phi_vec, r_outer, 'k', 2);
    plot_edge(theta_max*ones(size(phi_vec)), phi_vec, r_outer, 'k', 2);
    plot_edge(theta_vec, phi_start*ones(size(theta_vec)), r_outer, 'k', 2);
    plot_edge(theta_vec, phi_end*ones(size(theta_vec)), r_outer, 'k', 2);

    % Optional inner edges
    plot_edge(theta_min*ones(size(phi_vec)), phi_vec, r_inner, 'k--', 1);
    plot_edge(theta_max*ones(size(phi_vec)), phi_vec, r_inner, 'k--', 1);
    plot_edge(theta_vec, phi_start*ones(size(theta_vec)), r_inner, 'k--', 1);
    plot_edge(theta_vec, phi_end*ones(size(theta_vec)), r_inner, 'k--', 1);

    view(35, 25);
    camlight headlight;
    lighting gouraud;
    hold off;
end

function [x, y, z] = sph_west_to_cart(r, theta, phi)
    x = r .* sin(theta) .* cos(phi);
    y = -r .* sin(theta) .* sin(phi);
    z = r .* cos(theta);
end

function plot_edge(theta_vals, phi_vals, r, style, lineWidth)
    [x, y, z] = sph_west_to_cart(r, theta_vals, phi_vals);
    plot3(x, y, z, style, 'LineWidth', lineWidth);
end