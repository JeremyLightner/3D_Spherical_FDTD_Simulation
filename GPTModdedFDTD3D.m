%% 3D FDTD - Jeremiah Lightner 2026
% Ground conductivities derived from ITU-R P.832-2
% Rewritten to use a Yee grid and a first-order 3D ABC.
%
% Field array sizes:
%   Ex : (Nx-1, Ny,   Nz)
%   Ey : (Nx,   Ny-1, Nz)
%   Ez : (Nx,   Ny,   Nz-1)
%   Hx : (Nx,   Ny-1, Nz-1)
%   Hy : (Nx-1, Ny,   Nz-1)
%   Hz : (Nx-1, Ny-1, Nz)

clear;
clc;

%% Physical constants
c    = 299792458;
mu0  = 4*pi*1e-7;
eps0 = 1/(mu0*c^2);

num_dims = 3;

%% User input
f_max = NaN;
while isnan(f_max) || f_max <= 0
  f_max = str2double(input('Enter the maximum relevant frequency (Hz): ', 's'));
  if isnan(f_max) || f_max <= 0
    disp('Please enter a positive number.');
  end
end

%% Define grid size here
% These Nx, Ny, Nz are the main Yee-grid node counts used by Ex/Ey/Ez layouts.
Nx = ceil(401/2);
Ny = ceil(401/2);
Nz = ceil(401/2);

%% Define cell-centered material grids here
% These are cell/property grids that will be sampled onto Yee locations.
if num_dims == 1
  epsr   = ones(1, Nx);
  mur    = ones(1, Nx);
  sigmam = zeros(1, Nx);
  sigma  = zeros(1, Nx);

elseif num_dims == 2
  epsr   = ones(Nx, Ny);
  mur    = ones(Nx, Ny);
  sigmam = zeros(Nx, Ny);
  sigma  = zeros(Nx, Ny);

elseif num_dims == 3
  epsr   = ones(Nx, Ny, Nz);
  mur    = ones(Nx, Ny, Nz);
  sigmam = zeros(Nx, Ny, Nz);
  sigma  = zeros(Nx, Ny, Nz);

  % Example PEC-like block:
  % sigma(30:33, 30:33, 33:66) = 1e16;
end

%% Smallest wavelength and grid spacing
lam0_min = c / f_max;
lam_min  = lam0_min / sqrt(max(epsr(:)));
v_p      = lam_min * f_max;

% Use ~29 cells per minimum wavelength [ref pg 82 of taflove book]
ds = lam_min / 20;

if num_dims == 1
  dx = ds;
  max_dim_size = (Nx-1) * dx;

elseif num_dims == 2
  dx = ds;
  dy = ds;
  max_dim_size = max([(Nx-1)*dx, (Ny-1)*dy]);

elseif num_dims == 3
  dx = ds;
  dy = ds;
  dz = ds;
  max_dim_size = max([(Nx-1)*dx, (Ny-1)*dy, (Nz-1)*dz]);
end

disp("Max physical length of simulation: " + max_dim_size + " m")
disp("Max one-way propagation time at slowest phase velocity: " + max_dim_size/v_p + " s")
disp("Approximate round-trip propagation time: " + 2*max_dim_size/v_p + " s")

%% CFL timestep
if num_dims == 1
  dt = dx / v_p;

elseif num_dims == 2
  if dx ~= dy
    dt = 1 / (v_p * sqrt(1/dx^2 + 1/dy^2));
  else
    dt = ds / (v_p * sqrt(num_dims));
  end

elseif num_dims == 3
  if dx ~= dy || dy ~= dz || dx ~= dz
    dt = 1 / (v_p * sqrt(1/dx^2 + 1/dy^2 + 1/dz^2));
  else
    dt = ds / (v_p * sqrt(num_dims));
  end
end

disp("Time step dt: " + dt + " s")

%% Maximum suggested simulation time
max_steps = ceil((2 * max_dim_size / v_p) / dt * 100);  % generous upper limit
t_max = NaN;

while isnan(t_max) || t_max <= 0
  t_max = str2double(input('Enter the quantity of time to simulate (s): ', 's'));

  if isnan(t_max) || t_max <= 0
    disp('Please enter a positive number.');
    t_max = NaN;
  elseif t_max >= max_steps * dt
    disp('WOAH... that is way too long')
    t_max = NaN;
  end
end

num_steps = ceil(t_max / dt);
disp("Number of time steps: " + num_steps)

%% Create Yee-grid field arrays
Ex = zeros(Nx-1, Ny,   Nz);
Ey = zeros(Nx,   Ny-1, Nz);
Ez = zeros(Nx,   Ny,   Nz-1);

Hx = zeros(Nx,   Ny-1, Nz-1);
Hy = zeros(Nx-1, Ny,   Nz-1);
Hz = zeros(Nx-1, Ny-1, Nz);

%% Material properties sampled to Yee locations
% Simple sampling from the cell-centered arrays.
% For more accuracy near interfaces, replace with averaging.

% H-field material arrays
mur_Hx    = mur(:,        1:Ny-1, 1:Nz-1);
sigmam_Hx = sigmam(:,     1:Ny-1, 1:Nz-1);

mur_Hy    = mur(1:Nx-1,   :,      1:Nz-1);
sigmam_Hy = sigmam(1:Nx-1,:,      1:Nz-1);

mur_Hz    = mur(1:Nx-1,   1:Ny-1, :);
sigmam_Hz = sigmam(1:Nx-1,1:Ny-1, :);

% E-field material arrays
epsr_Ex   = epsr(1:Nx-1, :,      :);
sigma_Ex  = sigma(1:Nx-1,:,      :);

epsr_Ey   = epsr(:,      1:Ny-1, :);
sigma_Ey  = sigma(:,     1:Ny-1, :);

epsr_Ez   = epsr(:,      :,      1:Nz-1);
sigma_Ez  = sigma(:,     :,      1:Nz-1);

%% Precompute update coefficients
C_Hx_H = (1 - (sigmam_Hx .* dt ./ (2 * mu0 .* mur_Hx))) ./ ...
         (1 + (sigmam_Hx .* dt ./ (2 * mu0 .* mur_Hx)));
C_Hx_E = dt ./ ((1 + (sigmam_Hx .* dt ./ (2 * mu0 .* mur_Hx))) .* mu0 .* mur_Hx .* ds);

C_Hy_H = (1 - (sigmam_Hy .* dt ./ (2 * mu0 .* mur_Hy))) ./ ...
         (1 + (sigmam_Hy .* dt ./ (2 * mu0 .* mur_Hy)));
C_Hy_E = dt ./ ((1 + (sigmam_Hy .* dt ./ (2 * mu0 .* mur_Hy))) .* mu0 .* mur_Hy .* ds);

C_Hz_H = (1 - (sigmam_Hz .* dt ./ (2 * mu0 .* mur_Hz))) ./ ...
         (1 + (sigmam_Hz .* dt ./ (2 * mu0 .* mur_Hz)));
C_Hz_E = dt ./ ((1 + (sigmam_Hz .* dt ./ (2 * mu0 .* mur_Hz))) .* mu0 .* mur_Hz .* ds);

C_Ex_E = (1 - (sigma_Ex .* dt ./ (2 * eps0 .* epsr_Ex))) ./ ...
         (1 + (sigma_Ex .* dt ./ (2 * eps0 .* epsr_Ex)));
C_Ex_H = dt ./ ((1 + (sigma_Ex .* dt ./ (2 * eps0 .* epsr_Ex))) .* eps0 .* epsr_Ex .* ds);

C_Ey_E = (1 - (sigma_Ey .* dt ./ (2 * eps0 .* epsr_Ey))) ./ ...
         (1 + (sigma_Ey .* dt ./ (2 * eps0 .* epsr_Ey)));
C_Ey_H = dt ./ ((1 + (sigma_Ey .* dt ./ (2 * eps0 .* epsr_Ey))) .* eps0 .* epsr_Ey .* ds);

C_Ez_E = (1 - (sigma_Ez .* dt ./ (2 * eps0 .* epsr_Ez))) ./ ...
         (1 + (sigma_Ez .* dt ./ (2 * eps0 .* epsr_Ez)));
C_Ez_H = dt ./ ((1 + (sigma_Ez .* dt ./ (2 * eps0 .* epsr_Ez))) .* eps0 .* epsr_Ez .* ds);

%% ABC initialization
Sc = c * dt / ds;              % for uniform cubic cells
abccoef = (Sc - 1) / (Sc + 1);

abcState = init_abc_3d(Nx, Ny, Nz);

%% Source setup
f0   = 0.25*f_max;
Nper = 10;
T0   = 1 / f0;
t0   = 0;

src_i = round(Nx/2);
src_j = 10;%round(Ny/2);
src_k = min(round(Nz/2), size(Ez,3));

%% Probe setup
probe_x = round(Nx/2);
probe_y = round(Ny/2);
probe_z = round(Nz/2);

srclog = zeros(1, num_steps+1);
EFieldProbe = zeros(3, num_steps+1);

%% Plot setup
figure;

%% Main time loop
if num_dims == 3
  for t = 0:num_steps
    % tic

    % =========================================================
    % Update H fields
    % =========================================================
    % Hx : (Nx, Ny-1, Nz-1)
    Hx(:, :, :) = ...
      C_Hx_H .* Hx + ...
      C_Hx_E .* ( ...
        (Ey(:, :, 2:Nz)   - Ey(:, :, 1:Nz-1)) - ...
        (Ez(:, 2:Ny, :)   - Ez(:, 1:Ny-1, :)) );

    % Hy : (Nx-1, Ny, Nz-1)
    Hy(:, :, :) = ...
      C_Hy_H .* Hy + ...
      C_Hy_E .* ( ...
        (Ez(2:Nx, :, :)   - Ez(1:Nx-1, :, :)) - ...
        (Ex(:, :, 2:Nz)   - Ex(:, :, 1:Nz-1)) );

    % Hz : (Nx-1, Ny-1, Nz)
    Hz(:, :, :) = ...
      C_Hz_H .* Hz + ...
      C_Hz_E .* ( ...
        (Ex(:, 2:Ny, :)   - Ex(:, 1:Ny-1, :)) - ...
        (Ey(2:Nx, :, :)   - Ey(1:Nx-1, :, :)) );

    % =========================================================
    % Update E fields
    % =========================================================
    % Ex interior : (Nx-1, Ny, Nz)
    Ex(:, 2:Ny-1, 2:Nz-1) = ...
      C_Ex_E(:, 2:Ny-1, 2:Nz-1) .* Ex(:, 2:Ny-1, 2:Nz-1) + ...
      C_Ex_H(:, 2:Ny-1, 2:Nz-1) .* ( ...
        (Hz(:, 2:Ny-1, 2:Nz-1) - Hz(:, 1:Ny-2, 2:Nz-1)) - ...
        (Hy(:, 2:Ny-1, 2:Nz-1) - Hy(:, 2:Ny-1, 1:Nz-2)) );

    % Ey interior : (Nx, Ny-1, Nz)
    Ey(2:Nx-1, :, 2:Nz-1) = ...
      C_Ey_E(2:Nx-1, :, 2:Nz-1) .* Ey(2:Nx-1, :, 2:Nz-1) + ...
      C_Ey_H(2:Nx-1, :, 2:Nz-1) .* ( ...
        (Hx(2:Nx-1, :, 2:Nz-1) - Hx(2:Nx-1, :, 1:Nz-2)) - ...
        (Hz(2:Nx-1, :, 2:Nz-1) - Hz(1:Nx-2, :, 2:Nz-1)) );

    % Ez interior : (Nx, Ny, Nz-1)
    Ez(2:Nx-1, 2:Ny-1, :) = ...
      C_Ez_E(2:Nx-1, 2:Ny-1, :) .* Ez(2:Nx-1, 2:Ny-1, :) + ...
      C_Ez_H(2:Nx-1, 2:Ny-1, :) .* ( ...
        (Hy(2:Nx-1, 2:Ny-1, :) - Hy(1:Nx-2, 2:Ny-1, :)) - ...
        (Hx(2:Nx-1, 2:Ny-1, :) - Hx(2:Nx-1, 1:Ny-2, :)) );

    % =========================================================
    % Source injection
    % =========================================================
    % tt = t * dt;
    % burst = (tt >= t0) && (tt < t0 + Nper*T0);
    % src = burst * sin(2*pi*f0*(tt - t0));
    % 
    % srclog(t+1) = src;
    % Ez(src_i, src_j, src_k) = Ez(src_i, src_j, src_k) + src;
    f0 = 0.2 * f_max;
    t0 = 4 / f0;
    
    tt = t * dt;
    arg = pi * f0 * (tt - t0);
    src = (1 - 2*arg^2) * exp(-arg^2);
    
    srclog(t+1) = src;
    Ez(src_i, src_j, src_k) = Ez(src_i, src_j, src_k) + src;

    % =========================================================
    % Apply first-order ABC to E field boundaries
    % =========================================================
    [Ex, Ey, Ez, abcState] = apply_abc_3d(Ex, Ey, Ez, abcState, abccoef);

    % =========================================================
    % Plot a mid-plane slice of Ex
    % =========================================================
    plot_k = min(round(Nz/2), size(Ez,3));
    imagesc(20*log10(abs(Ez(:,:,plot_k))));
    axis image;
    title(sprintf('t = %.6e s', tt));
    colorbar;
    clim([-150 0]);
    drawnow;

    % =========================================================
    % Probe logging
    % =========================================================
    ix_ex = min(probe_x, size(Ex,1));
    ix_ey = min(probe_x, size(Ey,1));
    ix_ez = min(probe_x, size(Ez,1));

    iy_ex = min(probe_y, size(Ex,2));
    iy_ey = min(probe_y, size(Ey,2));
    iy_ez = min(probe_y, size(Ez,2));

    iz_ex = min(probe_z, size(Ex,3));
    iz_ey = min(probe_z, size(Ey,3));
    iz_ez = min(probe_z, size(Ez,3));

    EFieldProbe(:, t+1) = [ ...
      Ex(ix_ex, iy_ex, iz_ex); ...
      Ey(ix_ey, iy_ey, iz_ey); ...
      Ez(ix_ez, iy_ez, iz_ez)];

    % toc
  end
end

%% Optional probe plot
time_axis = (0:num_steps) * dt;

figure;
plot(time_axis, EFieldProbe(1,:), 'DisplayName', 'Ex');
hold on;
plot(time_axis, EFieldProbe(2,:), 'DisplayName', 'Ey');
plot(time_axis, EFieldProbe(3,:), 'DisplayName', 'Ez');
xlabel('Time (s)');
ylabel('Field');
title('Electric Field Probe');
legend;
grid on;

figure;
plot(time_axis, srclog, 'DisplayName','Source Current');
hold on;
%diff does a forward differencing
srclog_pad = [0, srclog, 0]; %zero padding enables central difference
% if the source is nonzero at the edges this will be bad
diff_srclog = (srclog_pad(3:end) - srclog_pad(1:end-2));
plot(time_axis,diff_srclog./max(abs(diff_srclog)), 'DisplayName','dI/dt');
xlabel('Time (s)');
ylabel('Source amplitude');
title('Injected source');
legend;
grid on;

%% FFT Analysis of Probe vs Source

figure;
fs = 1/dt;

xdft_1 = fft(diff_srclog);
xdft_1 = xdft_1(1:(num_steps+1)/2+1);
psdx_1 = (1/(fs*(num_steps+1))) * abs(xdft_1).^2;
psdx_1(2:end-1) = 2*psdx_1(2:end-1);

freq = 0:fs/length(diff_srclog):fs/2;

xdft_2 = fft(EFieldProbe(3,:));
xdft_2 = xdft_2(1:(num_steps+1)/2+1);
psdx_2 = (1/(fs*(num_steps+1))) * abs(xdft_2).^2;
psdx_2(2:end-1) = 2*psdx_2(2:end-1);


plot(freq/1E6,pow2db(psdx_1)-max(pow2db(psdx_1)),'DisplayName','Injected Signal (dI/dt)')
hold on
plot(freq/1E6,pow2db(psdx_2)-max(pow2db(psdx_2)),'DisplayName','Ez Probe')
grid on
title("Normalized PSD: Ricker Wavelet Source")
xlabel("Frequency (MHz)")
ylabel("Power/Frequency (dB/Hz)")

%% Speed of Light in Sim Estimation

figure;
[r,lags] = xcorr(EFieldProbe(3,:),diff_srclog);
stem(lags,r/(max(abs(r))));
probe_delay = lags(find(max(abs(r))==r))*dt;
probe_dist = sqrt( (src_i-probe_x)^2 + (src_j-probe_y)^2 + (src_k-probe_z)^2) * dx;
sim_c = probe_dist/probe_delay;
title('Normalized Cross Correlation Between Source (dI/dt) and Probe Signal')
text(0.95, 0.95, 'Approximate Simulation Speed vs. Light Speed = ' + string(sim_c/c), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')


%% =======================================================================
%% Local functions
%% =======================================================================

function abcState = init_abc_3d(Nx, Ny, Nz)
% Initialize stored boundary-history terms for first-order ABC.

  abcState.Eyx0 = zeros(1,    Ny-1, Nz);
  abcState.Ezx0 = zeros(1,    Ny,   Nz-1);

  abcState.Eyx1 = zeros(1,    Ny-1, Nz);
  abcState.Ezx1 = zeros(1,    Ny,   Nz-1);

  abcState.Exy0 = zeros(Nx-1, 1,    Nz);
  abcState.Ezy0 = zeros(Nx,   1,    Nz-1);

  abcState.Exy1 = zeros(Nx-1, 1,    Nz);
  abcState.Ezy1 = zeros(Nx,   1,    Nz-1);

  abcState.Exz0 = zeros(Nx-1, Ny,   1);
  abcState.Eyz0 = zeros(Nx,   Ny-1, 1);

  abcState.Exz1 = zeros(Nx-1, Ny,   1);
  abcState.Eyz1 = zeros(Nx,   Ny-1, 1);
end

function [Ex, Ey, Ez, abcState] = apply_abc_3d(Ex, Ey, Ez, abcState, abccoef)
% First-order absorbing boundary condition for 3D Yee-grid E-fields.

  if abccoef == 0
    error('apply_abc_3d: abccoef is zero. Initialize abccoef before calling.');
  end

  [Nx_m1, Ny,    Nz   ] = size(Ex); %#ok<ASGLU>
  [Nx,    Ny_m1, Nz2  ] = size(Ey); %#ok<ASGLU>
  [Nx2,   Ny2,   Nz_m1] = size(Ez); %#ok<ASGLU>

  if Nx_m1 ~= Nx-1 || Nz2 ~= Nz || Nx2 ~= Nx || Ny2 ~= Ny || Nz_m1 ~= Nz-1 || Ny_m1 ~= Ny-1
    error('Field array sizes are inconsistent with a standard Yee grid.');
  end

  % ---------------------------------------------------------
  % ABC at x0
  % ---------------------------------------------------------
  mm = 1;

  Ey(mm, 1:Ny-1, 1:Nz) = abcState.Eyx0 + ...
    abccoef .* (Ey(mm+1, 1:Ny-1, 1:Nz) - Ey(mm, 1:Ny-1, 1:Nz));
  abcState.Eyx0 = Ey(mm+1, 1:Ny-1, 1:Nz);

  Ez(mm, 1:Ny, 1:Nz-1) = abcState.Ezx0 + ...
    abccoef .* (Ez(mm+1, 1:Ny, 1:Nz-1) - Ez(mm, 1:Ny, 1:Nz-1));
  abcState.Ezx0 = Ez(mm+1, 1:Ny, 1:Nz-1);

  % ---------------------------------------------------------
  % ABC at x1
  % ---------------------------------------------------------
  mm = Nx;

  Ey(mm, 1:Ny-1, 1:Nz) = abcState.Eyx1 + ...
    abccoef .* (Ey(mm-1, 1:Ny-1, 1:Nz) - Ey(mm, 1:Ny-1, 1:Nz));
  abcState.Eyx1 = Ey(mm-1, 1:Ny-1, 1:Nz);

  Ez(mm, 1:Ny, 1:Nz-1) = abcState.Ezx1 + ...
    abccoef .* (Ez(mm-1, 1:Ny, 1:Nz-1) - Ez(mm, 1:Ny, 1:Nz-1));
  abcState.Ezx1 = Ez(mm-1, 1:Ny, 1:Nz-1);

  % ---------------------------------------------------------
  % ABC at y0
  % ---------------------------------------------------------
  nn = 1;

  Ex(1:Nx-1, nn, 1:Nz) = abcState.Exy0 + ...
    abccoef .* (Ex(1:Nx-1, nn+1, 1:Nz) - Ex(1:Nx-1, nn, 1:Nz));
  abcState.Exy0 = Ex(1:Nx-1, nn+1, 1:Nz);

  Ez(1:Nx, nn, 1:Nz-1) = abcState.Ezy0 + ...
    abccoef .* (Ez(1:Nx, nn+1, 1:Nz-1) - Ez(1:Nx, nn, 1:Nz-1));
  abcState.Ezy0 = Ez(1:Nx, nn+1, 1:Nz-1);

  % ---------------------------------------------------------
  % ABC at y1
  % ---------------------------------------------------------
  nn = Ny;

  Ex(1:Nx-1, nn, 1:Nz) = abcState.Exy1 + ...
    abccoef .* (Ex(1:Nx-1, nn-1, 1:Nz) - Ex(1:Nx-1, nn, 1:Nz));
  abcState.Exy1 = Ex(1:Nx-1, nn-1, 1:Nz);

  Ez(1:Nx, nn, 1:Nz-1) = abcState.Ezy1 + ...
    abccoef .* (Ez(1:Nx, nn-1, 1:Nz-1) - Ez(1:Nx, nn, 1:Nz-1));
  abcState.Ezy1 = Ez(1:Nx, nn-1, 1:Nz-1);

  % ---------------------------------------------------------
  % ABC at z0
  % ---------------------------------------------------------
  pp = 1;

  Ex(1:Nx-1, 1:Ny, pp) = abcState.Exz0 + ...
    abccoef .* (Ex(1:Nx-1, 1:Ny, pp+1) - Ex(1:Nx-1, 1:Ny, pp));
  abcState.Exz0 = Ex(1:Nx-1, 1:Ny, pp+1);

  Ey(1:Nx, 1:Ny-1, pp) = abcState.Eyz0 + ...
    abccoef .* (Ey(1:Nx, 1:Ny-1, pp+1) - Ey(1:Nx, 1:Ny-1, pp));
  abcState.Eyz0 = Ey(1:Nx, 1:Ny-1, pp+1);

  % ---------------------------------------------------------
  % ABC at z1
  % ---------------------------------------------------------
  pp = Nz;

  Ex(1:Nx-1, 1:Ny, pp) = abcState.Exz1 + ...
    abccoef .* (Ex(1:Nx-1, 1:Ny, pp-1) - Ex(1:Nx-1, 1:Ny, pp));
  abcState.Exz1 = Ex(1:Nx-1, 1:Ny, pp-1);

  Ey(1:Nx, 1:Ny-1, pp) = abcState.Eyz1 + ...
    abccoef .* (Ey(1:Nx, 1:Ny-1, pp-1) - Ey(1:Nx, 1:Ny-1, pp));
  abcState.Eyz1 = Ey(1:Nx, 1:Ny-1, pp-1);
end