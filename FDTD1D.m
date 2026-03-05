%% 1D FDTD - Jeremiah Lightner 2026
clear all
clc

eps0 = 8.854E-12;
mu0 = 4*pi*10^-7;
c = 1/sqrt(mu0*eps0);

num_dims = 1;
f_max = 1E6;

lam_min = c/f_max;

del_x = lam_min/10; %Grid Dispersion Error req ds << lambda
del_t = 0.99 * del_x/(c*sqrt(num_dims)); % CFL Stability Criterion (dt < 1/(ds/(c*sqrt(n))) for cubic grid

num_t = 10000;

Nx = 201;
Ny = 201;
Nz = 201;

eps = eps0; %this will need to become a 2d matrix to capture information
mu = mu0;

Ey = zeros(1,Ny);
Hz = zeros(1,Nz);

%Precompute the update coefficients
mEy = del_t/(eps*del_x); % multiplying c normalizes the fields (should reduce floating point error)
mHz = del_t/(mu*del_x);
% mEy = c*del_t/eps;
% mHz = c*del_t/mu;

%Precompute injected source then step through it
source = gauspuls(0:del_t:(num_t*del_t)-100*del_t,3E6,0.5);
source = [zeros(1,100),source];
% t_range = 0:20;
% source =  sin(2*pi*f_max*del_t*t_range);

%prog loop
for t = 1:num_t

    % Hz(1:Nz-1) = Hz(1:Nz-1) + mHz*(Ey(2:Ny)-Ey(1:Ny-1));
    % Ey(2:Ny-1) = Ey(2:Ny-1) + mEy*(Hz(1:Nz-2)-Hz(2:Nz-1))  ;

    %% H field loop
    Hz(1:Nz-1) = Hz(1:Nz-1)+(del_t/(eps0*del_x))*diff(Ey);
    Hz(Nz) = Hz(Nz) + (del_t/(eps0*del_x))*(-Ey(Nz));
    % for ii = 1:Nz-1
    %     Hz(ii) = Hz(ii)+mHz*(Ey(ii+1)-Ey(ii));
    % end   
    % Hz(Nz) = Hz(Nz) + mHz*(0-Ey(Nz));

    %% E field loop
    Ey(1) = Ey(1) + (del_t/(mu0*del_x))*Hz(1);
    Ey(2:Ny) = Ey(2:Ny)+(del_t/(mu0*del_x))*diff(Hz(1:Ny));
    % Ey(1) = Ey(1) + mHz*(Hz(1));
    % for ii = 2:Ny
    %     Ey(ii) = Ey(ii)+mHz*(Hz(ii)-Hz(ii-1));
    % end


    % if t<=length(source)
    Ey(round(Ny/2)) = Ey(round(Ny/2)) + exp(-.5*((t-40)/8)^2);
    % end

    plot(1:Ny,Ey); title(t*del_t);drawnow;
    % pause(0.2)
end

