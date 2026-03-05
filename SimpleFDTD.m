%Constants
mu0 = (4*pi)*10^-7;
ep0 = 8.854E-12;
c = 1/sqrt(ep0*mu0);
f_max = 500E3;
lam_max = c/f_max; % grid spacing must be an order of magnitude smaller

h_ionosphere = 150E3;
max_dist = 2000E3;
max_width = 100E3;

%Grid Dispersion Error req ds << lambda
% Max freq of interest is 500 kHz so 60 meter max grid spacing
%IN THE FUTURE WE MAY CREATE A CUSTOM GRID THAT WE CHECK AGAINST THIS
%REQUIREMENT
ds = lam_max/10; %meters

dx = ds; %meters
dy = ds; %meters
dz = ds; %meters

%Grid Size Parameters: For testing this is medium size (This will need to
%go to the full size of the system eventually)

Nx = round(max_dist/dx);
Ny = ceil(max_width/dy);
Nz = ceil(h_ionosphere/dz);

%CFL Stability Criterion (dt < 1/(ds/(c*sqrt(n))) for cubic grid

dt = 0.9 * dx/(c*sqrt(2));

dt_steps = 30000; % 30k steps at dt of 6.4E-08 is about 600km net propagation



%Material properties - Currently set to free space

eps_r = ones(Nx,Ny);
mu_r = ones(Nx,Ny);

%%%% TEXTBOOK USES Face center E field and the notes i had use face center
%%%% H field
%Initialize fields in 2D

Ez = zeros(Nx,Ny); %go to the edges of the box in yee lattice

% H are face centered so they will be evaluated at half steps and thus cannot exceed boundary
Hx = zeros(Nx,Ny-1); % evaluated about the offset y dim so they are missing a cell in that dim
Hy = zeros(Nx-1,Ny); % evaluated about the offset x dim 

% Initialize Source Loc (center of environment)

is = round(Nx/2);
js = round(Ny/2);

% --- Simple absorbing sponge layer (loss term on Ez) ---
nSponge = 20;               % sponge thickness [cells]
sigma = zeros(Nx,Ny);       % conductivity-like damping term


%Start Loop

for t = 1:dt_steps
    
    %calc H fields a half step ahead using time t E field
    Hx = 
    Hy = 

    %calc E fields a step ahead using H field a half step ahead
    Ez_Now() = Ez_Prev(previous) + (dt./epsilon())