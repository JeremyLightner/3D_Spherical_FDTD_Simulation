%% setup_clutter_and_sweep.m
% Feature 1: build 11 theta positions for the target sweep at fixed phi=45 km,
%            6th point centered, spanning just past the low-theta keep-out to
%            just before the high-theta keep-out.
% Feature 2: generate random clutter that AVOIDS the target's traversal corridor
%            (re-placing any ship landing within 1 cell of it).
% Feature 3: save the clutter (and sweep) so all 11 runs reuse the SAME clutter.
%
% Assumes these are already defined (as in example_ext_cfg.m):
%   ds, R_Earth, del_the, del_phi, lat0_rad, phi_start, Nr, Nthe, Nphi,
%   cpml_cells, fraunhofer_cells, margin_cells, phi_min, phys2box

rng(98457);                      % FIXED seed -> reproducible clutter (Feature 3)

ds        = (299792458/500e3)/20;      % lambda/20 at f_max=500 kHz
R_Earth   = 6371e3;
Nr=168; Nthe=235; Nphi=1671;
% Nr = 92; Nthe = 112; Nphi = 1336;      % <-- SET to your grid
del_the   = ds/(R_Earth+5e3);        % <-- SET to sim's del_the
del_phi   = ds/(R_Earth+5e3);                % <-- SET to sim's del_phi
lat0_rad  = 0; phi_start = 0;          % <-- SET to sim's values

% ---- keep-out margin from each CPML wall ----
% Boxes must sit at least (CPML thickness + Fraunhofer distance) inside every
% lateral wall, so the scatterer's reactive near field has decayed to the far
% field before reaching the PML (otherwise the PML absorbs a near field it was
% not designed for -> spurious reflection).
cpml_cells      = 10;                  % <-- SET to sim's CPML thickness
fraunhofer_m    = 147;                 % far-field distance of a ship this size [m]
fraunhofer_cells = ceil(fraunhofer_m / ds);           % ~5 cells
margin_cells    = cpml_cells + fraunhofer_cells;      % total keep-out per wall
fprintf('Keep-out per wall: %d CPML + %d Fraunhofer = %d cells (%.0f m)\n', ...
        cpml_cells, fraunhofer_cells, margin_cells, margin_cells*ds);

phys2box = @(h_m, Lthe_m, Lphi_m, range_km, jc) deal_box(h_m, Lthe_m, Lphi_m, ...
             range_km, jc, ds, R_Earth, del_the, del_phi, lat0_rad, phi_start, ...
             Nr, Nthe, Nphi, margin_cells);
phi_min = 235;

opts.cpml_cells       = cpml_cells;      % red CPML band
opts.fraunhofer_cells = fraunhofer_cells; % orange ship-margin band (or inferred: margin - cpml)
opts.margin_cells     = margin_cells;    % total wall margin
opts.phi_min          = phi_min; % MTLE cutoff: blue zone + labeled line
opts.Nthe             = Nthe;        
opts.Nphi             = Nphi;          
opts.ds               = ds;           
% opts.src_j = round(Nthe/2); opts.src_k = 117;   % <-- optional: draw the source star

%% ---- target geometry (must match the target you sweep) ----
tgt_h_m    = 30;   tgt_the_m = 200;  tgt_phi_m = 30;
tgt_range_km = 45;                                  % fixed range of the sweep
Ltgt_the = max(1, round(tgt_the_m / (R_Earth*del_the)));   % target theta cells
Ltgt_phi = max(1, round(tgt_phi_m / (R_Earth*del_phi)));   % target phi cells

%% ============================================================
%% Feature 1: 11 theta positions, 6th centered, spanning keep-outs
%% ============================================================
% Allowed theta-CENTER range so the whole target footprint clears the margins:
jc_lo = margin_cells + 1 + ceil(Ltgt_the/2);        % just past low-theta keep-out
jc_hi = Nthe - margin_cells - ceil(Ltgt_the/2);     % just before high-theta keep-out

nPos = 11;
theta_centers = round(linspace(jc_lo, jc_hi, nPos));

% Force the 6th point to be the exact center (linspace centers it for odd nPos,
% but enforce explicitly so it is exact even after rounding).
center_jc = round((jc_lo + jc_hi)/2);
theta_centers(6) = center_jc;

fprintf('Target theta sweep (phi=%g km, %d positions):\n', tgt_range_km, nPos);
fprintf('  jc range [%d, %d], 6th (center) = %d\n', jc_lo, jc_hi, center_jc);
disp(theta_centers);

%% ---- the target's phi extent at 45 km (fixed for all positions) ----
dphi = tgt_range_km*1e3 / (R_Earth*cos(lat0_rad));
phi_node = phi_start + (0:Nphi-1)*del_phi;
[~, tgt_kc] = min(abs(phi_node - (phi_start + dphi)));
tgt_k1 = tgt_kc - floor(Ltgt_phi/2);
tgt_k2 = tgt_k1 + Ltgt_phi;

%% ============================================================
%% Feature 2: target-traversal corridor (keep-out for clutter)
%% ============================================================
% As the target sweeps theta at fixed phi, it occupies the phi band
% [tgt_k1, tgt_k2] across the FULL theta span of the sweep. Clutter must avoid
% this corridor (+1 cell on every side).
corridor_j1 = min(theta_centers) - floor(Ltgt_the/2) - 1;   % lowest target edge -1
corridor_j2 = max(theta_centers) + ceil (Ltgt_the/2) + 1;   % highest target edge +1
corridor_k1 = tgt_k1 - 1;
corridor_k2 = tgt_k2 + 1;
fprintf('Target corridor (clutter keep-out): theta [%d,%d], phi [%d,%d]\n', ...
        corridor_j1, corridor_j2, corridor_k1, corridor_k2);

%% ============================================================
%% Random clutter avoiding the corridor (re-place on overlap)
%% ============================================================
numClutter = 20;
ship_h_m   = 30;  ship_the_m = 200;  ship_phi_m = 30;
Lthe = max(1, round(ship_the_m / (R_Earth*del_the)));
Lphi = max(1, round(ship_phi_m / (R_Earth*del_phi)));
i1 = 1; i2 = min(Nr, 1 + max(1,round(ship_h_m/ds)));
gap = 1;                                            % ship-to-ship / corridor gap
max_tries = 5000;

clutter = zeros(0,6);
for s = 1:numClutter
    placed = false;
    for tries = 1:max_tries
        if randi([0 1]); Lj=Lthe; Lk=Lphi; else; Lj=Lphi; Lk=Lthe; end
        jc = randi([margin_cells+1+ceil(Lj/2), Nthe-margin_cells-ceil(Lj/2)]);
        % phi lower bound respects BOTH the wall margin and the MTLE cutoff
        kc_lo = max(margin_cells+1+ceil(Lk/2), phi_min+ceil(Lk/2));
        kc    = randi([kc_lo, Nphi-margin_cells-ceil(Lk/2)]);
        j1=jc-floor(Lj/2); j2=j1+Lj; k1=kc-floor(Lk/2); k2=k1+Lk;

        % --- reject if within 1 cell of the target corridor ---
        if ~(j2 < corridor_j1-gap || j1 > corridor_j2+gap || ...
             k2 < corridor_k1-gap || k1 > corridor_k2+gap)
            continue;                               % overlaps corridor -> re-place
        end
        % --- reject if overlapping the source column (if defined) ---
        if exist('src_j','var') && exist('src_k','var')
            if k1<=src_k && src_k<=k2 && j1<=src_j && src_j<=j2, continue; end
        end
        % --- reject if overlapping another ship (with gap) ---
        ok = true;
        for q = 1:size(clutter,1)
            b = clutter(q,:);
            if ~(j2 < b(3)-gap || j1 > b(4)+gap || k2 < b(5)-gap || k1 > b(6)+gap)
                ok = false; break;
            end
        end
        if ~ok, continue; end

        clutter(end+1,:) = [i1 i2 j1 j2 k1 k2]; %#ok<SAGROW>
        placed = true; break;
    end
    if ~placed
        warning('Placed only %d of %d ships (corridor + packing too tight).', ...
                size(clutter,1), numClutter);
        break;
    end
end
numClutter = size(clutter,1);
fprintf('Placed %d clutter ships (all clear of the target corridor).\n', numClutter);

%% ============================================================
%% Feature 3: save the fixed clutter + sweep for all 11 runs
%% ============================================================
save('clutter_sweep.mat', 'clutter', 'theta_centers', 'tgt_k1', 'tgt_k2', ...
     'Ltgt_the', 'Ltgt_phi', 'i1', 'i2', 'tgt_range_km', ...
     'corridor_j1','corridor_j2','corridor_k1','corridor_k2', ...
     'Nthe','Nphi','margin_cells','phi_min','cpml_cells','fraunhofer_cells','ds');
fprintf('Saved clutter + sweep to clutter_sweep.mat (reused by all runs).\n');

%% ---- visualize the layout (clutter + corridor + center target position) ----
jc = theta_centers(6);
tj1 = jc-floor(Ltgt_the/2); tj2 = tj1+Ltgt_the;
target_center = [i1 i2 tj1 tj2 tgt_k1 tgt_k2];

vopts.Nthe = Nthe; vopts.Nphi = Nphi;
vopts.cpml_cells = cpml_cells; vopts.fraunhofer_cells = fraunhofer_cells;
vopts.margin_cells = margin_cells; vopts.phi_min = phi_min; vopts.ds = ds;
vopts.theta_centers = theta_centers;      % draws the target travel path
vopts.target_k1 = tgt_k1; vopts.target_k2 = tgt_k2;
vopts.target_rows = size(clutter,1)+1;    % highlight the target (last box)

visualize_boxes([clutter; target_center], vopts);
title(sprintf('Fixed clutter (%d ships) + target sweep, %d theta positions', ...
      numClutter, nPos));
