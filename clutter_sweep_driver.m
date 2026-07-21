%% run_sweep_driver.m
% Run the FDTD suite for the target-in-clutter theta sweep.
%
% Loads the FIXED clutter + sweep from clutter_sweep.mat (built by
% setup_clutter_and_sweep.m), then for each of the 11 theta positions appends
% the target box to the SAME clutter and runs one full FDTD. Also runs a
% background (empty) and a clutter-only (no target) case for the integration
% post-processor.
%
% Only the target moves between runs; the clutter is identical every run --
% which is what lets the later coherent integration separate target from clutter.
%
% Requires: spherical_fdtd_cpml_allwalls_SIBC_v2 on the path, and clutter_sweep.mat (run setup_clutter_and_sweep.m).
%
% NOTE: each position is a full FDTD run. Set run_positions to a subset
% (e.g. 6) to validate the pipeline before the full sweep.

clear; clc;

%% -------------------- load the fixed clutter + sweep --------------------
S = load('clutter_sweep.mat');
clutter       = S.clutter;            % Nships x 6  (FIXED across all runs)
theta_centers = S.theta_centers;      % 11 target theta-centers
Ltgt_the      = S.Ltgt_the;
tgt_k1        = S.tgt_k1;
tgt_k2        = S.tgt_k2;
i1 = S.i1; i2 = S.i2;
nPos = numel(theta_centers);

fprintf('Loaded clutter_sweep.mat: %d clutter ships, %d target positions.\n', ...
        size(clutter,1), nPos);

%% -------------------- which runs to do --------------------
run_positions  = 1:nPos;     % <-- set to e.g. 6 to test just the center first
do_background  = true;       % empty-domain run (no clutter, no target)
do_clutteronly = true;       % fixed clutter, NO target (background for subtraction)

out_dir = 'clutter_runs';
if ~exist(out_dir,'dir'); mkdir(out_dir); end

%% -------------------- background run (empty) --------------------
if do_background
    cfg = struct();
    cfg.box_list   = zeros(0,6);
    cfg.viz_enable = false;
    cfg.save_path  = fullfile(out_dir, 'run_background.mat');
    fprintf('\n=== background run (empty domain) ===\n');
    visualize_boxes(cfg.box_list)
    spherical_fdtd_cpml_allwalls_SIBC(cfg);
end
close all
%% -------------------- clutter-only run (fixed clutter, no target) --------------------
if do_clutteronly
    cfg = struct();
    cfg.box_list   = clutter;                      % clutter only
    cfg.viz_enable = false;
    cfg.save_path  = fullfile(out_dir, 'run_clutteronly.mat');
    fprintf('\n=== clutter-only run (fixed clutter, no target) ===\n');
    visualize_boxes(cfg.box_list)
    spherical_fdtd_cpml_allwalls_SIBC(cfg);
end
close all
%% -------------------- target-in-clutter, one run per theta position --------------------
for p = run_positions
    jc  = theta_centers(p);
    tj1 = jc - floor(Ltgt_the/2);
    tj2 = tj1 + Ltgt_the;
    target_box = [i1 i2 tj1 tj2 tgt_k1 tgt_k2];

    cfg = struct();
    cfg.box_list   = [clutter; target_box];        % SAME clutter + moved target
    cfg.viz_enable = false;
    cfg.save_path  = fullfile(out_dir, sprintf('run_pos%02d_j%04d.mat', p, jc));
    visualize_boxes(cfg.box_list)
    fprintf('\n=== position %d/%d: target theta-center j=%d ===\n', p, nPos, jc);
    spherical_fdtd_cpml_allwalls_SIBC(cfg);
end
close all
fprintf('\nSweep complete. Runs saved to %s/.\n', out_dir);
fprintf('Post-process with integrate_clutter_suite.m\n');
