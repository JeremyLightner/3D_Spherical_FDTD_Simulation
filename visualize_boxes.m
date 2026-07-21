%% visualize_boxes.m
% Top-down (theta-phi plan view) of a PEC box_list, showing footprints, the CPML
% walls, the 150 m guard margin, the source column, and (optionally) the target
% travel path. Use this to eyeball clutter/target placement before running.
%
% USAGE:
%   visualize_boxes(box_list)                       % boxes only
%   visualize_boxes(box_list, opts)                 % with grid/source/target
%
% box_list : Nbox x 6, rows [i1 i2 j1 j2 k1 k2] (i=radial, j=theta, k=phi).
% opts (struct, all optional):
%   .Nthe,.Nphi          grid extent (draws the domain box)
%   .cpml_cells          CPML thickness (draws PML walls)
%   .cpml_cells          CPML thickness (drawn as the red wall band)
%   .fraunhofer_cells    ship-Fraunhofer margin beyond CPML (light-orange band);
%                        if omitted, inferred as margin_cells - cpml_cells
%   .margin_cells        total wall margin = CPML + ship-Fraunhofer
%   .phi_min             MTLE source-Fraunhofer cutoff: shaded blue zone + line
%                        at phi=phi_min (no clutter east of it)
%   .guard_cells         (legacy alias; use margin_cells instead)
%   .src_j,.src_k        source location (drawn as a star)
%   .target_rows         indices into box_list that are the TARGET (highlighted)
%   .theta_centers       target sweep theta-centers (draws the travel path)
%   .target_k1,.target_k2  target phi extent (for the path band)
%   .ds                  cell size [m] for a secondary axis in meters
%
% Axes: horizontal = phi index (k, range/propagation), vertical = theta index
% (j, cross-range). This matches a plan view looking down on the ocean surface.

function visualize_boxes(box_list, opts)
if nargin < 2, opts = struct(); end

figure('Color','w'); hold on;

%% ---- domain + keep-out zones ----
if isfield(opts,'Nthe') && isfield(opts,'Nphi')
    Nthe = opts.Nthe; Nphi = opts.Nphi;

    % Force the axes to show the WHOLE simulation space (1..Nphi, 1..Nthe).
    xlim([1 Nphi]); ylim([1 Nthe]);

    % full domain outline
    rectangle('Position',[1 1 Nphi-1 Nthe-1], 'EdgeColor',[0 0 0], 'LineWidth',1.2);

    c = 0; frn = 0;
    if isfield(opts,'cpml_cells') && ~isempty(opts.cpml_cells), c = opts.cpml_cells; end
    if isfield(opts,'fraunhofer_cells') && ~isempty(opts.fraunhofer_cells)
        frn = opts.fraunhofer_cells;
    elseif isfield(opts,'margin_cells') && ~isempty(opts.margin_cells) && c>0
        frn = opts.margin_cells - c;                 % ship-Fraunhofer = margin - CPML
    end
    m = c + frn;                                     % total wall margin

    % --- CPML walls (innermost-red band, all four sides) ---
    if c > 0
        pml_col = [1 0.80 0.80];
        rectangle('Position',[1 1 c Nthe-1],       'FaceColor',pml_col,'EdgeColor','none'); % east
        rectangle('Position',[Nphi-c 1 c Nthe-1],  'FaceColor',pml_col,'EdgeColor','none'); % west
        rectangle('Position',[1 1 Nphi-1 c],       'FaceColor',pml_col,'EdgeColor','none'); % south
        rectangle('Position',[1 Nthe-c Nphi-1 c],  'FaceColor',pml_col,'EdgeColor','none'); % north
    end

    % --- ship-Fraunhofer band (between CPML edge and full margin) ---
    if frn > 0
        frn_col = [1 0.92 0.80];                     % lighter orange band
        % east / west strips
        rectangle('Position',[c 1 frn Nthe-1],        'FaceColor',frn_col,'EdgeColor','none');
        rectangle('Position',[Nphi-m 1 frn Nthe-1],   'FaceColor',frn_col,'EdgeColor','none');
        % south / north strips
        rectangle('Position',[1 c Nphi-1 frn],        'FaceColor',frn_col,'EdgeColor','none');
        rectangle('Position',[1 Nthe-m Nphi-1 frn],   'FaceColor',frn_col,'EdgeColor','none');
    end

    % --- allowed-region outline (inside all wall margins) ---
    if m > 0
        jlo=m+1; jhi=Nthe-m; klo=m+1; khi=Nphi-m;
        rectangle('Position',[klo jlo khi-klo jhi-jlo], ...
                  'EdgeColor',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
    end

    % --- MTLE source-Fraunhofer cutoff (phi_min): no clutter east of this ---
    if isfield(opts,'phi_min') && ~isempty(opts.phi_min)
        pm = opts.phi_min;
        mtle_col = [0.82 0.87 1.0];                  % blue-ish source zone
        % shade the whole strip from the domain start up to phi_min (over full theta)
        rectangle('Position',[1 1 max(pm-1,1) Nthe-1], ...
                  'FaceColor',mtle_col,'EdgeColor','none');
        plot([pm pm],[1 Nthe],'-','Color',[0.1 0.2 0.7],'LineWidth',1.5);
        text(pm, Nthe*0.5, sprintf(' MTLE source-Fraunhofer cutoff (\\phi=%d) ', pm), ...
             'Color',[0.1 0.2 0.7], 'Rotation',90, 'FontSize',8, ...
             'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end

%% ---- target travel path (theta sweep at fixed phi band) ----
if isfield(opts,'theta_centers') && isfield(opts,'target_k1') && isfield(opts,'target_k2')
    tc = opts.theta_centers;
    kmid = 0.5*(opts.target_k1 + opts.target_k2);
    % vertical dotted line marking the sweep, with tick marks at each position
    plot([kmid kmid], [min(tc) max(tc)], ':', 'Color',[0.2 0.2 0.8], 'LineWidth',1.2);
    plot(kmid*ones(size(tc)), tc, '.', 'Color',[0.2 0.2 0.8], 'MarkerSize',8);
    text(kmid, max(tc), '  target path', 'Color',[0.2 0.2 0.8], ...
         'Rotation',90, 'VerticalAlignment','bottom', 'FontSize',8);
end

%% ---- boxes ----
target_rows = [];
if isfield(opts,'target_rows'), target_rows = opts.target_rows; end
for b = 1:size(box_list,1)
    r = box_list(b,:);
    j1=r(3); j2=r(4); k1=r(5); k2=r(6);
    w = k2-k1; h = j2-j1;
    is_target = ismember(b, target_rows);
    if is_target
        rectangle('Position',[k1 j1 w h], 'FaceColor',[1 0.6 0.1], ...
                  'EdgeColor',[0.6 0.3 0],'LineWidth',1.5);
    else
        rectangle('Position',[k1 j1 w h], 'FaceColor',[0.6 0.75 0.95], ...
                  'EdgeColor',[0.2 0.35 0.6],'LineWidth',0.8);
    end
end

%% ---- source ----
if isfield(opts,'src_j') && isfield(opts,'src_k')
    plot(opts.src_k, opts.src_j, 'p', 'MarkerSize',16, ...
         'MarkerFaceColor',[1 0.9 0], 'MarkerEdgeColor','k');
    text(opts.src_k, opts.src_j, '  source', 'FontSize',8, 'VerticalAlignment','middle');
end

%% ---- labels / cosmetics ----
% Show the whole simulation space. Do NOT use 'axis equal' -- the domain is very
% elongated (Nphi >> Nthe), so equal aspect would make it an unreadable sliver.
if isfield(opts,'Nthe') && isfield(opts,'Nphi')
    xlim([1 opts.Nphi]); ylim([1 opts.Nthe]);
else
    axis tight;
end
xlabel('\phi index (k)  -- range / propagation');
ylabel('\theta index (j)  -- cross-range');
ttl = sprintf('PEC box layout: %d boxes', size(box_list,1));
if ~isempty(target_rows), ttl = sprintf('%s (%d target)', ttl, numel(target_rows)); end
title(ttl);
grid on; box on;

% optional secondary axes in meters
if isfield(opts,'ds') && ~isempty(opts.ds)
    ax = gca;
    xl = xlim; yl = ylim;
    fprintf('Plan view spans phi %.0f m x theta %.0f m (ds=%.2f m)\n', ...
            (xl(2)-xl(1))*opts.ds, (yl(2)-yl(1))*opts.ds, opts.ds);
end

hold off;
end
