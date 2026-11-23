function animate_trajectory(fh, params, opts)
% Animation of a trajectory segment (VTOL take-off, time-dependent cruise,
% or a multiphase history) from its flight history (fh from flight_history).
% Left: the aircraft, pitched by theta with the nacelles tilted by
% epsilon1/2 (CAD side views if fuselage_sideview.png and
% nacelle_sideview.png are on the path, grey boxes otherwise), with every
% force drawn at its point of application and the pitching moments as arcs
% (aero moment of each surface at its station, net moment at the CG).
% Right: the trajectory with the current position. Bottom: force and
% moment time histories with a time cursor.
% Arrows share one scale over the whole trajectory, and so do arcs, so their
% evolution is comparable between frames. Saved as Images/<name>.mp4.
%   opts: fps (25, or 10 for segments longer than 60 s),
%         speed (1 = real time, 0.5 = half speed; default 1),
%         duration [s] of the video for segments longer than 60 s (e.g. a
%         cruise), which sets the speed (default 12 s),
%         save (true), play (true: real-time pauses when the figure is
%         visible), name (default [fh.segment '_animation'])
% Every frame is rendered and captured (~0.7 s each in batch), so long
% segments default to 12 s at 10 fps: 120 frames, ~1.5 min.

if nargin < 3, opts = struct(); end
long = fh.tf > 60;
if ~isfield(opts, 'fps'),      opts.fps      = 25 - 15*long; end
if ~isfield(opts, 'duration'), opts.duration = 12;           end
if ~isfield(opts, 'speed')
    opts.speed = 1;
    if long, opts.speed = fh.tf/opts.duration; end
end
if ~isfield(opts, 'save'),  opts.save  = true; end
if ~isfield(opts, 'play'),  opts.play  = true; end
if ~isfield(opts, 'name'),  opts.name  = [fh.segment '_animation']; end

%% Data at the frame times

tq = fh.t(1):opts.speed/opts.fps:fh.t(end);
if tq(end) < fh.t(end), tq(end+1) = fh.t(end); end
n_frames = numel(tq);

S = struct();
for f = {'x', 'h', 'theta', 'epsilon1', 'epsilon2', 'V', 'gamma', 'T1', 'T2', 'H1', 'H2', ...
         'D_nac1', 'D_nac2', 'alpha_nac1', 'alpha_nac2', 'L_wing', 'D_wing', 'L_fus', 'D_fus', ...
         'L_tail', 'D_tail', 'MA_wing', 'MA_fus', 'MA_tail', 'M'}
    S.(f{1}) = interp1(fh.t, fh.(f{1}), tq);
end
W = params.mass*params.g;

% Time in seconds, or minutes for long segments (cruise)
t_div = 1; t_unit = 's';
if fh.t(end) > 120, t_div = 60; t_unit = 'min'; end
% Path in metres, or kilometres when it is long
x_div = 1; x_unit = 'm';
if range(fh.x) > 2e3, x_div = 1e3; x_unit = 'km'; end

%% Points of application (body axes: x forward, z up, origin at the CG)
% As used by the moment equation of vtolto_forces_3GDL (r x F of each
% term): rotor 1 at (xw, ztw), rotor 2 at (-xtw, ztw), wing at (xw, zw),
% tail at (-xtw, ztw). Nacelle drag is drawn at its rotor and fuselage
% forces at (xfus, zfus).

g  = params.geo;
st.cg    = [0; 0];
st.wing  = [g.xw;   g.zw];
st.fus   = [g.xfus; g.zfus];
st.tail  = [-g.xtw; g.ztw];
st.prop1 = [g.xw;   g.ztw];
st.prop2 = [-g.xtw; g.ztw];

%% Forces in inertial axes at frame k: {label, station, colour, @(k) [Fx; Fz]}

wind_L = @(k) [-sin(S.gamma(k)); cos(S.gamma(k))];   % Lift direction
wind_D = @(k) -[cos(S.gamma(k)); sin(S.gamma(k))];   % Drag direction
rotor  = @(T, H, ang) T*[cos(ang); sin(ang)] + H*[-sin(ang); cos(ang)];
nac    = @(D, ang) -D*[cos(ang); sin(ang)];

forces = {
    'Weight',          'cg',    [0 0 0],       @(k) [0; -W]
    'Rotor 1 (T, H)',  'prop1', [0.80 0 0.80], @(k) rotor(S.T1(k), S.H1(k), S.theta(k) + S.epsilon1(k))
    'Rotor 2 (T, H)',  'prop2', [0.45 0 0.60], @(k) rotor(S.T2(k), S.H2(k), S.theta(k) + S.epsilon2(k))
    'Nacelle 1 drag',  'prop1', [0.85 0.45 0], @(k) nac(S.D_nac1(k), S.theta(k) + S.epsilon1(k) - S.alpha_nac1(k))
    'Nacelle 2 drag',  'prop2', [0.85 0.60 0.1], @(k) nac(S.D_nac2(k), S.theta(k) + S.epsilon2(k) - S.alpha_nac2(k))
    'Wing lift',       'wing',  [0 0.55 0],    @(k) S.L_wing(k)*wind_L(k)
    'Wing drag',       'wing',  [0.85 0 0],    @(k) S.D_wing(k)*wind_D(k)
    'Fuselage lift',   'fus',   [0.3 0.75 0.3], @(k) S.L_fus(k)*wind_L(k)
    'Fuselage drag',   'fus',   [0.95 0.35 0.35], @(k) S.D_fus(k)*wind_D(k)
    'Tail lift',       'tail',  [0.1 0.45 0.3], @(k) S.L_tail(k)*wind_L(k)
    'Tail drag',       'tail',  [0.6 0.15 0.15], @(k) S.D_tail(k)*wind_D(k)
    };
n_f = size(forces, 1);

% One arrow scale for the whole animation: the largest force is 1.2 m long
F_max = W;
for k = 1:n_frames
    for j = 1:n_f
        F_max = max(F_max, norm(forces{j,4}(k)));
    end
end
f_scale = 1.2/F_max;

% Moments: {label, station, colour, values}; arc sweep 270 deg at the largest |M|
moments = {
    'Wing aero moment',     'wing', [0 0.55 0],     S.MA_wing
    'Fuselage aero moment', 'fus',  [0.3 0.75 0.3], S.MA_fus
    'Tail aero moment',     'tail', [0.1 0.45 0.3], S.MA_tail
    'Net moment at CG',     'cg',   [0 0 0],        S.M
    };
M_max = max(cellfun(@(m) max(abs(m)), moments(:,4)));
m_scale = 1.5*pi/max(M_max, eps);
m_radius = [0.15 0.15 0.15 0.28];

%% Figure layout

animFig = figure('Name', ['Animation: ' opts.name], 'Color', 'w', 'Position', [60 60 1400 860]);
axBody = axes(animFig, 'Position', [0.05 0.33 0.58 0.62]);
axTraj = axes(animFig, 'Position', [0.70 0.33 0.27 0.62]);
axF    = axes(animFig, 'Position', [0.05 0.06 0.40 0.18]);
axM    = axes(animFig, 'Position', [0.55 0.06 0.42 0.18]);

% Aircraft view (CG fixed at the origin, inertial axes)
hold(axBody, 'on'); axis(axBody, 'equal'); grid(axBody, 'on'); box(axBody, 'on');
xlim(axBody, [-2.2 2.2]); ylim(axBody, [-1.8 1.8]);
xlabel(axBody, 'x [m]'); ylabel(axBody, 'h [m]');
plot(axBody, [-2.2 2.2], [0 0], 'k:', 'HandleVisibility', 'off');   % Horizon

t_body = hgtransform('Parent', axBody);
[t_nac1, t_nac2] = draw_aircraft(t_body, params, st);

h_force = gobjects(n_f, 1);
for j = 1:n_f
    h_force(j) = quiver(axBody, 0, 0, 0, 0, 0, 'LineWidth', 2, 'MaxHeadSize', 0.6, ...
                        'Color', forces{j,3}, 'DisplayName', forces{j,1});
end
h_moment = gobjects(size(moments, 1), 2);
for j = 1:size(moments, 1)
    h_moment(j,1) = plot(axBody, NaN, NaN, '-', 'LineWidth', 1.8, 'Color', moments{j,3}, ...
                         'DisplayName', moments{j,1});
    h_moment(j,2) = patch(axBody, NaN, NaN, moments{j,3}, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
h_moment(4,1).LineStyle = '--';
legend(axBody, [h_force; h_moment(:,1)], 'Location', 'eastoutside', 'FontSize', 8);
h_title = title(axBody, '');

% Trajectory
% Trajectory (x in m or km); equal axes only when x and h are comparable
xs = fh.x/x_div;
hold(axTraj, 'on'); grid(axTraj, 'on'); box(axTraj, 'on');
plot(axTraj, xs, fh.h, ':', 'Color', [0.6 0.6 0.6]);
h_path  = plot(axTraj, NaN, NaN, 'b-', 'LineWidth', 1.5);
h_glyph = plot(axTraj, NaN, NaN, 'k-', 'LineWidth', 3);
h_pos   = plot(axTraj, NaN, NaN, 'ro', 'MarkerFaceColor', 'r');
if x_div == 1
    axis(axTraj, 'equal');
    x_pad = max(3, 0.15*range(fh.h));
    xlim(axTraj, [min(xs) - x_pad, max(xs) + x_pad]);
    ylim(axTraj, [min(fh.h) - 1, max(fh.h) + 2]);
    glyph_x = 0.08*max(range(fh.h), 10);
    glyph_h = glyph_x;
else
    x_pad = 0.05*range(xs);
    h_pad = max(10, 0.15*range(fh.h));
    xlim(axTraj, [min(xs) - x_pad, max(xs) + x_pad]);
    ylim(axTraj, [min(fh.h) - h_pad, max(fh.h) + h_pad]);
    glyph_x = 0.06*range(xs);                     % Glyph drawn in axis units, so its
    glyph_h = 0.06*diff(ylim(axTraj));            % angle reads as theta on screen
end
xlabel(axTraj, ['x [' x_unit ']']); ylabel(axTraj, 'h [m]'); title(axTraj, 'Trajectory');

% Time histories (at the nodes)
ts = fh.t/t_div;
hold(axF, 'on'); grid(axF, 'on'); box(axF, 'on');
plot(axF, ts, fh.T1, 'Color', forces{2,3}, 'LineWidth', 1.5);
plot(axF, ts, fh.T2, 'Color', forces{3,3}, 'LineWidth', 1.5);
plot(axF, ts, fh.L_aero, 'Color', [0 0.55 0], 'LineWidth', 1.5);
plot(axF, ts, fh.D_aero + fh.D_nac1 + fh.D_nac2, 'Color', [0.85 0 0], 'LineWidth', 1.5);
yline(axF, W, 'k--');
legend(axF, 'T_1', 'T_2', 'Aero lift', 'Aero + nacelle drag', 'Weight', 'Location', 'eastoutside', 'FontSize', 8);
xlabel(axF, ['t [' t_unit ']']); ylabel(axF, 'Force [N]');
cur_F = xline(axF, ts(1), 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');

hold(axM, 'on'); grid(axM, 'on'); box(axM, 'on');
plot(axM, ts, fh.M_Eng1 + fh.M_Eng2, 'Color', forces{2,3}, 'LineWidth', 1.5);
plot(axM, ts, fh.MA + fh.M_Wing + fh.M_Fus + fh.M_Tail, 'Color', [0 0.55 0], 'LineWidth', 1.5);
plot(axM, ts, fh.M_Nac1 + fh.M_Nac2, 'Color', forces{4,3}, 'LineWidth', 1.5);
plot(axM, ts, fh.M, 'k--', 'LineWidth', 1.5);
legend(axM, 'Rotors', 'Aero', 'Nacelles', 'Net (I_{yy} dq/dt)', 'Location', 'eastoutside', 'FontSize', 8);
xlabel(axM, ['t [' t_unit ']']); ylabel(axM, 'Pitching moment [N m]');
cur_M = xline(axM, ts(1), 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');

%% Video

if opts.save
    targetFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Images');
    if ~isfolder(targetFolder), mkdir(targetFolder); end
    video = VideoWriter(fullfile(targetFolder, opts.name), 'MPEG-4');
    video.FrameRate = opts.fps;
    open(video);
end
pause_frames = opts.play && strcmp(animFig.Visible, 'on');

%% Frames

for k = 1:n_frames
    th = S.theta(k);
    R  = [cos(th), -sin(th); sin(th), cos(th)];

    set(t_body, 'Matrix', makehgtform('zrotate', th));
    set(t_nac1, 'Matrix', makehgtform('translate', [st.prop1; 0]) * makehgtform('zrotate', S.epsilon1(k)));
    set(t_nac2, 'Matrix', makehgtform('translate', [st.prop2; 0]) * makehgtform('zrotate', S.epsilon2(k)));

    for j = 1:n_f
        p  = R*st.(forces{j,2});
        Fv = forces{j,4}(k)*f_scale;
        set(h_force(j), 'XData', p(1), 'YData', p(2), 'UData', Fv(1), 'VData', Fv(2));
    end
    for j = 1:size(moments, 1)
        p = R*st.(moments{j,2});
        [arc, head] = moment_arc(p, m_radius(j), moments{j,4}(k)*m_scale);
        set(h_moment(j,1), 'XData', arc(1,:),  'YData', arc(2,:));
        set(h_moment(j,2), 'XData', head(1,:), 'YData', head(2,:));
    end

    set(h_title, 'String', sprintf(['t = %.2f %s    V = %.1f m/s    h = %.1f m    ' ...
        '\\theta = %.1f\\circ    \\epsilon_1 = %.0f\\circ    \\epsilon_2 = %.0f\\circ'], ...
        tq(k)/t_div, t_unit, S.V(k), S.h(k), rad2deg(th), rad2deg(S.epsilon1(k)), rad2deg(S.epsilon2(k))));

    xk = S.x/x_div;
    set(h_path,  'XData', xk(1:k), 'YData', S.h(1:k));
    set(h_pos,   'XData', xk(k),   'YData', S.h(k));
    set(h_glyph, 'XData', xk(k) + glyph_x*[-1 1]*cos(th)/2, 'YData', S.h(k) + glyph_h*[-1 1]*sin(th)/2);
    cur_F.Value = tq(k)/t_div;
    cur_M.Value = tq(k)/t_div;

    drawnow;
    if opts.save, writeVideo(video, getframe(animFig)); end
    if pause_frames, pause(opts.speed/opts.fps); end
end

if opts.save
    close(video);
    fprintf('Animation saved to %s.mp4\n', fullfile(targetFolder, opts.name));
end

end

function [t_nac1, t_nac2] = draw_aircraft(t_body, params, st)
% Fuselage, wing and tail plates on t_body; nacelles with rotor disks on
% child transforms, placed and tilted every frame
L = params.fus_Lref;
c = params.geo.c_wing;
D = params.prop.diameter;
nac_w = 0.35; nac_h = 0.08;
grey = [0.55 0.55 0.55];

try
    [img_fuse, ~, a_fuse] = imread('fuselage_sideview.png');
    [img_nac,  ~, a_nac]  = imread('nacelle_sideview.png');
    use_img = true;
catch
    use_img = false;
end

if use_img
    image('Parent', t_body, 'XData', [-L/2 L/2], 'YData', [0.2 -0.2], 'CData', img_fuse, 'AlphaData', a_fuse);
else
    patch('Parent', t_body, 'XData', [-L/2 L/2 L/2 -L/2], 'YData', [-0.1 -0.1 0.1 0.1], ...
          'FaceColor', grey, 'EdgeColor', 'k', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');
end
for p = {st.wing, st.tail}
    patch('Parent', t_body, 'XData', p{1}(1) + [-c/2 c/2 c/2 -c/2], 'YData', p{1}(2) + [-0.015 -0.015 0.015 0.015], ...
          'FaceColor', grey*0.7, 'EdgeColor', 'k', 'HandleVisibility', 'off');
end

t_nac1 = hgtransform('Parent', t_body);
t_nac2 = hgtransform('Parent', t_body);
for t_nac = [t_nac1, t_nac2]
    if use_img
        image('Parent', t_nac, 'XData', [-nac_w/2 nac_w/2], 'YData', [nac_h/2 -nac_h/2], 'CData', img_nac, 'AlphaData', a_nac);
    else
        patch('Parent', t_nac, 'XData', [-nac_w/2 nac_w/2 nac_w/2 -nac_w/2], 'YData', [-nac_h/2 -nac_h/2 nac_h/2 nac_h/2], ...
              'FaceColor', grey*0.8, 'EdgeColor', 'k', 'FaceAlpha', 0.8, 'HandleVisibility', 'off');
    end
    % Rotor disk, normal to the thrust axis at the nacelle front
    line('Parent', t_nac, 'XData', [nac_w/2 nac_w/2], 'YData', [-D/2 D/2], ...
         'Color', [0.3 0.3 0.3], 'LineWidth', 2.5, 'HandleVisibility', 'off');
end
end

function [arc, head] = moment_arc(p, r, sweep)
% Arc of angle sweep (positive = nose up, counter-clockwise) around p,
% starting at the top, with a triangular arrowhead at its end
ang = pi/2 + linspace(0, sweep, 40);
arc = p + r*[cos(ang); sin(ang)];
if abs(sweep) < 1e-3
    head = [NaN; NaN];
    return
end
tip  = arc(:,end);
tang = sign(sweep)*[-sin(ang(end)); cos(ang(end))];   % Direction of travel
nrm  = [-tang(2); tang(1)];
s = 0.06;
head = [tip, tip - s*tang + 0.5*s*nrm, tip - s*tang - 0.5*s*nrm];
end
