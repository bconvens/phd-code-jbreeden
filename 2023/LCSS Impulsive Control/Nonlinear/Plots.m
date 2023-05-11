function Plots
sim_lin15 = load('Runs/linear15.mat');
sim_lin30 = load('Runs/linear30.mat');
sim_lin45 = load('Runs/linear45.mat');
sim_lin60 = load('Runs/linear60.mat');
sim_nl15 = load('Runs/nonlin15.mat');
sim_nl30 = load('Runs/nonlin30.mat');
sim_nl45 = load('Runs/nonlin45.mat');
sim_nl60 = load('Runs/nonlin60.mat');
sim_nl90 = load('Runs/nonlin90.mat');
sim_nl120 = load('Runs/nonlin120.mat');
sim_nl180 = load('Runs/nonlin180.mat');
sim_nl240 = load('Runs/nonlin240.mat');
sim_nl300 = load('Runs/nonlin300.mat');
sim_nl360 = load('Runs/nonlin360.mat');
sim_nl420 = load('Runs/nonlin420.mat');
sim_nl480 = load('Runs/nonlin480.mat');
sim_nl540 = load('Runs/nonlin540.mat');

figure(11); clf;
theta = linspace(0, 2*pi, 100);
cx = cos(theta);
cy = sin(theta);
rho = 1;
xc = get_center(0);
obs = obstacle_location(1,0)-xc(1:2); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r'); hold on;
obs = obstacle_location(2,0)-xc(1:2); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(3,0)-xc(1:2); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(4,0)-xc(1:2); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
fill([-10 -10 10 10 -10], [6 0 0 6 6], 'r');
xlabel 'x_1 (km)'; ylabel 'x_2 (km)';

color1 = [0; 0.7; 0];
color2 = [0.2; 0.2; 1];
color3 = [0.9; 0; 0.9];
% p1 = plot(sim_lin15.x_lin(1,:)/1e3, sim_lin15.x_lin(2,:)/1e3, 'b', 'LineWidth', 2)
p2 = plot(sim_lin30.x_lin(1,:)/1e3, sim_lin30.x_lin(2,:)/1e3, '-', 'LineWidth', 2, 'Color', color1);
p3 = plot(sim_lin45.x_lin(1,:)/1e3, sim_lin45.x_lin(2,:)/1e3, '--', 'LineWidth', 2, 'Color', color1);
p4 = plot(sim_lin60.x_lin(1,:)/1e3, sim_lin60.x_lin(2,:)/1e3, ':', 'LineWidth', 3, 'Color', color1);
% p5 = plot(sim_nl15.x_lin(1,:)/1e3, sim_nl15.x_lin(2,:)/1e3, 'g', 'LineWidth', 2)
p6 = plot(sim_nl30.x_lin(1,:)/1e3, sim_nl30.x_lin(2,:)/1e3, '-', 'LineWidth', 2, 'Color', color2);
p7 = plot(sim_nl45.x_lin(1,:)/1e3, sim_nl45.x_lin(2,:)/1e3, '--', 'LineWidth', 2, 'Color', color2);
p8 = plot(sim_nl60.x_lin(1,:)/1e3, sim_nl60.x_lin(2,:)/1e3, ':', 'LineWidth', 3, 'Color', color2);
% p9 = plot(sim_nl90.x_lin(1,:)/1e3, sim_nl90.x_lin(2,:)/1e3, 'g', 'LineWidth', 2)
% p10 = plot(sim_nl120.x_lin(1,:)/1e3, sim_nl120.x_lin(2,:)/1e3, 'r', 'LineWidth', 2)
p11 = plot(sim_nl180.x_lin(1,:)/1e3, sim_nl180.x_lin(2,:)/1e3, '-', 'LineWidth', 2, 'Color', color3);
% p12 = plot(sim_nl240.x_lin(1,:)/1e3, sim_nl240.x_lin(2,:)/1e3, 'c', 'LineWidth', 2)
p13 = plot(sim_nl300.x_lin(1,:)/1e3, sim_nl300.x_lin(2,:)/1e3, '--', 'LineWidth', 2, 'Color', color3);
% p14 = plot(sim_nl360.x_lin(1,:)/1e3, sim_nl360.x_lin(2,:)/1e3, 'g', 'LineWidth', 2)
p15 = plot(sim_nl420.x_lin(1,:)/1e3, sim_nl420.x_lin(2,:)/1e3, ':', 'LineWidth', 3, 'Color', color3);
% p16 = plot(sim_nl480.x_lin(1,:)/1e3, sim_nl480.x_lin(2,:)/1e3, 'g', 'LineWidth', 2)
% p17 = plot(sim_nl540.x_lin(1,:)/1e3, sim_nl540.x_lin(2,:)/1e3, 'g', 'LineWidth', 2)

plot(sim_lin15.x0(1)/1e3, sim_lin15.x0(2)/1e3, 'ko', 'MarkerFaceColor', 'k');
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');

ax = gca;
legend(ax,[p2 p3 p4], {'$$\psi_h$$, 30 sec', '$$\psi_h$$, 45 sec', '$$\psi_h$$, 60 sec'}, ...
    'FontSize', 13, 'Location', 'SouthEast', 'interpreter', 'latex');
ax2 = axes('position',get(gca,'position'),'visible','off');
l2 = legend(ax2,[p6 p7 p8 p11 p13 p15], {'$$\psi_h^*$$, 30 sec', '$$\psi_h^*$$, 45 sec', '$$\psi_h^*$$, 60 sec', ...
    '$$\psi_h^*$$, 180 sec', '$$\psi_h^*$$, 300 sec', '$$\psi_h^*$$, 420 sec'}, ...
    'FontSize', 13, 'Location', 'SouthWest', 'interpreter', 'latex');

set(l2, 'Position', [0.205357142857143 0.129365079365079 0.240648319691257 0.312619041488284]);
axis(ax,'equal');
axis(ax,[-8 8 -14 1]);
axis(ax2,'equal');
axis(ax2,[-8 8 -14 1]);

set(gcf, 'Position', [1972 271 560 420]);
text(0.25, -10.04, '$$x_0$$', 'interpreter', 'latex', 'FontSize', 14)
text(0.1, 0.55, '$$x_{target}$$', 'interpreter', 'latex', 'FontSize', 14)

%%
figure(12); clf;
i_end = 2000;
p1 = plot(sim_lin45.t(1:i_end)/60, sim_lin45.u(1,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [0;0.7;0], 'Color', color1); hold on;
p2 = plot(sim_nl45.t(1:i_end)/60, sim_nl45.u(1,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [0;0;1], 'Color', color2);
p3 = plot(sim_nl180.t(1:i_end)/60, sim_nl180.u(1,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [1;0;1], 'Color', color3);
set(gca, 'Xlim', [0 30], 'YLim', [-12 30]);
set(gcf, 'Position', [2600 900 560 160]);
l = legend([p1 p2 p3], {'$$\psi_h$$, 45 sec', '$$\psi_h^*$$, 45 sec', '$$\psi_h^*$$, 180 sec'}, ...
    'FontSize', 13, 'Location', 'NorthEast', 'interpreter', 'latex');
xlabel 'Time (minutes)';
ylabel 'u_1 (m/s)';
set(l, 'Position', [0.664113585070648 0.526666674713293 0.240648319691257 0.421249991953373]);


figure(13); clf;
plot(sim_lin45.t(1:i_end)/60, sim_lin45.u(2,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [0;0.7;0], 'Color', color1); hold on;
plot(sim_nl45.t(1:i_end)/60, sim_nl45.u(2,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [0;0;1], 'Color', color2);
plot(sim_nl180.t(1:i_end)/60, sim_nl180.u(2,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [1;0;1], 'Color', color3);
set(gca, 'Xlim', [0 30], 'YLim', [-30 15]);
set(gcf, 'Position', [2600 600 560 160]);
xlabel 'Time (minutes)';
ylabel 'u_2 (m/s)';

%%
figure(14); clf;
sim = sim_lin45;
p1a = plot(sim.t/60, sim.V, 'Color', color1); hold on;
p1b = plot(sim.t(sim.jumps)/60, sim.V(sim.jumps), '.', 'Color', color1, 'MarkerSize', 10);
sim = sim_nl45;
p2a = plot(sim.t/60, sim.V, 'Color', color2); hold on;
p2b = plot(sim.t(sim.jumps)/60, sim.V(sim.jumps), '.', 'Color', color2, 'MarkerSize', 10);
sim = sim_nl180;
p3a = plot(sim.t/60, sim.V, 'Color', color3); hold on;
p3b = plot(sim.t(sim.jumps)/60, sim.V(sim.jumps), '.', 'Color',  color3, 'MarkerSize', 10);
p4 = plot(-1, 1, 'k.', 'MarkerSize', 10);
set(gcf, 'Position', [2600 300 560 180]);
set(gca, 'Xlim', [0 30]);
xlabel 'Time (minutes)';
ylabel 'V';
legend([p1a p2a p3a p4], ...
    {'$$\psi_h$$, 45 sec, Flows', '$$\psi_h^*$$, 45 sec, Flows', ...
    '$$\psi_h^*$$, 18 sec, Flows', 'Jumps'}, ...
    'FontSize', 13, 'Location', 'NorthEast', 'interpreter', 'latex');

end