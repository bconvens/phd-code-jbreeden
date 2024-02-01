%% CBFs Slide
t = 0:.1:15;
N = length(t);
x1 = zeros(4, N);
u1 = zeros(2, N-1);
x2 = zeros(4, N);
u2 = zeros(2, N-1);

x1(:,1) = [-3; -0.2; 0; 0];
for i=1:N-1
    u1(:,i) = CalculateU_nom(t(i),x1(:,i));
    x1(:,i+1) = UpdateX(t(i+1)-t(i), x1(:,i), u1(:,i));
end

x2(:,1) = x1(:,1);
for i=1:N-1
    u2(:,i) = CalculateU(t(i),x2(:,i));
    x2(:,i+1) = UpdateX(t(i+1)-t(i), x2(:,i), u2(:,i));
end

% figure(2); clf;
% plot(t, x(1:2,:));

%%
f1 = figure(1); clf;
plot(x1(1,:), x1(2,:), '--', 'Color', [0.6 0.6 0.6]); hold on; axis equal;
axis([-4 4 -2 2]);
px1 = plot(x1(1,1), x2(2,1), 'ko', 'MarkerFaceColor', 'k');
set(gca, 'XTick', [], 'YTick', []);
title('Without Obstacle','FontSize', 16);
set(f1, 'Position', [1900 660 560 300], 'Color', 'w');

f2 = figure(2); clf;
cx = -1:.01:1;
cy = sqrt(1-cx.^2);
cx = [cx, fliplr(cx)];
cy = [cy, -cy];
fill(cx*.95, cy*.95, 'r', 'EdgeColor', [0.5 0 0]); hold on; axis equal;
axis([-4 4 -2 2]);
plot(x2(1,:), x2(2,:), '--', 'Color', [0.6 0.6 0.6]);
px2 = plot(x2(1,1), x2(2,1), 'ko', 'MarkerFaceColor', 'k');
set(gca, 'XTick', [], 'YTick', []);
title('With Obstacle and CBF','FontSize', 16);
set(f2, 'Position', [2500 660 560 300], 'Color', 'w');

%%
record = 1;
moviename = 'Comparison';
frame_rate = 30;
vidObj = VideoWriter(moviename);
vidObj.FrameRate = frame_rate;
if record, open(vidObj); end
clear this_frame;
% this_frame.cdata = zeros(300, 1120, 3, 'uint8');
this_frame.cdata = zeros(750, 1400*2, 3, 'uint8');
this_frame.colormap = [];

for i=1:N
    set(px1, 'XData', x1(1,i), 'YData', x1(2,i));  
    set(px2, 'XData', x2(1,i), 'YData', x2(2,i));
    drawnow;
    if record
        frame1 = getframe(f1);
        frame2 = getframe(f2);
%         this_frame.cdata(:,1:560,:) = frame1.cdata;
%         this_frame.cdata(:,561:end,:) = frame2.cdata;
        this_frame.cdata(:,1:1400,:) = frame1.cdata;
        this_frame.cdata(:,1401:end,:) = frame2.cdata;
        writeVideo(vidObj,this_frame);
    end
end

if record, close(vidObj); end

function u = CalculateU_nom(t,x)
kp = 0.6;
kd = 1.2;
xdes = [3; -0.2];
dx = (xdes - x(1:2));
x_max = 2;
if norm(dx) > x_max, dx = dx/norm(dx)*x_max; end
u = kp*dx - kd*x(3:4);
end

function u = CalculateU(t,x)
u_nom = CalculateU_nom(t,x);

h = @(x) 1 - norm(x);
dhdx = @(x) -x'/norm(x);
hdot = @(x,y) dhdx(x)*y;
mu = 0.5;
H = @(x,y) h(x) + hdot(x,y)*abs(hdot(x,y))/(2*mu);
LfH = @(x,y) hdot(x,y) - abs(hdot(x,y))/mu*(x(1)*y(2) - x(2)*y(1))^2/norm(x)^3;
LgH = @(x) dhdx(x);
alpha = 1;

A = abs(hdot(x(1:2),x(3:4)))*LgH(x(1:2));
b = -alpha*H(x(1:2), x(3:4)) - LfH(x(1:2), x(3:4));
u = quadprog(eye(2), -u_nom, A, b, [], [], [], [], [], optimset('Display', 'off'));
end

function x2 = UpdateX(dt,x,u)
Ad = [eye(2), dt*eye(2); zeros(2), eye(2)];
Bd = [dt^2/2*eye(2); dt*eye(2)];
x2 = Ad*x + Bd*u;
end