clc; clear; close all;

t = linspace(0, 2*pi, 1000);
cx = 1; cy = 3; r = 5;
pts = [cx + r*cos(t); cy + r*sin(t)];

theta = 55 * pi/180;
Rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];

m = 2;
Refl = (1/(1+m^2)) * [1-m^2, 2*m; 2*m, m^2-1];

S = [2 0; 0 1];

rot_pts   = Rot  * pts;
refl_pts  = Refl * pts;
str_pts   = S    * pts;

lx = linspace(-10, 10, 100);
ly = m * lx;

figure('Position', [100 100 1200 900]);

subplot(2,2,1);
plot(pts(1,:), pts(2,:), 'b', 'LineWidth', 1.5);
axis equal; grid on;
title('Original Circle'); xlabel('x'); ylabel('y');

subplot(2,2,2);
plot(pts(1,:), pts(2,:), 'b--', 'LineWidth', 1);
hold on;
plot(rot_pts(1,:), rot_pts(2,:), 'r', 'LineWidth', 1.5);
axis equal; grid on;
legend('Original', 'Rotated 55°');
title('Rotation about Origin (55°)'); xlabel('x'); ylabel('y');

subplot(2,2,3);
plot(pts(1,:), pts(2,:), 'b--', 'LineWidth', 1);
hold on;
plot(refl_pts(1,:), refl_pts(2,:), 'g', 'LineWidth', 1.5);
plot(lx, ly, 'k--', 'LineWidth', 1);
axis equal; grid on; xlim([-10 10]); ylim([-10 10]);
legend('Original', 'Reflected', 'y=2x');
title('Reflection about y=2x'); xlabel('x'); ylabel('y');

subplot(2,2,4);
plot(pts(1,:), pts(2,:), 'b--', 'LineWidth', 1);
hold on;
plot(str_pts(1,:), str_pts(2,:), 'm', 'LineWidth', 1.5);
axis equal; grid on;
legend('Original', 'Stretched');
title('Stretch 2x along x-axis'); xlabel('x'); ylabel('y');

saveas(gcf, 'transformations.png');
close all;

fprintf('Q1 done\n');
