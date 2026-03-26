clc; clear; close all;

t = linspace(0, 2*pi, 1000);
cx = 1; cy = 3; r = 5;

orig = [cx + r*cos(t); cy + r*sin(t); zeros(1, length(t))];

n = [2; 3; -1];
nn = dot(n, n);

refl = zeros(size(orig));
for i = 1:size(orig, 2)
    p = orig(:, i);
    refl(:, i) = p - 2 * (dot(n, p) / nn) * n;
end

[xx, yy] = meshgrid(linspace(-8, 10, 30), linspace(-4, 12, 30));
zz = 2*xx + 3*yy;

figure('Position', [100 100 900 700]);
surf(xx, yy, zz, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2, 'FaceColor', [0.5 0.7 0.9]);
hold on;
plot3(orig(1,:), orig(2,:), orig(3,:), 'b', 'LineWidth', 2);
plot3(refl(1,:), refl(2,:), refl(3,:), 'r', 'LineWidth', 2);

for k = 1:50:size(orig,2)
    plot3([orig(1,k) refl(1,k)], [orig(2,k) refl(2,k)], [orig(3,k) refl(3,k)], ...
        'k--', 'LineWidth', 0.5);
end

xlabel('x'); ylabel('y'); zlabel('z');
title('Reflection of Circle about Plane z=2x+3y');
legend('Plane z=2x+3y', 'Original Circle', 'Reflected Circle', 'Location', 'best');
grid on; view(35, 25);

saveas(gcf, 'reflection_3d.png');
close all;

fprintf('Q2 done\n');
