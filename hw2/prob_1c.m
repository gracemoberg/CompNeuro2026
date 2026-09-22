% 1c: active front speed c vs threshold a at D = 1
 
D = 1;

a = linspace(0.001, 0.999, 1000);
c = sqrt(D) .* (1 - 2*a) ./ sqrt(a .* (1 - a));
c_passive = 2 * sqrt(D);
 
% active speed = passive speed
a_cross = (2 - sqrt(2)) / 4;

figure;
plot(a, c, 'LineWidth', 2);
hold on;
 
yline(c_passive, '--', 'passive speed = 2\sqrt{D}', 'LineWidth', 1.5);
 
xline(a_cross, ':', sprintf('a = %.4f', a_cross), 'LineWidth', 1.5);
 
% crossing point
plot(a_cross, c_passive, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'auto');
% stall threshold
xline(0.5, '--', 'stall: a = 0.5', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.8);
xlabel('Threshold a');
ylabel('Wave speed c');
title('Traveling front speed c vs threshold a');
legend('c', 'c_{pass}', 'active/passive crossing a value', 'active/passive cross point', 'front stalls here', 'Location', 'best');
 
grid on;
hold off;
