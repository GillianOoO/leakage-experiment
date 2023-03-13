clear all;

start = 0.00001;
final = 0.001;
step = start;

size_eigs = length(start:step:final);

num_eigs = zeros(1, size_eigs); %not equal 0 and 1.

eigs = zeros(9, size_eigs);

k = 1;
for eps = start:step:final
    eigs(:, k) = eigs_iswap(eps);
    for t = 1 : 9
        if abs(eigs(t,k))<1e-15
            eigs(t,k) = 0;
        elseif abs(1 - eigs(t,k)) > 1e-5
            num_eigs(k) = num_eigs(k) + 1;
        end
    end
    k = k + 1;
end

eps = start:step:final;
%plot(eps, num_eigs, 'o');

% semilogx(eps, num_eigs, 'o');
ya = eigs(2,:);
yb = eigs(5,:);
% yc = eigs(5,:);
x = 1:length(ya);

intval = 5;
pa = polyfit(eps(1:intval), ya(1:intval),1);
pb = polyfit(eps(1:intval), yb(1:intval),1);
% pc = polyfit(eps(1:intval), yc(1:intval),1);

fa = polyval(pa, eps);
fb = polyval(pb, eps);
% fc = polyval(pc, eps);

% fb = fit(eps, yb, 'poly1');
% fc = fit(eps, yc, 'poly1');
% plot(eps, ya,'o', eps, yb,'x', eps, yc, '*');
figure
box on
hold on
 plot(eps, ya,'p', 'color', rand(1,3), 'MarkerSize',8);
%curve(1) = plot(eps, fa,  'color', rand(1,3), 'Linewidth', 1);
plot(eps, yb,'o', 'color', rand(1,3), 'MarkerSize',8);
%curve(2) = plot(eps, fb,  'color', rand(1,3), 'Linewidth', 1);
% plot(eps, yc,'d', 'color', rand(1,3), 'MarkerSize',8);
% curve(3) = plot(eps, fc, 'color', rand(1,3), 'Linewidth', 1);
hold off
%str1 = num2str(pa(1));
% str1 = strcat(str1, '\epsilon');
% str1 = strcat(str1, '+');
% str1 = strcat(str1, num2str(pa(2)));
% 
% % legend(curve2, str1);
% 
% str2 = num2str(pb(1));
% str2 = strcat(str2, '\epsilon');
% str2 = strcat(str2, '+');
% str2 = strcat(str2, num2str(pb(2)));

% legend(curve4, str2);

% str3 = num2str(pc(1));
% str3 = strcat(str3, '\epsilon');
% str3 = strcat(str3, '+');
% str3 = strcat(str3, num2str(pc(2)));

% legend(curve6, str3);

%legend(curve, {str1, str2}, 'FontSize',12, 'Location','SouthWest');



