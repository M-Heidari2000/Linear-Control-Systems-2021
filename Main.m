%% ------------------- Bicycle Dynamics and Control -------------------

%% Milad Heidari - 98101469

%% notes:
% run section1(or section2) to load parameters before running other sections 

%% Section1 - Physical parameters of the bicycle "with a rider"
% to consider the bike "with a rider" run this section otherwise run the
% next section
clear; clc; close all;

b = 1;  % wheel base (in meters)
c = 0.08;   % trail (in meters)
lambda = 70;    % head angle (in degrees)
Rrw = 0.35;    % rear wheel radius (in meters)
Rfw = 0.35;    % front wheel radius (in meters)
% The following values are the masses of the different parts of the bicycle
Mrf = 87;   % rear frame mass (kg)
Mff = 2;    % front frame mass (kg)
Mrw = 1.5;  % rear wheel mass (kg)
Mfw = 1.5;  % front wheel mass (kg)
m = Mrf + Mff + Mrw + Mfw;  % total mass of the system (kg)
% The follwing vectors are the center of mass coordinates in the form [x, z]
CMrf = [0.492, 1.028];  % CM of the rear frame
CMff = [0.866, 0.676];  % CM of the front frame
CMrw = [0, Rrw];    % CM of the rear wheel
CMfw = [b, Rfw];    % CM of the front wheel
% The following vectors are inertia tensor in the form [Jxx, Jxz, Jyy, Jzz]
ITrf = [3.28, -0.603, 3.88, 0.566];   % inertia tensor of the rear frame
ITff = [0.08, 0.02, 0.07, 0.02];    % inertia tensor of the front frame
ITrw = [0.07, 0, 0.14, 0.07];   % inertia tensor of the rear wheel
ITfw = [0.07,0, 0.14, 0.07];    % inertia tensor of the front frame
% Total center of mass:
CM = (Mrf*CMrf + Mff*CMff + Mrw*CMrw + Mfw*CMfw)/(Mrf + Mff + Mrw + Mfw);
a = CM(1);  % horizontal distance of center of mass from P1
h = CM(2);  % height of the center of mass

g = 9.81;   % gravity
lambda = lambda*pi/180;

%% Section2 - Physical parameters of the bicycle "without a rider"
clear; clc; close all;

b = 1;  % wheel base (in meters)
c = 0.08;   % trail (in meters)
lambda = 70;    % head angle (in degrees)
Rrw = 0.35;    % rear wheel radius (in meters)
Rfw = 0.35;    % front wheel radius (in meters)
% The following values are the masses of the different parts of the bicycle
Mrf = 12;   % rear frame mass (kg)
Mff = 2;    % front frame mass (kg)
Mrw = 1.5;  % rear wheel mass (kg)
Mfw = 1.5;  % front wheel mass (kg)
m = Mrf + Mff + Mrw + Mfw;  % total mass of the system (kg)
% The follwing vectors are the center of mass coordinates in the form [x, z]
CMrf = [0.439, 0.579];  % CM of the rear frame
CMff = [0.866, 0.676];  % CM of the front frame
CMrw = [0, Rrw];    % CM of the rear wheel
CMfw = [b, Rfw];    % CM of the front wheel
% The following vectors are inertia tensor in the form [Jxx, Jxz, Jyy, Jzz]
ITrf = [0.476, -0.274, 1.033, 0.527];   % inertia tensor of the rear frame
ITff = [0.08, 0.02, 0.07, 0.02];    % inertia tensor of the front frame
ITrw = [0.07, 0, 0.14, 0.07];   % inertia tensor of the rear wheel
ITfw = [0.07,0, 0.14, 0.07];    % inertia tensor of the front frame
% Total center of mass:
CM = (Mrf*CMrf + Mff*CMff + Mrw*CMrw + Mfw*CMfw)/(Mrf + Mff + Mrw + Mfw);
a = CM(1);  % horizontal distance of center of mass from P1
h = CM(2);  % height of the center of mass

g = 9.81;   % gravity
lambda = lambda*pi/180;

%% Section3 - Second-order model without a front fork - Stabilization
s = tf('s');
V = 1;  % velocity at the rear wheel (m/s)
Gpd = a*V*(s+V/a)/(b*h*(s^2-g/h));  % transfer function from steer angle to roll angle

% ------------------- stability without feedback ---------------------
figure(1);
ax = axes();
pzmap(Gpd)
title('Second-order model without a front fork', 'Interpreter', 'LaTeX', ...
       'FontSize', 15);
l_zero = findall(ax, 'tag', 'PZ_Zero');
l_pole = findall(ax, 'tag', 'PZ_Pole');
l_zero.MarkerSize = 10;
l_pole.MarkerSize = 10;
l_zero.LineWidth = 2;
l_pole.LineWidth = 2;
% not stable !

% ----------------- root locus with respect to velocity -----------------

k2 = 1;
V_vals = linspace(0.01, 10, 100);  % velocity at the rear wheel (m/s)
figure(2);
colormap jet
for i=1:length(V_vals)
    V = V_vals(i);
    CE = m*h^2*s^2+(m*a*h*V*k2/b)*s+((m*V^2*h*k2/b - m*g*h));
    CE_roots = zero(CE);
    scatter(real(CE_roots), imag(CE_roots), 20, [V, V], 'filled');
    hold on
end

colorbar
grid on
title('Root locus of characteristic equation', 'Interpreter', 'LaTeX' ...
       , 'FontSize', 15);
xlabel('Real(P)', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('Imag(P)', 'Interpreter', 'LaTeX', 'FontSize', 13);

%% Section4 - Second-order model with a front fork - Stabilization
s = tf('s');
T = 0;  % torque applied to handlebars (N.m)
Vc = sqrt(b*g*cot(lambda)); % critical velocity

% ------------------- stability without feedback ---------------------

% velocities less than critical velocity (Vc)
V_vals = linspace(0.01, Vc-0.01, 100);  % velocity at the rear wheel (m/s)
pole1_real = zeros(1, length(V_vals));
pole1_imag = zeros(1, length(V_vals));
pole2_real = zeros(1, length(V_vals));
pole2_imag = zeros(1, length(V_vals));
for i=1:length(V_vals)
    V = V_vals(i);
    K1 = b^2/(m*a*c*sin(lambda)*(V^2*sin(lambda)-b*g*cos(lambda)));
    K2 = b*g/(V^2*sin(lambda)-b*g*cos(lambda));
    Gpd = a*V*(s+V/a)/(b*h*(s^2-g/h));
    Gdt = K1/(1+K2*Gpd);    % transferfunction from steer torque to steer angle
    Gpt = Gpd * Gdt;    % transfer function from steer torque to roll angle
    poles = pole(minreal(Gpt));
    pole1_real(i) = real(poles(1));
    pole1_imag(i) = imag(poles(1));
    pole2_real(i) = real(poles(2));
    pole2_imag(i) = imag(poles(2));
end

figure(1);
subplot(2, 1, 1);
scatter(pole1_real, pole1_imag, 20 , V_vals, 'filled');
title('Trajectory of $P_1$', 'Interpreter', 'LaTeX', 'FontSize', 13);
xlabel('Real($P_1$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('Imag($P_1$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
colormap jet
colorbar
grid minor
hold on
ylim([-4, 4]);
subplot(2, 1, 2);
scatter(pole2_real, pole2_imag, 20, V_vals, 'filled');
title('Trajectory of $P_2$', 'Interpreter', 'LaTeX', 'FontSize', 13);
xlabel('Real($P_2$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('Imag($P_2$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
colormap jet
colorbar
grid minor
ylim([-4, 4]);

% velocities greater than critical velocity (Vc)
V_vals = linspace(Vc+0.01, 2*Vc, 1000);  % velocity at the rear wheel (m/s)
pole1_real = zeros(1, length(V_vals));
pole1_imag = zeros(1, length(V_vals));
pole2_real = zeros(1, length(V_vals));
pole2_imag = zeros(1, length(V_vals));
for i=1:length(V_vals)
    V = V_vals(i);
    K1 = b^2/(m*a*c*sin(lambda)*(V^2*sin(lambda)-b*g*cos(lambda)));
    K2 = b*g/(V^2*sin(lambda)-b*g*cos(lambda));
    Gpd = a*V*(s+V/a)/(b*h*(s^2-g/h));
    Gdt = K1/(1+K2*Gpd);    % transferfunction from steer torque to steer angle
    Gpt = Gpd * Gdt;    % transfer function from steer torque to roll angle
    poles = pole(minreal(Gpt));
    pole1_real(i) = real(poles(1));
    pole1_imag(i) = imag(poles(1));
    pole2_real(i) = real(poles(2));
    pole2_imag(i) = imag(poles(2));
end

figure(2);
subplot(2, 1, 1);
scatter(pole1_real, pole1_imag, 20, V_vals, 'filled');
title('Trajectory of $P_1$', 'Interpreter', 'LaTeX', 'FontSize', 13);
xlabel('Real($P_1$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('Imag($P_1$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
colormap jet
colorbar
grid minor
hold on
ylim([-4, 4]);
subplot(2, 1, 2);
scatter(pole2_real, pole2_imag, 20, V_vals, 'filled');
title('Trajectory of $P_2$', 'Interpreter', 'LaTeX', 'FontSize', 13);
xlabel('Real($P_2$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('Imag($P_2$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
colormap jet
colorbar
grid minor
ylim([-4, 4]);

%% Section5 - Second-order model with a front fork (Dynamic graphs) - Stabilization
s = tf('s');
T = 0;  % torque applied to handlebars (N.m)
Vc = sqrt(b*g*cot(lambda)); % critical velocity

% ------------------- stability without feedback ---------------------
V_vals = dataPoints(Vc, 0.3);  % velocity at the rear wheel (m/s)
figure(1);
title('Press any key to continue');
pause();
colormap jet
for i=1:length(V_vals)
    V = V_vals(i);
    K1 = b^2/(m*a*c*sin(lambda)*(V^2*sin(lambda)-b*g*cos(lambda)));
    K2 = b*g/(V^2*sin(lambda)-b*g*cos(lambda));
    Gpd = a*V*(s+V/a)/(b*h*(s^2-g/h));
    Gdt = K1/(1+K2*Gpd);    % transferfunction from steer torque to steer angle
    Gpt = Gpd * Gdt;    % transfer function from steer torque to roll angle
    poles = pole(minreal(Gpt));
    scatter(real(poles), imag(poles), 30, [V, V], 'filled');
    drawnow
    if i == 1
        axis([-10, 10, -4, 4]);
        set(gca, 'Color', 'black');
        set(gca, 'GridColor', 'white');
        set(gca, 'GridAlpha', 0.3);
        title('Pole trajectories', 'Interpreter', 'LaTeX', 'FontSize', 13);
        xlabel('Real($P$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
        ylabel('Imag($P$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
        grid on
    end
    hold on
    colorbar
    pause(0.01);
end


%% Section6 - Second order-model with rear-wheel steering - Stabilization
s = tf('s');
V = 1;  % velocity at the rear wheel (m/s)
V = -V; % reversing the sign of the velocity to account for rear-wheel steering
Gpd = a*V*(s+V/a)/(b*h*(s^2-g/h));  % transfer function from steer angle to roll angle

% ------------------- stability without feedback ---------------------
figure(1);
ax = axes();
pzmap(Gpd)
title('Second-order rear-wheel steering model', 'Interpreter', 'LaTeX', ...
       'FontSize', 15);
l_zero = findall(ax, 'tag', 'PZ_Zero');
l_pole = findall(ax, 'tag', 'PZ_Pole');
l_zero.MarkerSize = 10;
l_pole.MarkerSize = 10;
l_zero.LineWidth = 2;
l_pole.LineWidth = 2;
% not stable !

% ----------------- root locus with respect to velocity -----------------

k2 = 1;
V_vals = linspace(0.01, 10, 100);  % velocity at the rear wheel (m/s)
figure(2);  
colormap jet
for i=1:length(V_vals)
    V = V_vals(i); % reversing the sign of the velocity to account for rear-wheel steering
    V = -V;
    CE = m*h^2*s^2+(m*a*h*V*k2/b)*s+((m*V^2*h*k2/b - m*g*h));
    CE_roots = zero(CE);
    scatter(real(CE_roots), imag(CE_roots), 20, [-V, -V], 'filled');
    hold on
end

colorbar
grid on
title('Root locus of characteristic equation', 'Interpreter', 'LaTeX' ...
       , 'FontSize', 15);
xlabel('Real(P)', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('Imag(P)', 'Interpreter', 'LaTeX', 'FontSize', 13);



%% Section7 - Linear forth-order model "with a rider" - Stabilization
s = tf('s');
Vc = 1.74; % self-alignment velocity
T = 0;  % torque applied to handlebars(N.m)

% ------------------- stability without feedback ---------------------
V_vals = linspace(1, 20, 100);  % velocity at the rear wheel (m/s)
figure(1);
title('Press any key to continue');
pause();
colormap jet
for i=1:length(V_vals)
    V = V_vals(i);
    M = [96.8, -3.57; -3.57, 0.258];
    C = [0, -50.8; 0.436, 2.2];
    K0 = [-901, 35.17; 35.17, -12.03];
    K2 = [0, -87.06; 0, 3.5];
    outputs = inv(M.*s^2 + C.*V.*s + K0 + K2.*V^2);
    Gpt = outputs(1, 2);   % transfer function from steer torque to tilt angle
    Gdt = outputs(2, 2);   % transfer function from steer torque to steer angle
    poles = pole(minreal(Gpt));
    scatter(real(poles), imag(poles), 30, [V, V, V, V], 'filled');
    drawnow
    if i == 1
        axis([-40, 20, -30, 30]);
        set(gca, 'Color', 'black');
        set(gca, 'GridColor', 'white');
        set(gca, 'GridAlpha', 0.3);
        title('Pole trajectories', 'Interpreter', 'LaTeX', 'FontSize', 13);
        xlabel('Real($P$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
        ylabel('Imag($P$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
        grid on
    end
    hold on
    colorbar
    pause(0.01);
end

%% Section8 - Second-order model with a front fork - Maneuvering
s = tf('s');
T = 0;  % torque applied to handlebars (N.m)
Vc = sqrt(b*g*cot(lambda)); % critical velocity
V = Vc+4.5;  % velocity of the rear wheel (m/s)
K1 = b^2/(m*a*c*sin(lambda)*(V^2*sin(lambda)-b*g*cos(lambda)));
K2 = b*g/(V^2*sin(lambda)-b*g*cos(lambda));
Gpd = a*V*(s+V/a)/(b*h*(s^2-g/h));  % transfer function from steer angle to tilt angle
Gdt = K1/(1+K2*Gpd);    % transfer function from steer torque to steer angle
Ged = V^2/(b*s^2);  % transfer function from steer angle to path deviation
Get = Gdt * Ged;    % transfer function from torque to path deviation
Gpt = Gpd * Gdt;    % transfer function from torque to tilt angle

figure(1);
ax = axes();
pzmap(Get);
title(sprintf('$V= %0.2f$ $m/s$ ', V), 'Interpreter', 'LaTeX', 'FontSize', 14);
l_zero = findall(ax, 'tag', 'PZ_Zero');
l_pole = findall(ax, 'tag', 'PZ_Pole');
l_zero.MarkerSize = 10;
l_pole.MarkerSize = 10;
l_zero.LineWidth = 2;
l_pole.LineWidth = 2;

figure(2);
subplot(3, 1, 1)
[y, t] = step(Get);
plot(t, y, 'LineWidth', 1.5);
hold on
plot(t, zeros(1, length(t)), 'k-.');
xlabel('$t(s)$', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('$\eta(m)$', 'Interpreter', 'LaTeX', 'FontSize', 15);
title(sprintf('$V= %0.2f$ $m/s$ ', V), 'Interpreter', 'LaTeX', 'FontSize', 14);
grid minor
xlim([0, 1.5]);

subplot(3, 1, 2);
[y, t] = step(Gpt);
plot(t, y, 'LineWidth', 1.5, 'Color', 'magenta');
xlabel('$t(s)$', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('$\varphi(rad)$', 'Interpreter', 'LaTeX', 'FontSize', 15);
grid minor
hold on
plot(t, ones(1, length(t))*dcgain(Gpt), 'k-.');
xlim([0, 1.5]);

subplot(3, 1, 3);
[y, t] = step(Gdt);
plot(t, y, 'LineWidth', 1.5, 'Color', 'red');
xlabel('$t(s)$', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('$\delta(rad)$', 'Interpreter', 'LaTeX', 'FontSize', 15);
grid minor
hold on
plot(t, ones(1, length(t))*dcgain(Gdt), 'k-.');
xlim([0, 1.5]);

%% Section9 - Fourth-order model - Maneuvering
s = tf('s');
T = 0;  % torque applied to handlebars (N.m)
Vc = 5.96;  % critical velocity (m/s)
V = 7;  % velocity of the rear wheel (m/s)
M = [96.8, -3.57; -3.57, 0.258];
C = [0, -50.8; 0.436, 2.2];
K0 = [-901, 35.17; 35.17, -12.03];
K2 = [0, -87.06; 0, 3.5];outputs = inv(M.*s^2 + C.*V.*s + K0 + K2.*V^2);
Gpt = outputs(1, 2);   % transfer function from steer torque to tilt angle
Gdt = outputs(2, 2);   % transfer function from steer torque to steer angle
Ged = V^2/(b*s^2);  % transfer function from steer angle to path deviation
Get = Gdt * Ged;    % transfer function from torque to path deviation

figure(1);
ax = axes();
pzmap(Get);
title(sprintf('$V= %0.2f$ $m/s$ ', V), 'Interpreter', 'LaTeX', 'FontSize', 14);
l_zero = findall(ax, 'tag', 'PZ_Zero');
l_pole = findall(ax, 'tag', 'PZ_Pole');
l_zero.MarkerSize = 10;
l_pole.MarkerSize = 10;
l_zero.LineWidth = 2;
l_pole.LineWidth = 2;

figure(2);
subplot(3, 1, 1)
[y, t] = step(Get);
plot(t, y, 'LineWidth', 1.5);
hold on
plot(t, zeros(1, length(t)), 'k-.');
xlabel('$t(s)$', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('$\eta(m)$', 'Interpreter', 'LaTeX', 'FontSize', 15);
title(sprintf('$V= %0.2f$ $m/s$ ', V), 'Interpreter', 'LaTeX', 'FontSize', 14);
grid minor
xlim([0, 1.5]);

subplot(3, 1, 2);
[y, t] = step(Gpt);
plot(t, y, 'LineWidth', 1.5, 'Color', 'magenta');
xlabel('$t(s)$', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('$\varphi(rad)$', 'Interpreter', 'LaTeX', 'FontSize', 15);
grid minor
hold on
plot(t, ones(1, length(t))*dcgain(Gpt), 'k-.');
xlim([0, 10]);

subplot(3, 1, 3);
[y, t] = step(Gdt);
plot(t, y, 'LineWidth', 1.5, 'Color', 'red');
xlabel('$t(s)$', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('$\delta(rad)$', 'Interpreter', 'LaTeX', 'FontSize', 15);
grid minor
hold on
plot(t, ones(1, length(t))*dcgain(Gdt), 'k-.');
xlim([0, 10]);

%% Section10 - effects of rider lean (Dynamic graph)
s = tf('s');
T = 0;  % torque applied to handlebars (N.m)
Vc = sqrt(b*g*cot(lambda)); % critical velocity
mr = 75;    % mass of the rider
Jr = 2.804; % moment of inertia of the rider
hr = 0.3;   % distance between center of mass of rider and turning axis

Vc = sqrt(b*g*cot(lambda)); % critical velocity
V_vals = dataPoints(Vc, 0.3);  % velocity at the rear wheel (m/s)
figure(1);
title('Press any key to continue');
pause();
colormap jet
for i=1:length(V_vals)
    V = V_vals(i);  % velocity of the rear wheel (m/s)
    K1 = b^2/(m*a*c*sin(lambda)*(V^2*sin(lambda)-b*g*cos(lambda)));
    K2 = b*g/(V^2*sin(lambda)-b*g*cos(lambda));
    Gpp = (mr*g*hr-Jr*s^2)/(m*h^2*s^2 + m*a*h*K2*s/b + m*V^2*h*K2/b - m*g*h);
    poles = pole(minreal(Gpp));
    scatter(real(poles), imag(poles), 30, [V, V], 'filled');
    drawnow
    if i == 1
        axis([-20, 20, -10, 10]);
        set(gca, 'Color', 'black');
        set(gca, 'GridColor', 'white');
        set(gca, 'GridAlpha', 0.3);
        title('Pole trajectories', 'Interpreter', 'LaTeX', 'FontSize', 13);
        xlabel('Real($P$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
        ylabel('Imag($P$)', 'Interpreter', 'LaTeX', 'FontSize', 13);
        grid on
    end
    hold on
    colorbar
    pause(0.01);
end

%% 
s = tf('s');
T = 0;  % torque applied to handlebars (N.m)
Vc = sqrt(b*g*cot(lambda)); % critical velocity
mr = 75;    % mass of the rider
Jr = 2.804; % moment of inertia of the rider
hr = 0.3;   % distance between center of mass of rider and turning axis
V = 4;
K1 = b^2/(m*a*c*sin(lambda)*(V^2*sin(lambda)-b*g*cos(lambda)));
K2 = b*g/(V^2*sin(lambda)-b*g*cos(lambda));
Gpp = (mr*g*hr-Jr*s^2)/(m*h^2*s^2 + m*a*h*K2*s/b + m*V^2*h*K2/b - m*g*h);

figure(2);
[y, t] = step(0.1*Gpp);
plot(t, y, 'LineWidth', 1.5);
hold on
plot(t, ones(1, length(t))*dcgain(0.1*Gpp), 'k-.');
grid on
xlabel('$t(s)$', 'Interpreter', 'LaTeX', 'FontSize', 13);
ylabel('$\varphi(rad)$', 'Interpreter', 'LaTeX', 'FontSize', 15);
title('system response to $\phi(t)=0.1u(t)$', 'Interpreter', 'LaTeX', 'FontSize', 15);
axis([0, 25, -0.05, 0.15]);

