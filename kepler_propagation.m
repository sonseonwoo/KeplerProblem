% =========================================================================
% Title    : Orbital Propagation via COE (Kepler's Equation)
% Author   : Seonwoo Son (sonseonwoo)
% Source   : Vallado, *Fundamentals of Astrodynamics and Applications*, 5th Edition
% Date     : 2023-08-23
% Purpose  : 
%   - Convert RV to COE using RV2COE.m
%   - Propagate orbit using Kepler's Equation (Newton-Raphson)
%   - Convert COE back to RV using COE2RV.m
% =========================================================================

format long; clc; clear;
addpath(genpath(pwd));   % Include modules (RV2COE.m, COE2RV.m)

%% ---------------------- INITIAL VECTORS (ECI) ---------------------------
r0 = [9818.93 10324.744 9704.202];     % [km]
v0 = [3.717227 4.217780 -1.898238];    % [km/s]
mu = 3.986004418e5;                   % [km^3/s^2]

%% ---------------------- RV to COE --------------------------------------
coe0 = RV2COE([r0 v0]);
a = coe0(1); e = coe0(2); i = coe0(3); RAAN = coe0(4); AoP = coe0(5); TA0 = coe0(6);

%% ---------------------- TIME SETTINGS ----------------------------------
delta_T = 1:10:60000;          % [s]
n_mean = sqrt(mu / a^3);      % [rad/s]

%% ---------------------- INITIAL MEAN ANOMALY ---------------------------
E0 = 2 * atand(sqrt((1-e)/(1+e)) * tand(TA0/2));   % [deg]
M0 = deg2rad(E0) - e * sind(E0);                  % [rad]
M_vec = mod(M0 + n_mean * delta_T, 2*pi);         % [rad]

%% ---------------------- SOLVER SELECTION -------------------------------
% Choose solver: 'newton', 'bisection', 'regula-falsi'
solver_type = 'newton';

E_vec = zeros(size(M_vec));
switch lower(solver_type)
    case 'newton'
        fE  = @(E, M) E - e*sin(E) - M;
        dfE = @(E) 1 - e*cos(E);
        for k = 1:length(M_vec)
            E_vec(k) = newtraph(@(E) fE(E, M_vec(k)), dfE, pi/2);
        end

    case 'bisection'
        for k = 1:length(M_vec)
            fE = @(E) E - e*sin(E) - M_vec(k);
            E_vec(k) = bisect(fE, 0, 2*pi);
        end

    case 'regula-falsi'
        for k = 1:length(M_vec)
            fE = @(E) E - e*sin(E) - M_vec(k);
            E_vec(k) = lin_interp(fE, 0, 2*pi);
        end

    otherwise
        error('Invalid solver selected. Choose newton, bisection, or regula-falsi.');
end



%% ---------------------- TRUE ANOMALY & COE PROPAGATION ------------------
TA_vec = 2 * atan(sqrt((1+e)/(1-e)) .* tan(E_vec/2)); % [rad]

%% ---------------------- BACK TO RV USING COE2RV -------------------------
r_IJK = zeros(3, length(TA_vec));
v_IJK = zeros(3, length(TA_vec));

for k = 1:length(TA_vec)
    nu_deg_k = rad2deg(TA_vec(k));
    [r_IJK(:,k), v_IJK(:,k)] = COE2RV(a, e, i, RAAN, AoP, nu_deg_k, mu);
end
%% ---------------------- PLOTTING RESULTS -------------------------------
figure;
plot3(r_IJK(1,:), r_IJK(2,:), r_IJK(3,:), 'b');
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Orbital Trajectory in ECI Frame');
grid on; axis equal;

%% Optional sanity check
% fprintf('Initial h = %.6f km^2/s, Propagated h = %.6f km^2/s\n', h_norm, norm(cross(r_IJK(:,1), v_IJK(:,1))));
