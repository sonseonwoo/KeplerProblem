% =========================================================================
% Title    : ODE-Based Orbit Propagation (Two-Body Problem)
% Author   : Seonwoo Son (sonseonwoo)
% Source   : Vallado, *Fundamentals of Astrodynamics and Applications*, 5th Edition
% Date     : 2022-08-29
% Purpose  : 
%   - Propagate orbital motion using Newton's equation of motion
%   - Solves the two-body problem via ODE integration (e.g., ode45)
% =========================================================================

r0 = [9818.93; 10324.744; 9704.202];  % km
v0 = [3.717227; 4.217780; -1.898238]; % km/s
T_orbit = 68339.0582;                 % example orbital period (sec)

[t, r, v] = ode_2body_propagation(r0, v0, [0 1*T_orbit]);

plot3(r(:,1), r(:,2), r(:,3), 'g')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('ODE-Based Orbit Propagation');
grid on; axis equal;


function [t_out, r_out, v_out] = ode_2body_propagation(r0, v0, tspan, mu, reltol)
% Inputs:
%   r0     - Initial position vector [3x1] in km
%   v0     - Initial velocity vector [3x1] in km/s
%   tspan  - Time span for propagation [t0 tf] in seconds
%   mu     - Gravitational parameter (default: 398600.4415 km^3/s^2)
%   reltol - ODE solver relative tolerance (default: 1e-8)
%
% Outputs:
%   t_out  - Time vector from ODE solver
%   r_out  - Propagated position vectors [Nx3]
%   v_out  - Propagated velocity vectors [Nx3]

    if nargin < 4 || isempty(mu)
        mu = 398600.4415; % Earth's gravitational parameter [km^3/s^2]
    end
    if nargin < 5 || isempty(reltol)
        reltol = 1e-8;
    end

    y0 = [r0(:); v0(:)]; % 6x1 initial state vector
    options = odeset('RelTol', reltol);

    % ODE integration using ode45
    [t_out, y_out] = ode45(@(t, y) two_body_equation(t, y, mu), tspan, y0, options);

    % Output separation
    r_out = y_out(:, 1:3);  % position
    v_out = y_out(:, 4:6);  % velocity
end

function dydt = two_body_equation(~, y, mu)
% TWO_BODY_EQUATION
% Computes the derivative of the state vector under two-body gravity.

    r = y(1:3);
    v = y(4:6);

    r_norm = norm(r);
    a = -mu / r_norm^3 * r;

    dydt = [v; a];
end