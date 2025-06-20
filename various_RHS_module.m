% =========================================================================
% Title    : ODE-Based Orbit Propagation (Two-Body Problem)
% Author   : Seonwoo Son (sonseonwoo)
% Source   : Vallado, *Fundamentals of Astrodynamics and Applications*, 5th Edition
% Date     : 2023-11-10
% Purpose  : 
%   - Propagate orbital motion using Newton's equation of motion
%   - Solves the two-body problem via ODE integration (e.g., ode45)
% =========================================================================

r0 = [3950.48; 4152.89; 3905.89];   % km  (400 km 고도 방향 유지)
v0 = [-5.2810; 5.0242; 0.0328];     % km/s (수직 방향 속도)

T_orbit = 2*pi*sqrt(norm(r0)^3/398600.4418);                 % example orbital period (sec)


%% 1) 2-Body
P1.mu = 398600.4418;
[t1,r1,v1] = ode_propagation(r0,v0,[0 T_orbit],@rhs_2body,P1);

%% 2) Exponential Drag
P2 = struct('mu',398600.4418,'Re',6378.1363,'Cd',2.2,'ApM',0.01,...
            'H',60.828,'h_ref',450,'we',[0;0;7.2921e-5]);
[t2,r2,v2] = ode_propagation(r0,v0,[0 T_orbit],@rhs_drag_exp,P2);


%% 3) NRLMSISE Drag
P3 = P2;                                   % 기본값 재사용
P3.datetime0 = datetime(2018,1,4,'TimeZone','UTC');

[t3,r3,v3] = ode_propagation(r0,v0,[0 T_orbit],@rhs_drag_nrlmsise,P3);

%% 4) J2 + Exponential Drag
P4 = struct('mu',398600.4418,'Re',6378.1363,'Cd',2.2,'ApM',0.01,...
            'H',60.828,'h_ref',450,'we',[0;0;7.2921e-5],'J2',1.0826269e-3);
[t4, r4, v4] = ode_propagation(r0, v0, [0 T_orbit], @rhs_J2_drag_exp, P4);


%% 5) J2 + NRLMSISE Drag
P5 = P4; % 기본값 재사용
P5.datetime0 = datetime(2018,1,4,'TimeZone','UTC');
[t5, r5, v5] = ode_propagation(r0, v0, [0 T_orbit], @rhs_J2_drag_nrlmsise, P5);


figure;

% ---- 1) 2-Body
subplot(2,3,1); 
plot3(r1(:,1), r1(:,2), r1(:,3), 'k'); grid on; axis equal
title('2-Body'); xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');

% ---- 2) Exp Drag
subplot(2,3,2); 
plot3(r2(:,1), r2(:,2), r2(:,3), 'r'); grid on; axis equal
title('Exponential Drag'); xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');

% ---- 3) NRLMSISE Drag
subplot(2,3,3); 
plot3(r3(:,1), r3(:,2), r3(:,3), 'b'); grid on; axis equal
title('NRLMSISE Drag'); xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');

% ---- 4) J2 + Exp Drag
subplot(2,3,4); 
plot3(r4(:,1), r4(:,2), r4(:,3), 'c'); grid on; axis equal
title('J2 + Exp Drag'); xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');

% ---- 5) J2 + NRLMSISE Drag
subplot(2,3,5); 
plot3(r5(:,1), r5(:,2), r5(:,3), 'm'); grid on; axis equal
title('J2 + NRLMSISE'); xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');

% ---- 6) 전체 비교
subplot(2,3,6); hold on; grid on; axis equal
plot3(r1(:,1), r1(:,2), r1(:,3), 'k');
plot3(r2(:,1), r2(:,2), r2(:,3), 'r');
plot3(r3(:,1), r3(:,2), r3(:,3), 'b');
plot3(r4(:,1), r4(:,2), r4(:,3), 'c');
plot3(r5(:,1), r5(:,2), r5(:,3), 'm');
legend('2-Body','Exp','NRLMSISE','J2+Exp','J2+NRLMSISE');
title('All Models Comparison');
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');




function [t_out, r_out, v_out] = ode_propagation(r0, v0, tspan, rhs_func,Param, reltol)
% -------------------------------------------------------------------------
% Title    : ODE-Based Orbit Propagation Framework (Modular RHS)
% Author   : Seonwoo Son
% Date     : 2022-08-29
% Purpose  : Generalized ODE-based orbit propagation framework
% -------------------------------------------------------------------------
% Inputs:
%   r0       - Initial position vector [3x1] (km)
%   v0       - Initial velocity vector [3x1] (km/s)
%   tspan    - Time span [t0 tf] (sec)
%   rhs_func - Function handle to RHS dynamics (default: 2-body)
%   reltol   - ODE solver relative tolerance (default: 1e-8)
%
% Outputs:
%   t_out    - Time vector (s)
%   r_out    - Propagated position [Nx3] (km)
%   v_out    - Propagated velocity [Nx3] (km/s)




    % ode_propagation 내부 (기본값 덮어쓰기 제거)
    if nargin < 4 || isempty(rhs_func), rhs_func = @rhs_2body;   end
    if nargin < 5, Param = struct;                               end  % 비어 있으면 새 struct
    if nargin < 6 || isempty(reltol) , reltol = 1e-8;            end  % 기본값 지정
    if ~isfield(Param, 'mu'),    Param.mu    = 398600.4418;      end  % [km^3/s^2]
    if ~isfield(Param, 'Re'),    Param.Re    = 6378.1363;        end  % [km]
    if ~isfield(Param, 'ApM'),   Param.ApM   = 0.01;             end  % [m^2/kg]
    if ~isfield(Param, 'Cd'),    Param.Cd    = 2.2;              end  % [unitless]
    if ~isfield(Param, 'h_ref'), Param.h_ref = 450;              end  % [km]
    if ~isfield(Param, 'H'),     Param.H     = 60.828;           end  % [km]
    if ~isfield(Param, 'we'),    Param.we    = [0; 0; 7.2921e-5];end  % [rad/s]

    y0 = [r0(:); v0(:)];
    options = odeset('RelTol', reltol);
    [t_out, y_out] = ode45(@(t, y) rhs_func(t, y,Param), tspan, y0, options);

    r_out = y_out(:,1:3);  v_out = y_out(:,4:6);
end
%%
function dydt = rhs_2body(~, y, P)
% TWO_BODY_EQUATION
% Computes the derivative of the state vector under two-body gravity.
    r = y(1:3); v=y(4:6);
    a = -P.mu / norm(r)^3 * r;
    dydt = [v; a];
end

%%
function dXdt = rhs_drag_exp(~, X, P)
% Inputs:
%   t    - Time (not used explicitly)
%   X    - State vector [6x1] = [r; v] in km and km/s
%   Aux  - Struct with required constants:
%            Aux.mu    : gravitational parameter [km^3/s^2]
%            Aux.Re    : Earth radius [km]
%            Aux.ApM   : Area-to-mass ratio [m^2/kg]
%            Aux.Cd    : Drag coefficient
%            Aux.h_ref : Reference altitude [km]
%            Aux.H     : Scale height [km]
%            Aux.we    : Earth rotation vector [rad/s]

% Defaults (if Aux is partially missing)
if ~isfield(P, 'mu');     P.mu  = 398600.4418; end
if ~isfield(P, 'Re');     P.Re  = 6378.1363;   end
if ~isfield(P, 'ApM');    P.ApM = 0.01;        end
if ~isfield(P, 'Cd');     P.Cd  = 2.2;         end
if ~isfield(P, 'h_ref');  P.h_ref = 450;       end
if ~isfield(P, 'H');      P.H   = 60.828;      end
if ~isfield(P, 'we');     P.we  = [0; 0; 7.2921e-5]; end

% State extraction
r = X(1:3); v = X(4:6); R = norm(r);


% Gravity acceleration
a_grav = -P.mu / R^3 * r;

% Air density (Exponential Model)
h   = R - P.Re;
rho = 1.585e-12 * exp(-(h - P.h_ref)/P.H); % [kg/m^3]

% Relative velocity (w.r.t rotating atmosphere)
v_rel = v - cross(P.we, r); % [km/s]

% Drag acceleration
a_drag = -0.5*P.Cd*P.ApM*(rho*1e3)*norm(v_rel)*v_rel;% [km/s^2]

% Derivative
dXdt = [v; a_grav + a_drag];

end

%%
function dXdt = rhs_drag_nrlmsise(t, X, P)
% cfg fields:
%   .datetime0  : MATLAB datetime at t=0 [UTC]
%   .mu, .Re, .Cd, .ApM
%   .we         : Earth rotation (rad/s)
%   .f107, .ap  : space-weather indices OR function handle

% --- Constants with defaults
if ~isfield(P,'mu'),  P.mu  = 398600.4418;     end
if ~isfield(P,'Re'),  P.Re  = 6378.1363;       end
if ~isfield(P,'Cd'),  P.Cd  = 2.2;             end
if ~isfield(P,'ApM'), P.ApM = 0.01;            end
if ~isfield(P,'we'),  P.we  = [0 0 7.292115e-5]'; end
if ~isfield(P,'datetime0'), P.datetime0 = datetime('now','TimeZone','UTC'); end

% --- Current UTC time
utc = P.datetime0 + seconds(t);

% --- ECI → ECEF
[r_ecef, ~] = eci2ecef(utc, X(1:3), X(4:6));

% --- LLA & density
lla = ecef2lla(r_ecef'*1e3);           % m
% [f107, ap] = getSpaceWeather(P, utc);% 사용자 정의
f107 = 150; ap = 4; % 하드 코딩 
[~, rho]  = atmosnrlmsise00(lla(3), lla(1), lla(2), year(utc), ...
              day(utc,'dayofyear'), second(utc,'secondofday'), ...
              f107, f107, ap);         % ρ(6) kg/m³
rho = rho(6);

% --- Accelerations
r  = X(1:3);  v = X(4:6);     R  = norm(r);
a_grav = -P.mu/R^3 * r;
v_rel  = v - cross(P.we, r);
a_drag = -0.5*P.Cd*P.ApM*rho*1e3*norm(v_rel)*v_rel;
dXdt = [X(4:6); a_grav + a_drag];
end


%%
function dXdt = rhs_J2(~, X, P)

    % -------- Default parameters ---------------------------------------
    if ~isfield(P,'mu'),  P.mu  = 398600.4418;     end  % km^3/s^2
    if ~isfield(P,'Re'),  P.Re  = 6378.1363;       end  % km
    if ~isfield(P,'J2'),  P.J2  = 1.0826269e-3;    end

    r_vec = X(1:3);   v_vec = X(4:6);
    r_norm = norm(r_vec);
    x = r_vec(1);  y = r_vec(2);  z = r_vec(3);

    % --- Two‑body central gravity -------------------------------------
    a_grav = -P.mu / r_norm^3 * r_vec;

    % --- J2 perturbation (gravity potential degree‑2) ------------------
    Re2 = P.Re^2;
    z2  = z^2;   r2 = r_norm^2;   r5 = r_norm^5;
    factor = -1.5 * P.J2 * P.mu * Re2 / r5;

    a_J2 = factor * [ x*(1 - 5*z2/r2);
                       y*(1 - 5*z2/r2);
                       z*(3 - 5*z2/r2) ];

    % --- Derivative -----------------------------------------------------
    dXdt = [v_vec;
            a_grav + a_J2];
end

%%
function dXdt = rhs_J2_drag_exp(~, X, P)

    % -------- Default parameters ---------------------------------------
    if ~isfield(P,'mu'),  P.mu  = 398600.4418;     end
    if ~isfield(P,'Re'),  P.Re  = 6378.1363;       end
    if ~isfield(P,'J2'),  P.J2  = 1.0826269e-3;    end
    if ~isfield(P,'Cd'),  P.Cd  = 2.2;             end
    if ~isfield(P,'ApM'), P.ApM = 0.01;            end
    if ~isfield(P,'H'),   P.H   = 60.828;          end
    if ~isfield(P,'h_ref'),P.h_ref = 450;          end
    if ~isfield(P,'we'),  P.we  = [0; 0; 7.292115e-5]; end

    % --- State extraction
    r = X(1:3);   v = X(4:6);   R = norm(r);
    x = r(1);     y = r(2);     z = r(3);

    % --- Gravity (2-body + J2)
    a_grav = -P.mu / R^3 * r;

    Re2 = P.Re^2;   z2 = z^2;   R2 = R^2;   R5 = R^5;
    factor = -1.5 * P.J2 * P.mu * Re2 / R5;
    a_J2 = factor * [ x*(1 - 5*z2/R2);
                      y*(1 - 5*z2/R2);
                      z*(3 - 5*z2/R2) ];

    % --- Atmospheric Drag (Exponential Model)
    h   = R - P.Re;
    rho = 1.585e-12 * exp(-(h - P.h_ref)/P.H);   % kg/m^3
    v_rel = v - cross(P.we, r);                  % km/s
    a_drag = -0.5 * P.Cd * P.ApM * rho * 1e3 * norm(v_rel) * v_rel;  % km/s^2

    % --- Derivative
    dXdt = [v; a_grav + a_J2 + a_drag];
end

%%
function dXdt = rhs_J2_drag_nrlmsise(t, X, P)

    % -------- Default parameters ---------------------------------------
    if ~isfield(P,'mu'),  P.mu  = 398600.4418;     end
    if ~isfield(P,'Re'),  P.Re  = 6378.1363;       end
    if ~isfield(P,'J2'),  P.J2  = 1.0826269e-3;    end
    if ~isfield(P,'Cd'),  P.Cd  = 2.2;             end
    if ~isfield(P,'ApM'), P.ApM = 0.01;            end
    if ~isfield(P,'we'),  P.we  = [0; 0; 7.292115e-5]; end
    if ~isfield(P,'datetime0'), P.datetime0 = datetime('now','TimeZone','UTC'); end

    % --- State extraction
    r = X(1:3);   v = X(4:6);   R = norm(r);
    x = r(1);     y = r(2);     z = r(3);

    % --- UTC time
    utc = P.datetime0 + seconds(t);

    % --- ECI → ECEF
    [r_ecef, ~] = eci2ecef(utc, r, v);

    % --- ECEF → LLA
    lla = ecef2lla(r_ecef'*1e3);    % input in meters

    % --- Get Space Weather
    % [f107, ap] = getSpaceWeather(P, utc);
    f107 = 150; ap = 4; % 하드 코딩 

    % --- NRLMSISE00 Density
    [~, rho_vec] = atmosnrlmsise00(lla(3), lla(1), lla(2), year(utc), ...
                     day(utc,'dayofyear'), second(utc,'secondofday'), ...
                     f107, f107, ap);
    rho = rho_vec(6);  % kg/m³

    % --- Accelerations
    a_grav = -P.mu / R^3 * r;

    Re2 = P.Re^2;   z2 = z^2;   R2 = R^2;   R5 = R^5;
    factor = -1.5 * P.J2 * P.mu * Re2 / R5;
    a_J2 = factor * [ x*(1 - 5*z2/R2);
                      y*(1 - 5*z2/R2);
                      z*(3 - 5*z2/R2) ];

    v_rel = v - cross(P.we, r);
    a_drag = -0.5 * P.Cd * P.ApM * rho * 1e3 * norm(v_rel) * v_rel;

    dXdt = [v; a_grav + a_J2 + a_drag];
end
