% =========================================================================
% Title    : COE to RV Converter (PQW → ECI)
% Author   : Seonwoo Son (sonseonwoo)
% Source   : Vallado, *Fundamentals of Astrodynamics and Applications*, 5th Edition
% Date     : 2023-08-23
% Purpose  :
%   Given a set of classical orbital elements (COE), return the position
%   and velocity vectors in the Earth‑centred inertial (IJK/ECI) frame.
% =========================================================================
%
% [r_ijk, v_ijk] = coe2rv(a, e, i_deg, RAAN_deg, AoP_deg, nu_deg, mu)
%
% INPUTS
%   a         semi‑major axis [km]  (for parabolic orbits pass p instead
%             and set a = NaN; see note below)
%   e         eccentricity magnitude (0 <= e < 1   ellipse
%                                      e  = 1      parabola
%                                      e  > 1      hyperbola)
%   i_deg     inclination [deg]
%   RAAN_deg  right ascension of ascending node [deg]
%   AoP_deg   argument of periapsis [deg]
%   nu_deg    true anomaly [deg]
%   mu        gravitational parameter [km^3 s^-2] (default: 398600.4418)
%
% OUTPUTS
%   r_ijk     3×1 position vector in ECI [km]
%   v_ijk     3×1 velocity vector in ECI [km/s]
%
% NOTE
%   • For parabolic orbits (e ≈ 1) where the semi‑major axis is undefined,
%     supply the parameter p = h^2/μ in place of 'a' **and** set a = NaN.
% -------------------------------------------------------------------------
function [r_ijk, v_ijk] = coe2rv(a, e, i_deg, RAAN_deg, AoP_deg, nu_deg, mu)

    % ---------- defaults & tolerance ------------------------------------
    if nargin < 7 || isempty(mu)
        mu = 398600.4418;           % WGS‑84 Earth μ [km^3/s^2]
    end

    % ---------- parameter checks ----------------------------------------
    if ~isnumeric([a,e,i_deg,RAAN_deg,AoP_deg,nu_deg]) || numel([a,e])<2
        error('All COE parameters must be numeric scalars.');
    end

    % ---------- semilatus rectum p --------------------------------------
    if isnan(a) || abs(e-1) < 1e-12
        % Parabolic or p supplied directly as ''a''
        p = a;  % here user must pass p directly
    else
        p = a * (1 - e^2);
    end

    % ---------- angle conversions to rad --------------------------------
    i     = deg2rad(i_deg);
    RAAN  = deg2rad(RAAN_deg);
    AoP   = deg2rad(AoP_deg);
    nu    = deg2rad(nu_deg);

    % ---------- position & velocity in PQW ------------------------------
    r_pqw = [ p * cos(nu) / (1 + e*cos(nu));
              p * sin(nu) / (1 + e*cos(nu));
              0 ];

    v_pqw = [ -sqrt(mu/p) * sin(nu);
               sqrt(mu/p) * (e + cos(nu));
               0 ];

    % ---------- rotation matrices ---------------------------------------
    R3 = @(theta) [  cos(theta)  sin(theta) 0;
                    -sin(theta)  cos(theta) 0;
                          0          0      1 ];
    R1 = @(theta) [ 1      0           0;
                    0  cos(theta)  sin(theta);
                    0 -sin(theta)  cos(theta) ];

    Q_pX = R3(-RAAN) * R1(-i) * R3(-AoP);   % PQW → ECI

    % ---------- output vectors ------------------------------------------
    r_ijk = Q_pX * r_pqw;
    v_ijk = Q_pX * v_pqw;
end
