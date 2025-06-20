% =========================================================================
% Title    : RV-to-COE Converter with J2 Energy Term
% Author   : Seonwoo Son (sonseonwoo)
% Source   : Vallado, *Fundamentals of Astrodynamics and Applications*,
%            5th Edition
% Date     : 2023-08-23
% Purpose  :
%   Converts a 6-element inertial state vector [r; v] (km, km s⁻¹) into a
%   complete classical orbital-element set, including special-case angles
%   and specific mechanical energy with a first-order J₂ correction.
%
%   [a,e,i,RAAN,AoP,nu,lambda,longPer,u,xiJ2,xi] = rvtocoe(rv,mu,J2,Re)
%
% INPUTS
%   rv   : 6×1 or 1×6 state vector [x y z vx vy vz] (km, km s⁻¹)
%   mu   : (opt) gravitational parameter (km³ s⁻²)        (default Earth)
%   J2   : (opt) J₂ coefficient                           (default Earth)
%   Re   : (opt) equatorial radius (km)                   (default Earth)
%
% OUTPUTS (NaN for “undefined”)
%   a        semi-major axis (km)  (NaN for parabolic)
%   e        eccentricity magnitude
%   i        inclination (deg)
%   RAAN     right ascension of ascending node (deg)
%   AoP      argument of periapsis (deg)
%   nu       true anomaly (deg)
%   lambda   true longitude          (only circular-equatorial)
%   longPer  longitude of periapsis  (only elliptical-equatorial)
%   u        argument of latitude    (only circular-inclined)
%   xiJ2     specific mechanical energy incl. first-order J₂ (km² s⁻²)
%   xi       specific mechanical energy (Kepler, km² s⁻²)
% =========================================================================
function [a, e, i, RAAN, AoP, nu, lambda, longPer, u, xiJ2, xi] = ...
         RV2COE(rv, mu, J2, Re)

% ---------- default constants (WGS-84) -----------------------------------
if nargin < 2 || isempty(mu), mu  = 398600.4418;    end  % km³ s⁻²
if nargin < 3 || isempty(J2), J2  = 1.0826269e-3;   end
if nargin < 4 || isempty(Re), Re  = 6378.1363;      end % km

% ---------- input validation ---------------------------------------------
if ~isnumeric(rv) || numel(rv) ~= 6
    error('rv must be a numeric 6-element vector.');
end
rv = rv(:).';                       % force row

% ---------- primitive vectors --------------------------------------------
r = rv(1:3); v = rv(4:6);
rmag = norm(r); vmag = norm(v);
hvec = cross(r, v); hmag = norm(hvec);
khat = [0 0 1];
nvec = cross(khat, hvec); nmag = norm(nvec);
evec = ((vmag^2 - mu/rmag) * r - dot(r, v) * v) / mu;
e = norm(evec);

% ---------- energy -------------------------------------------------------
xi  = 0.5*vmag^2 - mu/rmag;
xiJ2 = xi + mu/rmag * (J2/2) * (Re/rmag)^2 * (-1 + 3*(r(3)/rmag)^2);

% ---------- semi-major axis & p -----------------------------------------
TOL_PARAB = 1e-10;
if abs(e - 1) > TOL_PARAB
    a = -mu / (2*xi);
    p = a * (1 - e^2);
else
    a = NaN;
    p = hmag^2 / mu;
end

% ---------- inclination --------------------------------------------------
i = acosd(hvec(3) / hmag);

% ---------- classification -----------------------------------------------
isCircular   = (e < 1e-12);
isEquatorial = (i < 1e-8);

% ---------- initialise outputs as NaN -----------------------------------
RAAN = NaN; AoP = NaN; nu = NaN; lambda = NaN; longPer = NaN; u = NaN;

% ---------- RAAN ---------------------------------------------------------
if ~isEquatorial
    RAAN = acosd(nvec(1)/nmag);
    if nvec(2) < 0, RAAN = 360 - RAAN; end
end

% ---------- AoP & related angles ----------------------------------------
if  isCircular &&  isEquatorial
    lambda = angleOf(r, [1 0 0]);
elseif  isCircular && ~isEquatorial
    RAAN = RAAN;
    u    = angleOf(nvec, r);
elseif ~isCircular &&  isEquatorial
    longPer = angleOf([1 0 0], evec);
    nu      = trueAnomaly(evec, r, v);
else
    AoP = angleOf(nvec, evec);
    nu  = trueAnomaly(evec, r, v);
    u      = AoP + nu;
    longPer= RAAN + AoP;
    lambda = RAAN + AoP + nu;
end

% ---------- Handle single-output case -----------------------------------
if nargout <= 1
    a = [a; e; i; RAAN; AoP; nu; lambda; longPer; u; xiJ2; xi];
end

% ---------- Nested functions ---------------------------------------------
    function ang = angleOf(vec1, vec2)
        ang = acosd(dot(vec1, vec2) / (norm(vec1)*norm(vec2)));
        if dot(cross(vec1, vec2), khat) < 0
            ang = 360 - ang;
        end
    end

    function nuDeg = trueAnomaly(evecLoc, rLoc, vLoc)
        nuDeg = acosd(dot(evecLoc, rLoc) / (norm(evecLoc)*norm(rLoc)));
        if dot(rLoc, vLoc) < 0
            nuDeg = 360 - nuDeg;
        end
    end
end
