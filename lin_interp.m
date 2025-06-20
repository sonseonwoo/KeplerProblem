function [root, fx, ea, iter] = lin_interp(func, xl, xu, es, maxit, varargin)
%LIN_INTERP   Root finding with the false-position (linear interpolation) method.
%
%   [root, fx, ea, iter] = lin_interp(func, xl, xu, es, maxit, p1, p2, ...)
%
%   Uses linear interpolation between the bounds xl and xu (Regula-Falsi)
%   to locate a root of the scalar function 'func'.
%
%   INPUT:
%       func    - Function handle (e.g., @(x) x.^3 - x - 2)
%       xl      - Lower guess
%       xu      - Upper guess
%       es      - Desired relative error (%)    [optional, default = 0.0001]
%       maxit   - Maximum iterations allowed    [optional, default = 50]
%       p1,...  - Extra parameters forwarded to 'func'
%
%   OUTPUT:
%       root    - Estimated root
%       fx      - Function value at the root
%       ea      - Approximate relative error (%)
%       iter    - Number of iterations performed
%
% -------------------------------------------------------------------------
%   Author  : Seonwoo Son (sonseonwoo)
%   Source  : Chapra, *Applied Numerical Methods with MATLAB*, 5th Edition
%   Date    : 2022-09-21
% -------------------------------------------------------------------------

    % ----- Input validation -----
    if nargin < 3
        error('At least three input arguments (func, xl, xu) are required.');
    end
    if nargin < 4 || isempty(es),    es    = 1e-4; end
    if nargin < 5 || isempty(maxit), maxit = 50;  end

    % Evaluate function at initial bounds
    fxl = func(xl, varargin{:});
    fxu = func(xu, varargin{:});
    if fxl * fxu > 0
        error('No sign change detected in the interval [%g, %g].', xl, xu);
    end

    % ----- Initialization -----
    iter = 0;
    xr   = xl;
    ea   = 100;

    % ----- Main iteration loop -----
    while true
        xrold = xr;
        % False-position formula (linear interpolation)
        xr = xu - fxu * (xl - xu) / (fxl - fxu);
        iter = iter + 1;

        if xr ~= 0
            ea = abs((xr - xrold) / xr) * 100;
        end

        fr = func(xr, varargin{:});

        % Bracketing update
        if fxl * fr < 0
            xu  = xr;
            fxu = fr;
        elseif fxl * fr > 0
            xl  = xr;
            fxl = fr;
        else
            ea = 0;  % Exact root found
        end

        % Termination condition
        if ea <= es || iter >= maxit
            break;
        end
    end

    % ----- Outputs -----
    root = xr;
    fx   = func(xr, varargin{:});
end
