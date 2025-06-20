function [root, fx, ea, iter] = bisect(func, xl, xu, es, maxit, varargin)
%BISECT   Root finding using the bisection method.
%
%   [root, fx, ea, iter] = bisect(func, xl, xu, es, maxit, p1, p2, ...)
%
%   This function finds a root of the function 'func' within the interval
%   [xl, xu] using the bisection method. It iteratively halves the interval
%   until the estimated relative error falls below a specified tolerance.
%
%   INPUT:
%       func    - Function handle (e.g., @(x) x^2 - 4)
%       xl      - Lower guess
%       xu      - Upper guess
%       es      - Desired relative error (%) [optional, default = 0.0001]
%       maxit   - Maximum iterations allowed [optional, default = 50]
%       p1,...  - Additional parameters passed to 'func'
%
%   OUTPUT:
%       root    - Estimated root location
%       fx      - Function value at root
%       ea      - Approximate relative error (%)
%       iter    - Number of iterations performed
%
% -------------------------------------------------------------------------
%   Author  : Seonwoo Son (sonseonwoo)
%   Source  : Chapra, *Applied Numerical Methods with MATLAB*, 5th Edition
%   Date    : 2022-09-21
% -------------------------------------------------------------------------

    % --- Input validation ---
    if nargin < 3
        error('At least three input arguments (func, xl, xu) are required.');
    end

    if nargin < 4 || isempty(es), es = 1e-4; end
    if nargin < 5 || isempty(maxit), maxit = 50; end

    % Check for sign change
    fxl = func(xl, varargin{:});
    fxu = func(xu, varargin{:});
    if fxl * fxu > 0
        error('No sign change detected in the interval [%g, %g].', xl, xu);
    end

    % --- Initialization ---
    iter = 0;
    xr = xl;
    ea = 100;

    % --- Main iteration loop ---
    while true
        xrold = xr;
        xr = (xl + xu) / 2;
        iter = iter + 1;

        if xr ~= 0
            ea = abs((xr - xrold) / xr) * 100;
        end

        fr = func(xr, varargin{:});

        % Bisection update
        if fxl * fr < 0
            xu = xr;
            fxu = fr;
        elseif fxl * fr > 0
            xl = xr;
            fxl = fr;
        else
            ea = 0;
        end

        % Exit condition
        if ea <= es || iter >= maxit
            break;
        end
    end

    % --- Final output ---
    root = xr;
    fx   = func(xr, varargin{:});
end
