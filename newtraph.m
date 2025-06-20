function [root, ea, iter] = newtraph(func, dfunc, xr, es, maxit, varargin)
%NEWTRAPH   Root finding using the Newton-Raphson method.
%
%   [root, ea, iter] = newtraph(func, dfunc, xr, es, maxit, p1, p2, ...)
%
%   This function uses the Newton-Raphson method to find a root of the
%   function 'func', given its derivative 'dfunc' and an initial guess.
%
%   INPUT:
%       func    - Function handle (e.g., @(x) x^3 - 2*x - 5)
%       dfunc   - Derivative of the function (e.g., @(x) 3*x^2 - 2)
%       xr      - Initial guess
%       es      - Desired relative error (%) [optional, default = 0.0001]
%       maxit   - Maximum iterations [optional, default = 50]
%       p1,...  - Additional parameters passed to 'func' and 'dfunc'
%
%   OUTPUT:
%       root    - Estimated root
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
        error('At least three input arguments (func, dfunc, xr) are required.');
    end
    if nargin < 4 || isempty(es),    es    = 1e-4; end
    if nargin < 5 || isempty(maxit), maxit = 50;  end

    % ----- Initialization -----
    iter = 0;
    ea   = 100;

    % ----- Iteration loop -----
    while true
        xrold = xr;

        fval  = func(xr, varargin{:});
        dfval = dfunc(xr, varargin{:});

        if dfval == 0
            error('Derivative is zero at x = %.6f. Method fails.', xr);
        end

        xr = xr - fval / dfval;
        iter = iter + 1;

        if xr ~= 0
            ea = abs((xr - xrold) / xr) * 100;
        end

        if ea <= es || iter >= maxit
            break;
        end
    end

    % ----- Output -----
    root = xr;
end
