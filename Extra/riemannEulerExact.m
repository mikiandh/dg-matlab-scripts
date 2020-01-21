function [r,u,p] = riemannEulerExact(t,x,rL,uL,pL,rR,uR,pR,x0,g)
% Solves the exact Riemann problem for the 1D Euler equations and a perfect
% gas. From Toro 2009, chapter 4. Miquel Herrera, 10-06-2019.
%
% Arguments:
%  t: sample time instants
%  x: sample locations
%  x0: position of the intial interface between L and R states
%  r,u,p: primitive variables: density, velocity and pressure (L: left, R:
%         right)
%
% Return:
%  r: row array of densities (column: position)
%  u: idem, velocity
%  p: idem, pressure

%% Some default arguments
if nargin < 9
    % ok
else
    x = x - x0; % map interface position
end
if nargin < 10
    g = 1.4;
end

%% Some preallocations
t = reshape(t,[],1);
x = reshape(x,1,[]);
I = length(t);
J = length(x);
r = zeros(I,J);
u = zeros(I,J);
p = zeros(I,J);

%% Trivial solution (t = 0)
r(t == 0,x <= 0) = rL;
r(t == 0,x > 0) = rR;
u(t == 0,x <= 0) = uL;
u(t == 0,x > 0) = uR;
p(t == 0,x <= 0) = pL;
p(t == 0,x > 0) = pR;

%% Some data-dependent constants
% Gamma-dependent constants (precomputed for speed):
g8 = g-1;
g1 = g8/(2*g);
g2 = (g+1)/(2*g);
g3 = 2*g/g8;
g4 = 2/g8;
g5 = 2/(g+1);
g6 = g8/(g+1);
g7 = 0.5*g8;
% Eq. 4.8:
AL = g5/rL;
BL = g6*pL;
AR = g5/rR;
BR = g6*pR;
% Speeds of sound:
aL = sqrt(g*pL/rL);
aR = sqrt(g*pR/rR);

%% Vacuum creation check
if g4*(aL + aR) <= uR - uL % eq. 4.40
    error('Initial conditions violate the pressure positivity condition.')
end

%% Newton-Raphson (pressure in the star region)
fun = @(p) funL(p) + funR(p) + uR - uL; % eq. 4.5
dFun = @(p) dFunL(p) + dFunR(p); % eq. 4.37
p0 = pressureInitialGuess; % initial guess
while 1
    pStar = p0 - fun(p0)/dFun(p0); % eq. 4.44
    if 2*abs(pStar - p0)/(pStar + p0) < 1e-6
        break % converged on tolerance
    end
    p0 = pStar;
end

%% Velocity in the star region (i.e. speed of the contact discontinuity)
uStar = 0.5*(uR + uL + funR(pStar) - funL(pStar)); % eq. 4.9

%% Sampling
s = x./t; % similarity variable
ijL = s < uStar & ~isnan(s) & ~isinf(s); % samples left of contact, t > 0
ijR = s >= uStar & ~isnan(s) & ~isinf(s); % samples right of contact, t > 0
% Left wave:
if pStar > pL % left wave is a shock
    pp = pStar/pL;
    rStarL = rL*(pp + g6)/(g6*pp + 1); % eq. 4.50
    sL = uL - aL*sqrt(g2*pp + g1); % eq. 4.52
    ij = ijL & s > sL; % between contact and left shock
    r(ij) = rStarL;
    u(ij) = uStar;
    p(ij) = pStar;
    ij = ijL & s < sL; % left of the left shock
    r(ij) = rL;
    u(ij) = uL;
    p(ij) = pL;
else % left wave is an expansion fan
    pp = pStar/pL;
    rStarL = rL*pp^(1/g); % eq. 4.53
    aStarL = aL*pp^g1; % eq. 4.54
    sHL = uL - aL; % eq. 4.55a
    sTL = uStar - aStarL; % eq. 4.55b
    ij = ijL & s > sTL; % between contact and tail of the fan
    r(ij) = rStarL;
    u(ij) = uStar;
    p(ij) = pStar;
    ij = ijL & s < sTL & s > sHL; % between tail and head of the fan
    r(ij) = rL*(g5 + g6/aL*(uL - s(ij))).^g4; % eq. 4.56a
    u(ij) = g5*(aL + g7*uL + s(ij)); % eq. 4.56b
    p(ij) = pL*(g5 + g6/aL*(uL - s(ij))).^g3; % eq. 4.56c
    ij = ijL & s < sHL; % left of the head of the left fan
    r(ij) = rL;
    u(ij) = uL;
    p(ij) = pL;
end
% Right wave:
if pStar > pR % right wave is a shock
    pp = pStar/pR;
    rStarR = rR*(pp + g6)/(g6*pp + 1); % eq. 4.57
    sR = uR + aR*sqrt(g2*pp + g1); % eq. 4.59
    ij = ijR & s < sR; % between contact and right shock
    r(ij) = rStarR;
    u(ij) = uStar;
    p(ij) = pStar;
    ij = ijR & s > sR; % right of the right shock
    r(ij) = rR;
    u(ij) = uR;
    p(ij) = pR;
else
    pp = pStar/pR;
    rStarR = rR*pp^(1/g); % eq. 4.60
    aStarR = aR*pp^g1; % eq. 4.61
    sHR = uR + aR; % eq. 4.62a
    sTR = uStar + aStarR; % eq. 4.62b
    ij = ijR & s < sTR; % between contact and tail of the fan
    r(ij) = rStarR;
    u(ij) = uStar;
    p(ij) = pStar;
    ij = ijR & s > sTR & s < sHR; % between tail and head of the fan
    r(ij) = rR*(g5 - g6/aR*(uR - s(ij))).^g4; % eq. 4.63a
    u(ij) = g5*(-aR + g7*uR + s(ij)); % eq. 4.63b
    p(ij) = pR*(g5 - g6/aR*(uR - s(ij))).^g3; % eq. 4.63c
    ij = ijR & s > sHR; % right of the head of the right fan
    r(ij) = rR;
    u(ij) = uR;
    p(ij) = pR;
end

%% Some auxiliary functions
% Functions connecting left, right and star regions:
    function f = funL(p) % eq. 4.6
        if p > pL % shock
            f = (p-pL)*sqrt(AL/(p+BL));
        else % expansion
            f = g4*aL*(power(p/pL,g1)-1);
        end
    end
    function f = funR(p) % eq. 4.7
        if p > pR % shock
            f = (p-pR)*sqrt(AR/(p+BR));
        else % expansion
            f = g4*aR*(power(p/pR,g1)-1);
        end
    end
% Derivatives of the above:
    function f = dFunL(p) % eq. 4.37a
        if p > pL % shock
            f = sqrt(AL/(p+BL))*(1-0.5*(p-pL)/(BL+p));
        else % expansion
            f = power(p/pL,-g2)/(rL*aL);
        end
    end
    function f = dFunR(p) % eq. 4.37b
        if p > pR % shock
            f = sqrt(AR/(p+BR))*(1-0.5*(p-pR)/(BR+p));
        else % expansion
            f = power(p/pR,-g2)/(rR*aR);
        end
    end
% Provide a good first guess for pressure in the star region:
    function p0 = pressureInitialGuess
        p0 = 0.5*(pL+pR)-0.125*(uR-uL)*(rL+rR)*(aL+aR);
        p0 = max(0,p0);
        pMin = min(pL,pR);
        pMax = max(pL,pR);
        if pMax/pMin <= 2.0 && pMin <= p0 && p0 <= pMax
            return % use primitive variables estimate (eq. 4.47)
        end
        if p0 < pMin % use two-rarefaction estimate (eq. 4.46)
            p0 = aL/pL^g1+aR/pR^g1;
            p0 = (aL + aR - g7*(uR-uL))/p0;
            p0 = p0^g3;
        else % use two-shock estimate (eq. 4.48)
            p0 = (pL*sqrt(AL/(p0 + BL)) + pR*sqrt(AR/(p0 + BR)) + ...
                uL - uR)/(sqrt(AL/(p0 + BL)) + sqrt(AR/(p0 + BR)));
            p0 = max(0,p0);
        end
    end
end