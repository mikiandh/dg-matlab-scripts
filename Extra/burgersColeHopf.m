function u = burgersColeHopf(t,x,f,nu,x0)
% Evaluates the solution of the viscous Burgers equation at given points in
% space-time. Might be unsafe, use at your own risk. Based on Salih (2016).
% Miquel Herrera, 27/11/2020.
%
% Arguments:
%  t: sample time instants
%  x: sample locations
%  F: initial condition function
%  nu: viscosity coefficient
%  x0: array of waypoints, including domain boundaries
%
% Return:
%  u: exact solution samples (row: instant; column: position)
%
%% Safety checks
if nargin < 5
    error('Not enough inputs arguments.')
end
if nu == 0
    nu = eps;
    warning('Vanishing viscosity (will use ~%g).',nu)
end
%% Rearrange sample locations
[T,X] = meshgrid(t,x);
X0 = X(T <= 0);
X = X(T > 0);
T = T(T > 0);
%% Periodic extension
    function u = F(x)
        % Given sample locations 'x', samples 'F', its periodic extension,
        % instead of 'f' itself. This defines 'F' over the real line.
        x = x0(1) + mod(x-x0(1),x0(end)-x0(1));
        u = f(x);
    end
%% Inner functional
    function y = G(eta)
        % Evaluates the functional defined by eq. 40 in Salih (2016) at
        % a single given 'eta' parameter (vector valued).
        y = integral(@F,0,eta,'Waypoints',x0);
        y = y./(2*nu);
        y = y + (X-eta).^2./(4*nu*T);
    end

%% Solution samples
U = integral(@(eta) exp(-G(eta)).*(X-eta)./T,-inf,inf,'ArrayValued',true);
U = U./integral(@(eta) exp(-G(eta)),-inf,inf,'ArrayValued',true);

%% Rearrange sample values
u = reshape([F(X0); U],numel(x),numel(t))';
end