%% SSP RK time update
function fun = ESSPRK(s,p,LIMITER)
%ESSPRK SSP RK explicit time-step update.
% returns a function handle that will update a semi-discrete ODE according
% to a selected RK multi-stage scheme of the SSP kind, if available. If a
% limiter function is provided, it will be used after every stage of the
% schemes.
%
if nargin == 2
    LIMITER = @(u) u; % no limiting (default)
end
switch p
    case 1
        fun = @(dt,u,f) EulerForward(dt,u,f,LIMITER);
    case 2
        fun = @(dt,u,f) SSPRK2(dt,u,f,LIMITER);
    case 3
        fun = @(dt,u,f) SSPRK33(dt,u,f,LIMITER);
    case 4
        fun = @(dt,u,f) SSPRK104(dt,u,f,LIMITER);
    otherwise
        error('Unknown time scheme!')
end
end

%% SSP RK Methods
function u = EulerForward(dt,u,f,LIMITER)
u = u + dt*f(u);
u = LIMITER(u);
end

function u = SSPRK2(dt,u,f,LIMITER)
u0 = u;
u = u + dt*f(u);
u = LIMITER(u);
u = 0.5*(u0 + u + dt*f(u));
u = LIMITER(u);
end

function u = SSPRK33(dt,u,f,LIMITER)
u0 = u;
u = u + dt*f(u);
u = LIMITER(u);
u = 0.75*u0 + 0.25*(u + dt*f(u));
u = LIMITER(u);
u = 1/3*u0 + 2/3*(u + dt*f(u));
u = LIMITER(u);
end

function u = SSPRK104(dt,u,f,LIMITER)
u0 = u; % u0 = u^{n}
u = u0+1/6*dt*f(u); % u1
u = LIMITER(u);
for i = 2:4
    u = u+1/6*dt*f(u); % u2, u3, u4
    u = LIMITER(u);
end
u4 = u; % non-optimal memory requirement (additional stage stored)
u = 0.6*u0+0.4*u+1/15*dt*f(u); % u5
u = LIMITER(u);
for i = 6:9
    u = u+1/6*dt*f(u); % u6, u7, u8, u9
    u = LIMITER(u);
end
u = 0.04*u0 + 0.36*u4 + 0.6*u + 0.06*dt*f(u4) + 0.1*dt*f(u); % u10 = u^{n+1}
u = LIMITER(u);
end