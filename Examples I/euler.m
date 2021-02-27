clc
clear
%close all
%path(pathdef)

% This script solves the Euler equations.

%% Conveniency stuff
sod_BCs = Farfield([1; 0; 2.5],[0.125; 0; 0.25]);
toro1_BCs = Farfield([1; 0.75; 2.7813],[1; 0; 0.25]);
toro2_BCs = Farfield([1; -2; 3],[1; 2; 3]);
toro3_BCs = Farfield([1; 0; 2500],[1; 0; 0.25]);
shuOsher_BCs = Farfield([3.8571; 10.1419; 39.1667],[0.9735; 0; 2.5]);

%% Discretization
% mesh = Mesh(DGSEM(2),[-5 5],shuOsher_BCs,200);
mesh = Mesh(DGIGA_AFC(1,3),[-5 5],shuOsher_BCs,300);
% mesh.bases.diffusionFun = @DGIGA_AFC.diffusionRobust;

%% Solver
solver = SSP_RK4_10(Euler('HLLC'),[0 1.8],...
    'exactSolution',@shuOsher,...
    ...'Limiter',Krivodonova('Sensor',KXRCF,'Stats',true),...
    ...'Limiters', [Krivodonova('Sensor',KXRCF,'Stats',true) EulerP1 EulerP0],...
    ...'Limiter',[AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true)  EulerP1_step('Stats',true)   EulerP0_step('Stats',true)],...
    'Limiter', AFC_2010('Sensor',KXRCF,'Control',[1 3 2],'Failsafe',3,'Stats',true),...
    'iterSkip',25,...
    'equations',1);
solver.courantNumber = .5*solver.optimizeCFL(mesh.bases);
% solver.timeDelta = 1e-3;
% solver.isTimeDeltaFixed = true;

%% Time-integration
% solver.initialize(mesh,'Method','interpolate')
solver.initialize(mesh)
solver.launch(mesh)

%% Exact solutions and/or initial conditions
% Simple density jump:
function y = jump(t,x)
[r,u,p] = riemannEulerExact(t,x,2,0,1,1,0,1,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
% Sod's shock tube:
function y = sod(t,x)
[r,u,p] = riemannEulerExact(t,x,1,0,1,0.125,0,0.1,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
% Leveque 2002 (designed to avoid non-positive density due to undershoots):
function y = leveque(t,x)
[r,u,p] = riemannEulerExact(t,x,3,0.9,3,1,0.9,1,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
% Toro (2009)'s collection (table 11.1):
function y = toro1(t,x)
% tEnd = 0.2
[r,u,p] = riemannEulerExact(t,x,1,0.75,1,0.125,0,0.1,0.3);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
function y = toro2(t,x)
% tEnd = 0.15
[r,u,p] = riemannEulerExact(t,x,1,-2,0.4,1,2,0.4,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
function y = toro2bis(t,x)
% Adapted to avoid non-positive density due to undershoots.
% tEnd = 0.15
[r,u,p] = riemannEulerExact(t,x,1,-.5,0.4,1,.5,0.4,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
function y = toro3(t,x)
% tEnd = 0.012
[r,u,p] = riemannEulerExact(t,x,1,0,1000,1,0,0.01,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
function y = toro4(t,x)
% tEnd = 0.035
[r,u,p] = riemannEulerExact(t,x,5.99924,19.5975,460.894,5.99242,-6.19633,46.095,0.4);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
function y = toro5(t,x)
% tEnd = 0.012
[r,u,p] = riemannEulerExact(t,x,1,-19.59745,1000,1,-19.59745,0.01,0.8);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
% Woodward and Colella:
function y = woodwardColella(~,x)
% tEnd = 0.038
% L = [0 1]
% Reflective(0,0)
y = ones(3,length(x));
y1 = Euler.primitiveToState([1 0 1000]');
y2 = Euler.primitiveToState([1 0 0.01]');
y3 = Euler.primitiveToState([1 0 100]');
y(:,x < 0.1) = y1.*y(:,x < 0.1);
y(:,x >= 0.1 & x < 0.9) = y2.*y(:,x >= 0.1 & x < 0.9);
y(:,x >= 0.9) = y3.*y(:,x >= 0.9);
end
% Shu and Osher (mimicks shock-turbulence interaction):
function y = shuOsher(~,x)
% tEnd = 1.8
% L = [-5 5]
y = ones(3,length(x));
y(:,x < -4) = [3.857143; 2.629369; 10.33333].*y(:,x < -4);
y(1,x >= -4) = 1+0.2*sin(5*x(x >= -4));
y(2,x >= -4) = 0.*y(2,x >= -4);
y = Euler.primitiveToState(y);
end
% Lax:
function y = lax(t,x)
% tEnd = 1.3
% L = [-5 5]
[r,u,p] = riemannEulerExact(t,x,0.445,0.698,3.528,0.5,0,0.571,0);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
% Smooth sine wave:
function y = sineSpeed(~,x)
y = ones(3,length(x));
y(2,:) = y(2,:).*sin(2*pi*x);
end
% Riemann "1-2-3" problem
function y = riemann123(t,x)
[r,u,p] = riemannEulerExact(t,x,1,2,.4,1,-2,.4,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
% Riemann "blast waves" (left) problem
function y = riemannBlastL(t,x)
[r,u,p] = riemannEulerExact(t,x,1,-0,1000,1,0,1000,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end
% Riemann "blast waves" (right) problem
function y = riemannBlastR(t,x)
[r,u,p] = riemannEulerExact(t,x,1,0,100,1,-0,100,0.5);
y = [r; r.*u; p/0.4 + 0.5*r.*u.^2];
end