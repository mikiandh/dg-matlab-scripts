clc
clear
close all

% This script plots a 1D (nonrational) B-spline curve using a B-spline 
% basis and a set of control points.

%% Parameters
N = 250; % arc length samples
p = 3; % degree
knots = [0 1 2 2 3 3 3 4]; % knot vector
cPoints = logspiral(2,.1,10,1.5); % control points

%% B-spline curve
% Open the knot vector:
reps = [p zeros(1,numel(knots)-2) p];
knots = repelem(knots,reps+1);
% Check sizes:
if size(cPoints,1) ~= length(knots)-p-1
    error(['The number of control points must match the number of basis functions (' num2str(length(knots)-p-1) ').'])
end
% Sample locations (along the arc length):
x = [knots linspace(knots(1),knots(end),N)];
% Sample the curve:
tic
C = zeros(length(x),size(cPoints,2));
for i = 1:length(x)
    C(i,:) = curvePoint(knots,cPoints,p,x(i));
end
toc
disp(['Done (' num2str(toc) ').'])

%% Plot
figure()
plot(cPoints(:,1),cPoints(:,2),'--ob','MarkerFaceColor','b','MarkerSize',10) % control points
hold on
plot(C(length(knots)+1:end,1),C(length(knots)+1:end,2),'-k') % B-spline curve
plot(C(1:length(knots),1),C(1:length(knots),2),'sr','MarkerFaceColor','r','MarkerSize',10) % break points
%legend('Control polygon','B-spline curve','Location','Best')
axis equal

%% Functions
% Piegl & Tiller, 1997; A2.1:
function n = findSpan(X,p,x)
% Returns the knot span to which the given sample point
% belongs. Not vectorized. Only for open knot vectors.
%
% Arguments:
%  p: polynomial degree
%  X: knot vector
%  x: sample location
% Output:
%  n: knot span that contains the sample point x
%
% Special case:
n = length(X)-p-1; % id of the knot at the left of the right breakpoint
if (x == X(end))
    return;
end
% Binary search:
low = p+1;
high = n+1;
n = floor(0.5*(low + high));
while x < X(n) || x >= X(n+1)
    if x < X(n)
        high = n;
    else
        low = n;
    end
    n = floor(0.5*(low + high));
end
end

% Piegl & Tiller, 1997; A2.2:
function B = basisFuns(X,p,n,x)
% Sample all basis functions of given degree that are non-zero at a given 
% location (in a given knot span of a given knot vector).
%
% Arguments
%  X: knot vector
%  p: basis functions degree
%  n: knot span index, such that: X(n) <= x < X(n+1)
%  x: sample location
% Output
%  B: row array of p+1 basis functions of degree p, sampled at x; all
%     others (of degree p) are zero.
%
% Preallocations:
left = zeros(1,p);
right = zeros(1,p);
B = ones(1,p+1); % trivial case (p = 0)
% Loop over p > 1 degrees:
for j=1:p
    left(j) = x-X(n+1-j);
    right(j) = X(n+j)-x;
    saved = 0;
    for r=1:j
        temp = B(r)/(right(r)+left(j-r+1));
        B(r) = saved+right(r)*temp;
        saved = left(j-r+1)*temp;
    end
    B(j+1) = saved;
end
end

% Piegl & Tiller, 1997; A3.1:
function C = curvePoint(X,P,p,x)
% Evaluate a B-spline curve at a location by sampling its B-spline basis
% according to a set of control points.
%
% Arguments
%  X: knot vector
%  P: 2D array of control points (row: control point; column: space dimension)
%  p: basis functions degree
%  x: sample location in "knot vector space"
% Output
%  C: curve sample (row array; column: space dimension)
%
% Find the knot span index such that X(n) <= x < X(n+1):
n = findSpan(X,p,x);
% Sample the basis at x:
B = basisFuns(X,p,n,x);
% Sample the B-spline curve:
C = B*P(n-p:n,:); % if x is in the n-th knot span, only the n-th basis function and the p previous ones will be non-zero ("triangle of influence").
end

% Logarithmic spiral:
function z = logspiral(a,k,n,m)
% Distributes n points along a logatithmic spiral of parameters 'a' and 'k'
t = linspace(0,m*2*pi,n)';
z = a*exp(k*t).*[cos(t) sin(t)];
end