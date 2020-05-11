clc
clear
%close all

% This script plots a 1D (nonrational) B-spline curve using a B-spline 
% basis and a set of control points.

%% Parameters
N = 50; % arc length samples
p = 3; % degree
knots = [0 1 2 2 3 3 3 4]; % knot vector
cPoints = logspiral(1,.3063489,10,1); % control points

%% Sample
% Preprocess:
reps = [p zeros(1,numel(knots)-2) p]; % knot multiplicities
knots = repelem(knots,reps+1); % clamped knot vector
breaks = unique(knots); % breakpoint vector
% Check sizes:
if size(cPoints,1) ~= length(knots)-p-1
    error(['The number of control points must match the number of basis functions (' num2str(length(knots)-p-1) ').'])
end
% Plot the control polygon:
plot(cPoints(:,1),cPoints(:,2),'--ob','MarkerFaceColor','b','MarkerSize',10); % control points
hold all
% Plot the B-spline curve:
C = zeros(N*(numel(breaks)-1),2); % B-spline sample coordinates
B = zeros(numel(breaks),2); % breakpoint coordinates
for i = 1:numel(breaks)-1
    xi = linspace(breaks(i),breaks(i+1),N);
    for j = 1:numel(xi)
         C(j+(i-1)*N,:) = curvePoint(knots,cPoints,p,xi(j));
    end
    B(i,:) = C(1+(i-1)*N,:);
end
B(end,:) = C(end,:);
plot(B(:,1),B(:,2),'sr','MarkerFaceColor','r','MarkerSize',10) % breakpoints
plot(C(:,1),C(:,2),'-k') % B-spline curve
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