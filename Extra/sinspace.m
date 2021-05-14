function y = sinspace(d1, d2, n, factor)
%SINSPACE sine spaced vector.
%   SINSPACE(X1, X2) generates a row vector of 100 sine spaced points
%   between X1 and X2.
%
%   SINSPACE(X1, X2, N) generates N points between X1 and X2.
%
%   A sine spaced vector clusters the elements toward X2:
%    X1 |      |      |      |     |    |   |  | || X2
%
%   Make n negative to cluster the elements toward X1:
%    X1 || |  |   |    |     |      |      |      | X2
%
%   For -2 < N < 2, SINSPACE returns X2.
%
%   SINSPACE(X1, X2, N, W) clusters the elements to a lesser degree as
%   dictated by W. W = 0 returns a normal sine spaced vector. W = 1 is the
%   same as LINSPACE(X1, X2, N). Experiment with W < 0 and W > 1 for
%   different clustering patterns.
%
%   Author:     Sky Sartorius
%
%   See also COSSPACE, LINSPACE, LOGSPACE.

% Copyright (c) 2010, Sky Sartorius
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if nargin == 2
    n = 100;
end
if nargin < 4
    factor = false;
end
if n<0
    n = floor(double(-n));
    y = d2 + (d1-d2)*sin(pi/2*(1-(0:(n-1))/(n-1)));
else
    n = floor(double(n));
    y = d1 + (d2-d1)*sin(pi/(2*(n-1))*(0:n-1));
end
if factor
    y = (1-factor)*y+factor*[d1+(0:n-2)*(d2-d1)/(n-1) d2];
end
end