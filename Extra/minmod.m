function out = minmod(varargin)
% Minmod function. Fully vectorized. Accepts any number of inputs, each of
% them having an arbitrary number of dimensions (equal to that of the 
% others!). Very generic, but slower than a more specific implementation.
%
dims = size(varargin{1}); % array of size per dimension of the arguments
d = length(dims)+1; % concatenation dimension
% Arrange arguments in as a 2D array:
M = cat(d,varargin{:});
n = prod(dims); % total number of entries in each argument
M = reshape(M,n,nargin);
% Identify columns of the final M matrix with all equal signs:
S = sign(M); % 2D array of signs
plusCols = ismember(S,ones(1,nargin),'rows'); % columns where all entries are positive
minusCols = ismember(S,-ones(1,nargin),'rows'); % columns where all entries are positive
% Final form of M matrix:
M = M';
out = zeros(1,n);
% Find minimum absoute value along columns, keeping the sign:
out(plusCols) = min(M(:,plusCols));
out(minusCols) = max(M(:,minusCols));
% Recover original shape:
dimsCellArray = num2cell(dims);
out = reshape(out,dimsCellArray{:});
end