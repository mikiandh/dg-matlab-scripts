function data = assembleTestTuple(varargin)
% Assembles all combinations of an arbitrary number of cell arrays of
% any number of variables each into a 2D cell array, such that each row
% is a combination and each column is an input cell array.
%
% Based on the "fundamental principle of counting":
% https://stackoverflow.com/questions/3536833/arbitrary-number-of-nested-loops
%
% Preallocation:
N = cellfun('length',varargin); % number of data points in each category
M = prod(N); % total number of combinations taking one data point from each category
data = cell(M,1+nargin);
% Loop over combinations:
n = ones(size(N)); % set each category's counter
for i = 1:M
    % Store each combination's counter:
    data{i,1} = i;
    % Loop over categories:
    for j = 1:nargin
        data{i,1+j} = varargin{j}{n(j)}; % select one cell from each category
    end
    % Advance each category's counter:
    for k = 1:nargin
        n(k) = n(k) + 1;
        if n(k) > N(k)
            n(k) = 1; % if it surpasses its maximum, reset it and advance next category
        else
            break % if not, do not advance the remaining categories
        end
    end
end