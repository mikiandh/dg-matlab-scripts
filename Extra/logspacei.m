function x = logspacei(a,b,N)
% Returns a 1D array of n quasi-logarithmically spaced integers between
% values a and b.
narginchk(3,3)
if b-a < N-1
    error('Interval [%d,%d] contains fewer than %d integers.',a,b,N)
end
x = ones(1,N); % preallocation
r = (b/a)^(1/N); % initial growth ratio
for n = 1:N-1 % linear portion
    if r*x(n) - x(n) >= 1
        break % skip to log portion
    else
        x(n+1) = x(n) + 1; % linear portion
        r = (b/x(n+1))^(1/(N-n));
    end
end
x(n:end) = logspace(log10(x(n)),log10(b),N-n+1); % log portion
x = round(x); % guaranteed to be unique
end
