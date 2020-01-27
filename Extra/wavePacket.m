function f = wavePacket(x,L,n,A,f0)
% Returns a wave packet made up of a linear combination of pure tones, 
% concentrated within the span L, centered in the domain x.
%
%  inputs: grid, domain length, mode array, amplitude matrix, offset level
%  output: signal
%
x0 = 0.5*(x(1) + x(end) - L);
L0 = x(1) - x(end);
n = n(:);
if size(A,2) ~= 1 && size(A,2) ~= length(x)
    A = A'; % place modes along columns
end
f = sum(A.*sin(n.*2*pi/L0*x),1);
f = f0 + f.*(heaviside(x-x0) - heaviside(x-x0-L));
end