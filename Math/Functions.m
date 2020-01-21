classdef Functions < handle
    methods (Static)
        %% Gaussian pulse
        function fval = gauss(x,xL,xR,x0)
            % Gaussian pulse; centered at x0, amplitude (in some sense) is
            % 1/3 of [xL,xR]. See A5-VI, 4/10/2019.
            % Defaults:
            if nargin < 2, xL = 0; end
            if nargin < 3, xR = xL + 1; end
            if nargin < 4, x0 = .5*(xL + xR); end
            % Sample:
            fval = exp(-18*(.5*sqrt(2*pi))^2*(x-x0).^2/(xR-xL)^2);
        end
        %% Discontinuous pulse
        function fval = jump(x,xL,xR)
            % Rectangular pulse; central 1/3 portion of [xL,xR].
            % Defaults:
            if nargin < 2, xL = 0; end
            if nargin < 3, xR = xL + 1; end
            % Sample:
            fval = heaviside(x-2/3*xL-1/3*xR) - heaviside(x-1/3*xL-2/3*xR);
        end
        %% Polychromatic plane wave
        function fval = tones(x,n,A,xL,xR)
            % Linear superposition of monochromatic planar waves i.e.
            % tones.
            %
            % Defaults:
            if nargin < 2, n = 1; end
            if nargin < 3, A = ones(size(n))./length(n); end
            if nargin < 4, xL = 0; end
            if nargin < 5, xR = xL + pi; end
            % Re-arrange input:
            L = xR - xL;
            x = reshape(x,1,[]);
            n = reshape(n,[],1);
            A = reshape(A,[],1);
            if size(A,2) ~= 1 && size(A,2) ~= length(x)
                A = A'; % place modes along columns
            end
            fval = sum(A.*sin(n.*2*pi/L*x),1);
        end
    end
end