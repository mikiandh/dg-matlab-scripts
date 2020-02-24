classdef TVDM < Limiter
    %
    % Slope limiter proposed by Cockburn and Shu in 1989. Proven to be TVD
    % in the element-wise means. Known to cause loss of accuracy at smooth
    % extrema. Only compatible with a Legendre basis, single element per 
    % patch, p > 0.
    % 
    % For p > 1, the slope of the fictituous line joining the two edge
    % values of the cell is used in the detection step (a la Zhu et al., 
    % 2013).
    %
    methods
        %% Constructor
        function this = TVDM(physics,varargin)
            this.physics = physics;
        end
        %% Apply limiter
        function apply(this,mesh,varargin)
            % Compatibility check:
            troubled = this.isCompatible(mesh.elements);
            % Apply element-wise:
            I = this.physics.equationCount;
            for element = mesh.elements(troubled)
                % Interpolate solution at element edges:
                element.interpolateStateAtEdges;
                % Aliases:
                N = element.dofCount;
                uL = element.edgeL.elementL.states(:,1); % left cell means
                uR = element.edgeR.elementR.states(:,1); % right cell means
                u = element.states(:,1); % current cell means
                vL = element.stateL; % unlimited left edge intercept
                vR = element.stateR; % unlimited right edge intercept
                % Limited edge intercepts:
                wL = u - this.minmod(u - vL,u-uL,uR-u); % left edge
                wR = u + this.minmod(vR - u,u-uL,uR-u); % right egde
                % Determine troubled components:
                eqs = find(wL ~= vL & wR ~= vR);
                % Overwrite troubled slopes with limited slopes:
                element.isLimited = false(I,N);
                for i = eqs'
                    s = .5*(vR(i) - vL(i)); % unlimited slope
                    sL = u(i) - uL(i); % left-limited slope
                    sR = uR(i) - u(i); % right-limited slope
                    element.states(i,2) = this.minmod(s,sL,sR); % most conservative slope
                    element.states(i,3:end) = 0; % reduce to linear basis
                    element.isLimited(i,2:N) = true;
                end
            end
        end
    end
    methods (Static)
        %% Compatibilty check
        function tagged = isCompatible(elements)
            % Untags any element not using a modal basis.
            tagged = Limiter.isCompatible(elements);
            tagged(tagged) = ...
                cellfun(@(x) x.isModal, {elements(tagged).basis});
        end
        %% Minmod (3 column array arguments)
        function d = minmod(a,b,c)
            d = zeros(length(a),1);
            M = abs(horzcat(a,b,c));
            a = sign(a);
            ids = a == sign(b) & a == sign(c);
            d(ids) = min(M(ids,:),[],2);
            d = a.*d;
        end
    end
end