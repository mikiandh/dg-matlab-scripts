classdef TVDM < Limiter
    %
    % Slope limiter proposed by Cockburn and Shu in 1989. Proven to be TVD
    % in the element-wise means. Known to cause loss of accuracy at smooth
    % extrema. Only compatible with a Legendre basis, single element per 
    % patch, p > 0.
    % 
    % For p > 1, the slope of the fictituous line joining the two edge
    % values of the cell is used in the TVD detection step (a la Zhu et al. 
    % 2013).
    %
    methods
        %% Limiting
        function apply(this,mesh,~)
            % Apply base class limiter:
            apply@Limiter(this,mesh);
            % Discard left/right-most elements (lacking a better option):
            set(mesh.elements([1 end]),'isSensed',false);
            % Retrieve troubled elements:
            elements = findobj(mesh.elements,'isSensed',true);
            % Discard any element with a non-modal basis:
            elements(arrayfun(@(x) ~x.basis.isModal,elements)) = [];
            % Apply on every remaining element:
            for element = elements'
                this.applyOnElement(element);
            end
        end
    end
    methods (Access = protected)
        function applyOnElement(this,element)
            % Interpolate solution at element edges:
            element.interpolateStateAtEdges;
            % Aliases:
            uL = element.edgeL.elementL.states(:,1); % left cell means
            uR = element.edgeR.elementR.states(:,1); % right cell means
            u = element.states(:,1); % current cell means
            vL = element.stateL; % unlimited left edge intercept
            vR = element.stateR; % unlimited right edge intercept
            % Limited edge intercepts:
            wL = u - this.minmod(u-vL,u-uL,uR-u); % left edge
            wR = u + this.minmod(vR-u,u-uL,uR-u); % right edge
            % Determine troubled components:
            eqs = find(wL ~= vL | wR ~= vR);
            % Overwrite troubled slopes with limited slopes:
            for i = eqs'
                s = .5*(vR(i) - vL(i)); % unlimited slope
                sL = u(i) - uL(i); % left-limited slope
                sR = uR(i) - u(i); % right-limited slope
                element.states(i,2) = this.minmod(s,sL,sR); % most conservative slope
                element.states(i,3:end) = 0; % reduce to linear basis
                element.isLimited(i,2:element.dofCount) = true;
            end
        end
    end
    methods (Static, Access = protected)
        %% Minmod (3 scalar arguments)
        function d = minmod(a,b,c)
            d = sign(a);
            if d ~= sign(b) || d ~= sign(c)
                d = 0;
            else
                d = d*min(abs([a b c]));
            end
        end
    end
end