classdef TVB < Limiter
    %
    % Slope limiter as reported in Cockburn and Shu, 2001. Proven to be TVB
    % in the element-wise means (TVD if M = 0). Known to cause loss of 
    % accuracy at smooth extrema unless M is set appropriately.
    %
    properties
        M % user-definable sensitivity parameter (see Cockburn & Shu, 2001)
    end
    methods
        %% Constructor
        function this = TVB(varargin)
            % Superclass constructor:
            this = this@Limiter(varargin{:});
            % Initialize an input parser:
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'M',0);
            % Parse the M parameter:
            parse(p,varargin{:});
            this.M = p.Results.M;
        end
        %% Apply
        function apply(this,mesh,~)
            % Apply default limiting:
            apply@Limiter(this,mesh);
            % Retrieve troubled elements:
            elements = findobj(mesh.elements,'isTroubled',true)';
            % Apply on every remaining element:
            for element = elements
                this.applyOnElement(element);
            end
        end
    end
    methods (Access = protected)
        function applyOnElement(this,element)
            % Applies TVB limiting (based on the TVD limiter of Osher) a la
            % Cockburn & Shu, 2001.
            %            
            % Cell averages:
            uL = element.basis.getLegendre(element.edgeL.elementL,1); % left cell
            uR = element.basis.getLegendre(element.edgeR.elementR,1); % right cell
            u = element.basis.getLegendre(element,1); % current cell
            % Unlimited edge intercepts:
            element.interpolateStateAtEdges;
            vL = element.stateL; % left edge
            vR = element.stateR; % right edge
            % Limited edge intercepts:
            wL = u - this.minmod(u-vL,u-uL,uR-u); % left edge
            wR = u + this.minmod(vR-u,u-uL,uR-u); % right edge
            % Determine troubled components:
            rows = find(wL ~= vL | wR ~= vR);
            cols = 3:element.dofCount; % higher order modes (above linear)
            % Overwrite troubled slopes with limited ones:
            for i = rows'
                s = .5*(vR(i) - vL(i)); % unlimited slope
                sL = u(i) - uL(i); % left-limited slope
                sR = uR(i) - u(i); % right-limited slope
                s = this.minmod(s,sL,sR); % most conservative limited slope
                element.basis.setLegendre(element,s,2,i); % override slope
                element.basis.setLegendre(element,0,cols,i); % truncate to linear
            end
        end
    end
    methods (Static)
        %% Minmod
        function [d,flag] = minmod(a,b,c)
            % Minmod function for 3 column array arguments. Returns the
            % column array of results and a logical 2D array that is
            % true at the position of each returned value.
            %
            d = zeros(length(a),1);
            M = abs([a b c]);
            a = sign(a);
            flag = a == sign(b) & a == sign(c);
            d(flag) = min(M(flag,:),[],2);
            flag = d == M;
            d = a.*d;
        end
    end
end