classdef BDF < Limiter
    % Hierarchical moment limiter by Biswas, Devine and Flaherty (1994).
    % Free of user-defined parameters. Does not guarantee any TVD or TVB 
    % property (not even in the means).
    %
    properties (SetAccess = immutable)
        % Biswas et al., 1994 suggest carrying out 2 limiter sweeps each
        % time the limiter is called. This ensures that the modes used to
        % determine if limiting is required have themselves been limited.
        sweepCount
    end
    methods
        %% Constructor
        function this = BDF(varargin)
            % Superclass constructor:
            this = this@Limiter(varargin{:});
            % Initialize an input parser:
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'sweeps',2);
            % Parse the inputs:
            parse(p,varargin{:});
            this.sweepCount = p.Results.sweeps;
        end
        %% Apply (extension)
        function apply(this,mesh,~)
            % Default limiting:
            apply@Limiter(this,mesh);
            % Retrieve troubled elements:
            elements = findobj(mesh.elements,'isTroubled',true)';
            % Safety check:
            if isempty(elements)
                return 
            end
            % Limiter sweeps:
            for sweep = 1:this.sweepCount
                % Get/set limited Legendre coefficients of each element:
                coefs{2} = this.getLimitedLegendre(elements(1));
                for k = 2:length(elements)
                    coefs{1} = this.getLimitedLegendre(elements(k));
                    this.setLimitedLegendre(elements(k-1),coefs{2});
                    coefs{2} = coefs{1};
                end
                this.setLimitedLegendre(elements(end),coefs{2});
            end
        end
    end
    methods (Access = protected)
        %% Extract and limit coefficients
        function coefs = getLimitedLegendre(this,element)
            % Applies hierarchical limiting a la Biswas et al., 1994, on
            % the Legendre coefficients extracted from a given element. All
            % Legendre coefficients are limited in a single vectorized 
            % operation and then returned. The element's solution is not
            % modified, but the limited components are marked as such.
            %
            % Get unlimited legendre coefficents of current element:
            coefs = element.basis.getLegendre(element);
            % Aliases:
            j = 1:element.dofCount-1; % lower-order coefficient indices
            elementL = element.edgeL.elementL;
            elementR = element.edgeR.elementR;
            weights = 1./(2*j-1); % coefficient difference weights
            % Get left-wise differences:
            if elementL.dofCount >= j(end) % no padding
                diffsL = coefs(:,j) - elementL.basis.getLegendre(elementL,j);
            else % zero-padded
                aux = 1:elementL.dofCount;
                diffsL = coefs(:,j);
            	diffsL(:,aux) = diffsL(:,aux) - elementL.basis.getLegendre(elementL);
            end
            % Get righ-wise differences:
            if elementR.dofCount >= j(end) % no padding
                diffsR = elementR.basis.getLegendre(elementR,j) - coefs(:,j);
            else % zero-padded
                aux = 1:elementR.dofCount;
                diffsR = -coefs(:,j);
            	diffsR(:,aux) = elementR.basis.getLegendre(elementR) + diffsR(:,aux);
            end
            % Get mapping operators to/from local characteristic variables:
            [~,L,R] = this.physics.getEigensystemAt(coefs(:,1));
            % Apply the minmod operator coefficient-wise:
            aux = R*this.minmod(L*coefs(:,j+1),weights.*(L*diffsL),weights.*(L*diffsR));
            % Limit down to the highest untroubled mode:
            i = true(size(aux,1),1);
            for r = size(aux,2):-1:1
                i = abs(aux(i,r) - coefs(i,r+1)) > 1e-10;
                if all(i == 0)
                    break
                end
                element.isLimited(:,r) = i;
                coefs(i,r+1) = aux(i,r);
            end
        end
    end
    methods (Static, Access = protected)
        %% Set limited coefficients
        function setLimitedLegendre(element,coefs)
            % Find all rows of the given element's state matrix that have
            % been limited and replace them by the projection to the 
            % element's basis of the given matrix of Legendre coefficients.
            %
            rows = any(element.isLimited,2);
            element.basis.setLegendre(element,coefs(rows,:),rows);
        end
    end
    methods (Static)
        %% Information
        function info = getInfo
            info = 'BDF';
        end
        %% Minmod operator (Osher, 1984)
        function A = minmod(A,B,C)
            % Minmod function for 3 matrix inputs. Operates entry-wise.
            % "Inlined" for speed. All three matrices are assumed to have 
            % consistent sizes. 
            %
            A = (sign(A) == sign(B) & sign(A) == sign(C)).*sign(A).*min(abs(A),min(abs(B),abs(C)));
        end
    end
end