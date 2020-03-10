classdef BSB < BDF
    % Hierarchical moment limiter by Burbeau, Sagaut and Bruneau (2001).
    % Slight improvement over Biswas et al., 1994 (less diffusive near 
    % discontinuities, according to the authors).
    %
    methods
        %% Constructor
        function this = BSB(varargin)
            % Superclass constructor:
            this = this@BDF(varargin{:});
        end
    end
    methods (Access = protected)
        %% Apply equation-wise (override)
        function applyNoSync(this)
            % Applies the limiter asynchronously, i.e. on each 
            % characteristic variable separately.
            %
            % Loop over characteristic variables:
            for i = 1:this.I
                k = true(this.K,1); % flag every element for limiting
                % Limit hierarchically (top to bottom):
                for j = this.J:-1:1
                    % Get "minmodded" coefs.:
                    coefs_min = this.minmod((2*j-1).*this.coefs(k,j+1,i),...
                        this.coefs(k,j,i) - this.coefsL(k,j,i),...
                        this.coefsR(k,j,i) - this.coefs(k,j,i))./(2*j-1);
                    % Update list of elements flagged for limiting:
                    k0 = abs(this.coefs(k,j+1,i) - coefs_min) > 1e-10 | ~this.coefs(k,j+1,i);
                    k(k) = k0;
                    % Get "max" coefs.:
                    coefs_max = this.minmod((2*j-1).*this.coefs(k,j+1,i),...
                        this.coefs(k,j,i) - this.coefsL(k,j,i) - (2*j-1).*this.coefsL(k,j+1,i),...
                        this.coefsR(k,j,i) - (2*j-1).*this.coefsR(k,j+1,i) - this.coefs(k,j,i))./(2*j-1);
                    % Get "maxmodded" coefs.:
                    this.coefs(k,j+1,i) = this.maxmod(coefs_min(k0),coefs_max);
                end
            end
        end
    end
    methods (Static)
        %% Information
        function info = getInfo
            info = 'BSB';
        end
        %% Maxmod operator (Roe, 1985)
        function A = maxmod(A,B)
            % Maxmod function for 2 matrix inputs. Operates entry-wise.
            % "Inlined" for speed. The two matrices must have consistent 
            % sizes.
            %
            A = (sign(A) == sign(B)).*sign(A).*max(abs(A),abs(B));
        end
    end
end