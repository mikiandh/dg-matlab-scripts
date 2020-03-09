classdef BSB < BDF
    % Hierarchical moment limiter by Burbeau, Sagaut and Bruneau (2001).
    % Gradual improvement over Biswas et al., 1994 (less diffusive near 
    % discontinuities, according to the authors).
    %
    methods
        %% Constructor
        function this = BSB(varargin)
            % Superclass constructor:
            this = this@BDF(varargin{:});
        end
    end
    
end