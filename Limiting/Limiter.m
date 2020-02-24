classdef Limiter < handle
    %
    % Base class of all DG limiters, i.e. those which can be applied to a
    % mesh containing one or more DG patches in which p > 0. Any valid 
    % limiter must inherit this class.
    %
    properties (Access = protected)
        % Each limiter must be aware of the physics being solved; a 
        % compatibility check between the two is made at construction time.
        physics
    end
    methods (Abstract)
        % Method called by the time marching routine once per stage to 
        % limit the solutions in the entire mesh.
        apply(this,mesh,timeDelta,stage,stageCount)
    end
    methods (Static)
        %% Factory constructor
        function this = Limiter(varargin)
            % Given a number of input name-value pairs, decides which
            % limiter to instantiate and returns it.
            %
            %%% TO DO %%%
        end
        %% Default compatibility check
        function tagged = isCompatible(elements)
            % Given an array of element handles, returns an array of 
            % logicals that is true at all elements in which the current
            % limiter can be applied. 
            %
            % Discards all p = 0 elements (default implementation).
            %
            % Find p > 0 elements:
            tagged = [elements.dofCount] > 1;
            % Un-tag left/right-most (until BCs are implemented properly):
            tagged([1 end]) = false;
        end
    end
end