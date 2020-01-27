classdef (Abstract) Physics < handle
    properties
        equationCount
        boundaryConditionsFunction
    end
    methods (Abstract)
        flux(this,states)
        riemannFlux(this,stateL,stateR)
        applyBoundaryConditions(this,mesh)
        getConvectionMatrix(this,states,basis)
    end
    methods
        function displayData(this,rowNames,colNames,varargin)
            % Function that prints a table of data (with arbitrary number
            % of columns and rows) in which each position contains a column
            % array with as many positions as equations in the current PDE.
            %
            % Arguments:
            %  rowNames: 1D cell array with character arrays to use as 
            %            row headers.
            %  colNames: idem, column headers.
            %  varargin: matching number of cell arrays, each containing a 
            %            column array of values to display.
            %
            rows = length(rowNames);
            cols = length(colNames)-1;
            % Check arguments:
            if length(varargin) ~= cols
                error('Incompatible arguments.')
            end
            % Re-arrange variable argument list into a [rows,cols] array:
            data = reshape(cell2mat(varargin),this.equationCount*rows,cols);
            
            % Generate variable-length format specifiers:
            specs = ['%-14s\t' repmat('%-14.4g\t',1,cols) '\n'];
            fprintf(1,'%-14s\t',colNames{:})
            fprintf(1,'\n')
            for r = 1:rows
                for i = 1:this.equationCount
                    if i == 1
                        fprintf(1,'\n')
                        fprintf(1,specs,rowNames{r},data(sub2ind([this.equationCount rows],i,r),:));
                    else
                        fprintf(1,specs,' ',data(sub2ind([this.equationCount rows],i,r),:));
                    end
                end
            end
            disp(' ')
        end
    end
end