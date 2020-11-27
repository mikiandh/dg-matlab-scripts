clc
clear
close all

% This script computes norms of a numerical solution of the Advection
% equation as a function of time, using a monochromatic initial condition.

%% Setup
fileNameRoot = 'advection_sine';
Ndofs = 60; % e.g. 300, 180 or 60
dt = 1e-3; % e.g. 1e-3
T = 0:1:20; % instants when to sample the norms
bases = [
    % Baseline  % Optimal FR                       % Optimal DGIGA 	% Nodal DGIGA
    DGSEM(2)	FR({'eta',0.0227795117943271},2)   DGIGA(1,2,1)     DGIGA_nodal(1,2,1)
    DGSEM(3)	FR({'eta',0.0351005869100722},3)   DGIGA(1,3,2)     DGIGA_nodal(1,3,2)
    DGSEM(4)	FR({'eta',0.047091694344944},4)    DGIGA(1,4,3)     DGIGA_nodal(1,4,3)
    DGSEM(5)	FR({'eta',0.0580283806157087},5)   DGIGA(2,3,1)     DGIGA_nodal(2,3,1)
    DGSEM(7)	FR({'eta',0.0786884021016807},7)   DGIGA(2,5,3)     DGIGA_nodal(2,5,3)
    DGSEM(10)	FR({'eta',0.102887003589808},10)   DGIGA(2,7,4)     DGIGA_nodal(2,7,4)
    DGSEM(14)	FR({'eta',0.128444040756296},14)   DGIGA(2,10,6)    DGIGA_nodal(2,10,6)
    DGSEM(19)	FR({'eta',0.152203674556189},19)   DGIGA(2,14,9)    DGIGA_nodal(2,14,9)
    ];
export = struct('fig',true,'dat',true,'tbl',true); % yes or no?
targetWellResolvedFraction = .9; % try to nail the point of maximum effect (see combined-mode plots)

%% Batch run
tbl(1:size(bases,1)) = {table(T(:))};
tic
parfor i = 1:size(bases,1)
    % Some additional setup:
    bases_cpu = bases(i,:); % assign to worker
    norms = Norm({'TV','L1','L2','ErrorL2'});
    basisNames = arrayfun(@class,bases(i,:),'UniformOutput',false);
    colNames = compose('%s_%s',norms.string',string(basisNames));
    tbl{i} = array2table([T(:) zeros(numel(T),numel(colNames))]); %#ok<PFBNS>
    tbl{i}.Properties.VariableNames = ['t'; colNames(:)];
    for basis = bases_cpu
        % Mesh:
        mesh = Mesh(basis,[-1 1],Periodic(2),round(Ndofs/basis.basisCount));
        % Initial condition:
        nyquistWavemode = mesh.dofCount/2;
        resolvingWavemode = bases_cpu(1).getResolvingWavenumber*mesh.elementCount/(2*pi); % resolving eff. of the BASELINE
        targetWavemode = round(targetWellResolvedFraction*resolvingWavemode);
        fprintf(1,'Starting %d x %s\nTargeted wavemode: %g\n %g of well-resolved\n %g of resolved (Nyquist)\n\n',...
            mesh.elementCount,basis.getName,targetWavemode,targetWavemode/resolvingWavemode,targetWavemode/nyquistWavemode)
        exactSolution = @(t,x) advectionExact(t,x,@(x) 1+.1*sin(pi*targetWavemode*x),1,[-1,1]);
        % Solver:
        solver = SSP_RK3(Advection,[0 0],'timeDelta',dt,'norms',norms,'exactSolution',exactSolution);
        % Time-integration:
        solver.initialize(mesh)
        ticID = tic;
        for j = 1:numel(T)
            solver.timeStop = T(j); % set next stop
            solver.timeDelta = dt; % reset
            solver.launch(mesh) % re-launch
            tbl{i}{j,endsWith(tbl{i}.Properties.VariableNames,class(basis))} = [norms.vals]; % store results
        end
        fprintf(1,'Finished %d x %s\n Elapsed time: %g seconds\n\n',mesh.elementCount,basis.getName,toc(ticID))
        % Export figure:
        if export.fig %#ok<PFBNS>
            savefig(sprintf('%s_J=%d_%s',fileNameRoot,basis.basisCount,class(basis)))
        end
        % Export solution:
        if export.dat
            % Generate file:
            fileID = fopen(...
                sprintf('%s_t=%g_J=%d_%s.dat',fileNameRoot,solver.timeNow,basis.basisCount,class(basis)),...
                'wt');
            % Write info and header:
            fprintf(fileID,'# %s\n',solver.getInfo);
            fprintf(fileID,'# %d x %s; N_dofs = %d\n',mesh.elementCount,basis.getName,mesh.dofCount);
            fprintf(fileID,'# Targeted wavemode: %g\n',targetWavemode);
            fprintf(fileID,'# Resolving wavemode: %g\n',basis.getResolvingWavenumber*mesh.elementCount/(2*pi));
            fprintf(fileID,'# Nyquist wavemode: %g\n',mesh.dofCount/2);
            aux = [norms.cellstr; {norms.vals}];
            fprintf(fileID,'# %s: %g\n',aux{:});
            fprintf(fileID,'%s \t %s \t %s \t %s \n','x','approx','exact','exact_initial');
            % Write solution to file, element-wise:
            for element = mesh.elements
                x = linspace(element.xL,element.xR,250);
                fprintf(fileID,'%g \t %g \t %g \t %g \n',[x
                    element.interpolateStateAtCoords(x)
                    solver.exactSolution(solver.timeNow,x)
                    solver.exactSolution(0,x)]);
                fprintf(fileID,'\n');
            end
            fclose(fileID);
        end
    end
    % Export data:
    if export.dat
        writetable(tbl{i},sprintf('%s_J=%d.dat',fileNameRoot,basis.basisCount),'Delimiter','\t')
    end
end
toc

%% Export tables:
if export.tbl
    save(fileNameRoot,'tbl')
end

%% Plot tables:
for i = 1:size(bases,1)
    figure
    plot(tbl{i}.t,tbl{i}.ErrorL2_DGSEM,'^-','DisplayName',bases(i,1).getName)
    hold on
    plot(tbl{i}.t,tbl{i}.ErrorL2_FR,'s-','DisplayName',bases(i,2).getName)
    plot(tbl{i}.t,tbl{i}.ErrorL2_DGIGA,'d-','DisplayName',bases(i,3).getName)
    plot(tbl{i}.t,tbl{i}.ErrorL2_DGIGA_nodal,'o-','DisplayName',bases(i,4).getName)
    hold off
    legend('Location','Best')
    title(sprintf('Advection, N_{dofs} ~ %d, J = %d',Ndofs,bases(i,1).basisCount))
end