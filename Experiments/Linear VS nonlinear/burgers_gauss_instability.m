clc
clear
%close all

% This script is generates data I use for the 3D time-slice plots comparing
% DGSEM, DGIGA and DGIGA-AFC in regards to their tendency to undergo that
% "pesky instability" I don't quite understand.

%% Setup runs
testStruct = cell2struct({
    20  DGSEM(2)                    Limiter     [0,.1,.2,.3,.43]        6
    20  DGIGA(1,2)                  Limiter     [0,.1,.2,.3,.43]        6
    20  DGIGA_nodal(1,2)            Limiter     [0,.1,.2,.3,.43]        6
    20  DGIGA_AFC(1,2)              AFC_2010    [0,.1,.2,.3,.43]        6
    20  DGIGA_nodal_AFC(1,2)        AFC_2010    [0,.1,.2,.3,.43]        6
    3   DGSEM(19)                   Limiter     [0,.1,.16,.22,.33,.43]  3
    3   DGIGA(2,14,9)               Limiter     [0,.1,.16]              3
    3   DGIGA_nodal(2,14,9)         Limiter     [0,.1,.16,.22]          3
    3   DGIGA_AFC(2,14,9)           AFC_2010    [0,.1,.16,.22,.33,.43]  3
    3   DGIGA_nodal_AFC(2,14,9)     AFC_2010    [0,.1,.16,.22,.33,.43]  3
    },{'K','basis','limiter','t','n'},2);

%% Batch run
parfor i = 1:numel(testStruct)
    
    %% Discretization
    mesh = Mesh(testStruct(i).basis,[-1 1],Periodic(2),testStruct(i).K);
    
    %% Solver
    solver = SSP_RK3(Burgers,[0 0],...
        'timeDelta',1e-3,...
        'norm',Norm({'L2','TV','L1','ErrorL2'}),...
        ...'exactSolution',@(t,x) Burgers.MOC(t,x,@(x) (t < 0.2).*(1 + .1*sin(16*pi*x)),[-1,1]),...
        'exactSolution',@(t,x) Burgers.MOC(t,x,@(x) (t < 0.5).*exp(-9/4*pi*x.^2),[-1,1]),...
        ...'exactSolution',@(t,x) sawtoothBurgersExact(t,x,0),...
        'limiter',testStruct(i).limiter);
    
    %% Initial condition
    solver.initialize(mesh)
    
    %% Time-integration
    for t = testStruct(i).t
        solver.timeStop = t;
        dt = solver.timeDelta;
        solver.launch(mesh)
        solver.timeDelta = dt;
        filename = sprintf('burgers_gauss_J=%d_%s',mesh.bases(1).basisCount,class(mesh.bases(1)));
        solver.writeSolutionToFile(filename,testStruct(i).n);
    end
end