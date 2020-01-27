clc
clear
%close all
%path(pathdef)

% This script performs a numerical modified wavenumber analysis for a given 
% set of spatial semi-discretizations.

%% Dependencies
addpath('../../../../Extra')
addpath('../Discretization')
addpath('../Limiters')
addpath('../Physics')
addpath('../Solver')
addpath('../Grid')
addpath('../Math')

%% Test tuple setup
data = assembleTestTuple(...
    {DGIGA(1,2,1) DGIGA(2,2,1) DGIGA(3,2,1) DGIGA(2,2,0)}... semi-discretization (2)
    ,{2}... degree (3)
    ,{16}... #patches (4)
    ,{nan}... #spans (5)
    ,{nan}... continuity class within patch (6)
    ,{nan}... #DOFs (7)
    ,{nan}... wavenumbers (8)
    ,{nan}... modified wavenumbers (9)
    ,{nan}... WCT (10)
    ,{@MWA_eigen}... MWA routine to use (11)
    );

%% Batch run
I = size(data,1);
% parfor i = 1:I
for i = 1:I
    data_i = data(i,:);
    if data_i{2}.isHybrid
        data_i{5} = data_i{2}.nonzeroSpanCount;
        data_i{6} = data_i{2}.smoothness;
    else
        data_i{5} = 1;
        data_i{6} = inf;
    end
    if isnan(data_i{6}) || isinf(data_i{6}) && data_i{5} == 1
        data_i{6} = data_i{3}-1;
    end
    data_i{7} = data_i{4}*data_i{2}.basisCount; %%% Use this instead (?): data_i{8} = data_i{6}*(data_i{5}+data_i{7});
    tic
    [data_i{9},data_i{8}] = data_i{11}(data_i{4},data_i{3},data_i{2});
    data_i{10} = toc;
    % Gather data back:
    data(i,:) = data_i;
end

%% Combined plot
c = distinguishable_colors(I,{'w','k'});
figure('Renderer', 'painters', 'Position', [1000 100 700 700])
% Set up dispersion:
h(1) = subplot(2,3,[1 2]);
axis([0 pi -inf inf])
setFancyPlotSettings3
ylabel('$$\Re(\tilde{\kappa}) \frac{\Delta x}{N_b}$$','Interpreter','latex')
xlabel('$$\kappa \frac{\Delta x}{N_b}$$','Interpreter','latex')
% Set up diffusion:
h(2) = subplot(2,3,[4 5]);
axis([0 pi -inf inf])
setFancyPlotSettings3
ylabel('$$\Im(\tilde{\kappa}) \frac{\Delta x}{N_b}$$','Interpreter','latex')
xlabel('$$\kappa \frac{\Delta x}{N_b}$$','Interpreter','latex')
% Set up legend:
h(3) = subplot(2,3,[3 6]);
axis('off');
% Fill dispersion and diffusion plots:
%%% set(h,'LabelFontSizeMultiplier',1.75)
set(h,'NextPlot','add');
plot(h(1),[0 pi],[0 pi],'k--');
plot(h(2),[0 pi],[0 0],'k--');
hLeg = zeros(I,1);
for i = 1:I
    if isa(data{i,2}, 'FR')
        name = sprintf('%s(%s); p = %d',...
            class(data{i,2}),num2str(data{i,2}.param),data{i,3});
    elseif data{i,2}.isHybrid
        name = sprintf('%s; p = %d, N_\\sigma = %d, C^{%d}',...
            class(data{i,2}),data{i,3},data{i,5},data{i,6});
    else
        name = sprintf('%s; p = %d',...
            class(data{i,2}),data{i,3});
    end
    hAux = plot(h(1),data{i,8},real(data{i,9}),...
        '-','Color',c(i,:),'DisplayName',name); % dispersion
    plot(h(2),data{i,8},imag(data{i,9}),...
        '-','Color',c(i,:),'DisplayName',name); % diffusion
    hLeg(i) = hAux(1);
end
% Add the common legend:
legend(h(3),hLeg,'Location','Best')