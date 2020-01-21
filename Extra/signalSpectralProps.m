clc, clear, close all

%% Description
% Script that analyses the spectral properties of a given signal, on a
% given mesh (1D).
%
% See: Hirsch 2007, pp. 295
%

%% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid:
L = 1; % physical domain size
N = 128; % EVEN number of degrees of freedom (physical: nodes; spectral: modes)
x = linspace(0,L,N);

helpAboutMode(N/8,N,L);

% Signal matrix (rows => signals; columns => nodal values):
f = [
    %wavePacket(x,L/2,10,1,0);
    wavePacket(x,L/2,N/8,exp(-32*(x-0.5*L).^2),0);
    %exp(-32*(x-0.5*L).^2);
    wavePacket(x,L/2,[1 3 8 21 55],1/15*exp(-32*(x-0.5*L).^2).*[5 4 3 2 1]',0)
    ];

% Signal names:
names = {
    %'\phi = \pi/6.25 packet';
    'Tonal pulse';
    %'Gaussian';
    'Fibonacci pulse'
    };

% Signal colors:
colors = {'b','r','g','c','m','y'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Grid resolution bounds
% SCALES = {
%     'Wave length', Inf, 2*L, 2*L/N;
%     'Wave mode', 0, 1/2, N; % note that aliasing kicks in at N/2 (waves up to mode N are detected, but they are confused with modes < N/2 because of aliasing)
%     'Wave number', 0, pi/L, pi*N/L;
%     'Phase angle', 0, pi/N, pi;
%     };
% T = cell2table(SCALES,'VariableNames',{'Scale','Constant','Largest','Smallest'});
% disp('Grid resolution limits:');
% disp(T);
% disp(' ');

%% Physical domain
plotSignals(x,f,names,colors);

%% Spectral domain
n = 0:1:N/2;
F = fft(f,N,2); % Fourier transform of the signals
P2 = abs(F/N); % double sided spectrum
P1 = P2(:,1:N/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1); % single sided spectrum
plotSignalSpectra(n,P1,names,colors)


%% Functions
function plotSignals(x,signals,names,colors)
figure('Renderer', 'painters', 'Position', [100 100 600 300])
hold on
for i = 1:size(signals,1)
    plot(x,signals(i,:),colors{i},'DisplayName',names{i});
end
hold off
setFancyPlotSettings3
xlabel('Grid location, x');
ylabel('Signal, f(x)');
title('Signals in physical space')
legend('-DynamicLegend','Location','Best');
end

function plotSignalSpectra(n,spectra,names,colors)
figure('Renderer', 'painters', 'Position', [100 100 600 300])
b = bar(n,spectra');
for i = 1:size(spectra,1)
    b(i).FaceColor = colors{i};
    b(i).DisplayName = names{i};
    %plot(n,spectra(i,:),colors{i},'DisplayName',names{i});
end
setFancyPlotSettings1
xlabel('Mode, n');
ylabel('Amplitude, F_n');
title('Signals in spectral space')
legend('-DynamicLegend','Location','Best');
end