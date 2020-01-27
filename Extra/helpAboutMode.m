function helpAboutMode(n,varargin)
% Give information on where a given mode number lies in the spectrum of the
% given grid. Sets the conventionional definition of these quantities.
%
switch nargin
    case 1
        varargin = {32 pi};
    case 2
        varargin = [varargin {1}];
    case 3
        % ok
    otherwise
        error('Too many arguments!')
end
[N,L] = varargin{:};

VarNames = {'N','L'};

disp('Grid information:');
fprintf(1, '%s\t\t\t%s\n', VarNames{:});
fprintf(1, '%g\t\t\t%g\n', [N L]');
disp(' ');

dx = L/N; % grid size
lambda = L/n; % wavelength
k = 2*pi/lambda; % wavenumber
phi = k*dx/pi*180; % phase angle (in degrees)

Data = [n lambda/dx k/dx phi];
VarNames = {'Mode (-)', 'Wavelength (dx)', 'Wavenumber (1/dx)', 'Phase angle (º)'};

disp('Mode information:');
fprintf(1, '%s\t\t%s\t\t%s\t\t%s\n', VarNames{:});
fprintf(1, '%g\t\t\t\t%g\t\t\t\t\t%g\t\t\t\t\t%g\n', Data');
disp(' ');

SCALES = {
    'Wave length (m)', Inf, 2*L, 2*L/N, L/N;
    'Wave mode (-)', 0, 1/2, N/2, N; % note that aliasing kicks in at N/2 (waves up to mode N are detected, but they are confused with modes < N/2)
    'Wave number (1/m)', 0, pi/L, pi*N/L, 2*pi*N/L;
    'Phase angle (º)', 0, 180/N, 180, 360;
    };
T = cell2table(SCALES,'VariableNames',{'Scale','Constant','Largest','Aliasing','Smallest'});
disp('Grid resolution limits:');
disp(T);
disp(' ');
end