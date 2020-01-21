function varargout = MWA_eigen_full_GUI(varargin)
% MWA_EIGEN_FULL_GUI MATLAB code for MWA_eigen_full_GUI.fig
%      MWA_EIGEN_FULL_GUI, by itself, creates a new MWA_EIGEN_FULL_GUI or raises the existing
%      singleton*.
%
%      H = MWA_EIGEN_FULL_GUI returns the handle to a new MWA_EIGEN_FULL_GUI or the handle to
%      the existing singleton*.
%
%      MWA_EIGEN_FULL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MWA_EIGEN_FULL_GUI.M with the given input arguments.
%
%      MWA_EIGEN_FULL_GUI('Property','Value',...) creates a new MWA_EIGEN_FULL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MWA_eigen_full_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MWA_eigen_full_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MWA_eigen_full_GUI

% Last Modified by GUIDE v2.5 05-Dec-2019 14:40:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0; %%% edited this to 0 (was 1) - 2019/12/05
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MWA_eigen_full_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MWA_eigen_full_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before MWA_eigen_full_GUI is made visible.
function MWA_eigen_full_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MWA_eigen_full_GUI (see VARARGIN)

if nargin < 6
    error('Insufficient number of arguments.')
end
% Retrieve inputs:
[disc,k,kMod] = varargin{:};
% Normalize nondimensional wavenumbers to the range 0 - pi:
k = k/disc.basisCount;
kMod = kMod/disc.basisCount;
% Set up dispersion:
axes(handles.axes1);
set(handles.axes1,'NextPlot','add');
axis([0 pi -inf inf])
setFancyPlotSettings3
ylabel('$$\Re(\tilde{\kappa}) \frac{\Delta x}{N_b}$$','Interpreter','latex')
xlabel('$$\kappa \frac{\Delta x}{N_b}$$','Interpreter','latex')
% Set up diffusion:
axes(handles.axes2);
set(handles.axes2,'NextPlot','add');
axis([0 pi -inf inf])
setFancyPlotSettings3
ylabel('$$\Im(\tilde{\kappa}) \frac{\Delta x}{N_b}$$','Interpreter','latex')
xlabel('$$\kappa \frac{\Delta x}{N_b}$$','Interpreter','latex')
% Set up legend:
axes(handles.axes3);
set(handles.axes3,'NextPlot','add');
axis('off');
% Fill dispersion and diffusion plots:
%%% set(h,'LabelFontSizeMultiplier',1.75)
c = distinguishable_colors(disc.basisCount,{'w','k'});
plot(handles.axes1,[0 pi],[0 pi],'k--');
plot(handles.axes2,[0 pi],[0 0],'k--');
handles.lines1 = zeros(disc.basisCount,1); % dispersion
handles.lines2 = zeros(disc.basisCount,1); % diffusion
for i = 1:disc.basisCount
    handles.lines1(i) = plot(handles.axes1,k,real(kMod(i,:)),'-','Color',c(i,:),'DisplayName',sprintf('mode %d',i)); % dispersion
    handles.lines2(i) = plot(handles.axes2,k,imag(kMod(i,:)),'-','Color',c(i,:),'DisplayName',sprintf('mode %d',i)); % diffusion
end
% Add the common legend:
legend(handles.axes3,handles.lines1,'Location','Best')
% Add a title:
if disc.isHybrid
    contClass = disc.smoothness;
    if isinf(contClass) || isnan(contClass)
        contClass = disc.degree - 1;
    end
    title(sprintf('%s, $$N_\\sigma = %d$$, $$p = %d, C^{%d}$$',class(disc),disc.nonzeroSpanCount,disc.degree,contClass),'Interpreter','latex')
else
    title(sprintf('%s, $$p = %d$$',class(disc),disc.degree),'Interpreter','latex')
end

% Choose default command line output for MWA_eigen_full_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MWA_eigen_full_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MWA_eigen_full_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Button that toogles between single mode or all modes being shown ----
function toogleSingle_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toogleSingle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(hObject.State,'on')
    set(handles.lines1(2:end),'Visible','off')
    set(handles.lines2(2:end),'Visible','off')
else
    set(handles.lines1,'Visible','on')
    set(handles.lines2,'Visible','on')
end
