function varargout = diagnostic_signal(varargin)
% DIAGNOSTIC_SIGNAL MATLAB code for diagnostic_signal.fig
%      DIAGNOSTIC_SIGNAL, by itself, creates a new DIAGNOSTIC_SIGNAL or raises the existing
%      singleton*.
%
%      H = DIAGNOSTIC_SIGNAL returns the handle to a new DIAGNOSTIC_SIGNAL or the handle to
%      the existing singleton*.
%
%      DIAGNOSTIC_SIGNAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIAGNOSTIC_SIGNAL.M with the given input arguments.
%
%      DIAGNOSTIC_SIGNAL('Property','Value',...) creates a new DIAGNOSTIC_SIGNAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before diagnostic_signal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to diagnostic_signal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help diagnostic_signal

% Last Modified by GUIDE v2.5 22-May-2018 15:47:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @diagnostic_signal_OpeningFcn, ...
                   'gui_OutputFcn',  @diagnostic_signal_OutputFcn, ...
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


% --- Executes just before diagnostic_signal is made visible.
function diagnostic_signal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to diagnostic_signal (see VARARGIN)

% Choose default command line output for diagnostic_signal
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

addpath(genpath('.'));
try
    img = imread('src/ecg2.png');
    image(img);
catch
    textLabel = sprintf('Error');
    set(handles.text_main, 'String', textLabel);
end

% UIWAIT makes diagnostic_signal wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = diagnostic_signal_OutputFcn(hObject, eventdata, handles) 
% % varargout  cell array for returning output args (see VARARGOUT);
% % hObject    handle to figure
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in push.
function push_Callback(hObject, eventdata, handles)

[file,path] = uigetfile('*.mat', 'rt');
signal = load(fullfile(path, file));
handles.data = signal.ecg;
handles.Fs = signal.Fs;
handles.N = size(handles.data,2);
handles.time_axis = (1:handles.N)/handles.Fs;
plot(handles.time_axis, handles.data); grid on;
xlabel('Time (s)');
ylabel('Magnitude');
textLabel = sprintf('The electrocardiogram seems good.');
set(handles.display_error, 'String', textLabel);

try
[bpm,perc_of_p_value_AF,perc_sample_brady,perc_sample_tachy,percent_of_extopic_beat,gamma] = pathologies_detection( handles.data, handles.Fs)
    % Brady % 
  
    brady_display = sprintf('%f',perc_sample_brady);
    set(handles.brady_result, 'String', brady_display);
     if perc_sample_brady > 50
        text1 = sprintf('Warning : Bradycardia Risks. ');
        set(handles.display_error,'String',text1);
    end

    % Tachy % 
    tachy_display = sprintf('%f', perc_sample_tachy);
    set(handles.tachy_result, 'String', tachy_display);
    if perc_sample_tachy > 50
        text2 = sprintf('Warning : Tachycardia Risks. ');
        set(handles.display_error,'String',text2);
    end

    % Ectopic % 
    ectopic_display = sprintf('%f', percent_of_extopic_beat);
    set(handles.ectopic_result, 'String', ectopic_display);
     if percent_of_extopic_beat >= 25 
        text2 = sprintf(' Warning : significative ectopics beat ');
        set(handles.display_error,'String',text2);
     end

    % BPM 

    bpm_display = sprintf('%f ', bpm);
    set(handles.bpm_result, 'String', bpm_display);

    % P_value % 
    p_value_display = sprintf('%f ',perc_of_p_value_AF );
    set(handles.p_peaks_result, 'String', p_value_display);
    if perc_of_p_value_AF < 30 
        if perc_sample_tachy >50
            text = sprintf(' Warning : Atrial Fibrillation & Tachycardia Risks. Check AutoCov  ');
            set(handles.display_error,'String',text);
        else
            text = sprintf(' Warning : Arterial Fibrillation Risks.  ');
            set(handles.display_error,'String',text);
            
            
        end
    end
catch
     textLabel = sprintf('Error.');
     set(handles.display_error, 'String', textLabel);
end



guidata(hObject, handles);
function figure1_CreateFcn(hObject, eventdata,handles)


% --- Executes on button press in push_auto.
function push_auto_Callback(hObject, eventdata, handles)
% hObject    handle to push_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%  handles.gamma  = str2func(gamma); 
%  handles.time_axis = (1:handles.N)/handles.Fs;
[bpm,perc_of_p_value_AF,perc_sample_brady,perc_sample_tachy,percent_of_extopic_beat,gamma] = pathologies_detection( handles.data, handles.Fs)
try
    figure(1);
    plot(gamma);
    xlabel('data');
    ylabel('autocovariance');
    title('If a dirac happen, there is an Atrial Fibrillation.');
    
    
catch
     textLabel = sprintf('Error.');
     set(handles.display_error, 'String', textLabel);
end
guidata(hObject, handles);
