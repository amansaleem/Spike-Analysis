function varargout = TraceInspector(varargin)
% TRACEINSPECTOR inspects traces
% 
% TraceInspector( prot )
%
%      TRACEINSPECTOR, by itself, creates a new TRACEINSPECTOR or raises the existing
%      singleton*.
%
%      H = TRACEINSPECTOR returns the handle to a new TRACEINSPECTOR or the
%      handle to
%      the existing singleton*.
%
%      TRACEINSPECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACEINSPECTOR.M with the given input arguments.
%
%      TRACEINSPECTOR('Property','Value',...) creates a new TRACEINSPECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TraceInspector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TraceInspector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TraceInspector

% Last Modified by GUIDE v2.5 08-Aug-2010 20:35:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TraceInspector_OpeningFcn, ...
                   'gui_OutputFcn',  @TraceInspector_OutputFcn, ...
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


% --- Executes just before TraceInspector is made visible.
function TraceInspector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TraceInspector (see VARARGIN)

% Choose default command line output for TraceInspector
handles.output = hObject;

% Begin code by Matteo
SetDefaultDirs;
p = varargin{1};
set(handles.txtAnimal,'String',p.animal);
set(handles.txtSeries,'String',p.iseries);
set(handles.txtExperiment,'String',p.iexp);

handles.expt = ExptLoad(p);

if isempty(handles.expt)
    close(hObject);
    return;
end

nstim    = handles.expt.nstim;
nrepeats = handles.expt.nrepeats;

set(handles.lbChannel,'String',num2cell(handles.expt.channames));
set(handles.lbStim   ,'String',num2cell(1:nstim));
set(handles.lbRepeat ,'String',num2cell(1:nrepeats));

set(handles.editGain  ,'String',10);
set(handles.editOffset,'String', 0);

handles = TraceInspector('GetTraces',handles.figure1,handles);
% this fills handles.traces

handles = TraceInspector('ShowTraces',handles.figure1,handles);
% this shows handles.traces

% End code by Matteo

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TraceInspector wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TraceInspector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% MC corrected this as the figure may already hav closed...
if isempty(handles)
    varargout{1} = [];
else
    % Get default command line output from handles structure
    varargout{1} = handles.output;
end

% --- Executes on selection change in lbStim.
function lbStim_Callback(hObject, eventdata, handles)
% hObject    handle to lbStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lbStim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbStim

TraceInspector('ShowTraces',handles.figure1,handles);
% this shows handles.traces

% --- Executes during object creation, after setting all properties.
function lbStim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in lbChannel.
function lbRepeat_Callback(hObject, eventdata, handles)
% hObject    handle to lbRepeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lbRepeat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbRepeat

handles = TraceInspector('GetTraces',handles.figure1,handles);
% this fills handles.traces

handles = TraceInspector('ShowTraces',handles.figure1,handles);
% this shows handles.traces


% --- Executes during object creation, after setting all properties.
function lbRepeat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbRepeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lbChannel.
function lbChannel_Callback(hObject, eventdata, handles)
% hObject    handle to lbChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lbChannel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbChannel

handles = TraceInspector('GetTraces',handles.figure1,handles);
% this fills handles.traces

TraceInspector('ShowTraces',handles.figure1,handles);
% this shows handles.traces

% --- Executes during object creation, after setting all properties.
function lbChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editGain_Callback(hObject, eventdata, handles)
% hObject    handle to editGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGain as text
%        str2double(get(hObject,'String')) returns contents of editGain as a double

NewGain = str2double(get(hObject,'String'));
if ~isnumeric(NewGain) || isnan(NewGain) || isinf(NewGain)
    set(hObject,'String',10);
    return;
end

handles = TraceInspector('GetTraces',handles.figure1,handles);
% this fills handles.traces

handles = TraceInspector('ShowTraces',handles.figure1,handles);
% this shows handles.traces
    


% --- Executes during object creation, after setting all properties.
function editGain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editOffset_Callback(hObject, eventdata, handles)
% hObject    handle to editOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editOffset as text
%        str2double(get(hObject,'String')) returns contents of editOffset as a double

NewOffset = str2double(get(hObject,'String'));
if ~isnumeric(NewOffset) || isnan(NewOffset) || isinf(NewOffset)
    set(hObject,'String', 0);
    return;
end

handles = TraceInspector('GetTraces',handles.figure1,handles);
% this fills handles.traces

handles = TraceInspector('ShowTraces',handles.figure1,handles);
% this shows handles.traces


% --- Executes during object creation, after setting all properties.
function editOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbSave.
function pbSave_Callback(hObject, eventdata, handles)
% hObject    handle to pbSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SetDefaultDirs;
unit = UnitCreate( handles.expt );
unit.traces = handles.traces;
unit.sampledur = 1/(handles.expt.samplerate);
unit.ichan = get(handles.lbChannel,'Value');
unit.datatype = 'traces';
unit.icell    = 'traces';
UnitSave(unit,DIRS.spikes);

% send a wake-up call to FigPicker:
FigPicker_callbacks LoadUnits;



% --- Written by Matteo
function handles = GetTraces( hObject, handles)
% hObject    handle to figure (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% this function fills handles.traces

nstim    = handles.expt.nstim;
nrepeats = handles.expt.nrepeats;

gain    = str2double(get(handles.editGain  ,'String'));
offset  = str2double(get(handles.editOffset,'String'));

ichan = get(handles.lbChannel,'Value');

% load the traces (depends on gain and offset)
handles.traces = cell(nstim,nrepeats);
lims = [Inf -Inf];
for istim = 1:nstim
    for irpt = 1:nrepeats
        
        if ~isempty(handles.expt.data{ichan}{istim, irpt})
            vv = handles.expt.data{ichan}{istim, irpt};
            vv = (double(vv)*3276.8)/1000; 
            vv = ((vv*gain) + offset);
            handles.traces{istim,irpt} = vv;
            lims = [min(lims(1),min(vv)) max(lims(2),max(vv))];
        end
    end
end
handles.lims = lims;

% Update handles structure
guidata(hObject, handles);


% --- Written by Matteo
function handles = ShowTraces( hObject, handles)
% hObject    handle to figure (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% this function shows handles.traces

% show the traces
istim = get(handles.lbStim,'Value');
irpt  = get(handles.lbRepeat,'Value');
vv = handles.traces{istim,irpt};

if isempty(vv)
    cla(handles.axTrace);
else
    tt = (1:length(vv))/(handles.expt.samplerate);
    plot( handles.axTrace,tt, vv );
    set( handles.axTrace,'xlim',[-Inf Inf],'ylim',handles.lims );
end
