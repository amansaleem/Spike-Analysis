function varargout = ProtocolInspect(varargin)
% PROTOCOLINSPECT shows a table of the parameters values in a protocol
%
% ProtocolInspect inspects the protocol in the global PICK.protocol
%
% ProtocolInspect(p) inspects protocol p
%
% ProtocolInspect(animal,iseries,iexp)
%
% EXAMPLE:
% SetDefaultDirs; 
% ProtocolInspect('catz083',1,1)
%

% 2009-03 Matteo Carandini

global PICK


    
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProtocolInspect_OpeningFcn, ...
                   'gui_OutputFcn',  @ProtocolInspect_OutputFcn, ...
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


% --- Executes just before ProtocolInspect is made visible.
function ProtocolInspect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProtocolInspect (see VARARGIN)

% Choose default command line output for ProtocolInspect
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ProtocolInspect wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% added by MC from here
global PICK
switch length(varargin)
    case 0
        if ~isempty(PICK) && isfield(PICK,'protocol')
            p = PICK.protocol;
        end
    case 1
        p = varargin{1};
    case 3
        animal = varargin{1};
        iseries = varargin{2};
        iexp = varargin{3};
        p = ProtocolLoad(animal,iseries,iexp);
    otherwise
        error('Do not understand arguments');
end

set(handles.tblParameters,'Data',p.pars');
set(handles.tblParameters,'ColumnName',p.parnames);
set(hObject,'UserData',p);
set(handles.edit1,'String','Click on a parameter to see definition');
set(handles.txtAnimal,'String',p.animal);
set(handles.txtSeries,'String',p.iseries);
set(handles.txtExperiment,'String',p.iexp);
% ---- to here  

% --- Outputs from this function are returned to the command line.
function varargout = ProtocolInspect_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function tblParameters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tblParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected cell(s) is changed in tblParameters.
function tblParameters_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to tblParameters (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

ipar = eventdata.Indices(2);
p = get(handles.figure1,'UserData');
set(handles.edit1,'String',p.pardefs{ipar})



function txtAnimal_Callback(hObject, eventdata, handles)
% hObject    handle to txtAnimal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAnimal as text
%        str2double(get(hObject,'String')) returns contents of txtAnimal as a double


% --- Executes during object creation, after setting all properties.
function txtAnimal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAnimal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtSeries_Callback(hObject, eventdata, handles)
% hObject    handle to txtSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtSeries as text
%        str2double(get(hObject,'String')) returns contents of txtSeries as a double


% --- Executes during object creation, after setting all properties.
function txtSeries_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtExperiment_Callback(hObject, eventdata, handles)
% hObject    handle to txtExperiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtExperiment as text
%        str2double(get(hObject,'String')) returns contents of txtExperiment as a double


% --- Executes during object creation, after setting all properties.
function txtExperiment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtExperiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


