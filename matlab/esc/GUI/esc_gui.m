function varargout = esc_gui(varargin)
% ESC_GUI MATLAB code for esc_gui.fig
%      ESC_GUI, by itself, creates a new ESC_GUI or raises the existing
%      singleton*.
%
%      H = ESC_GUI returns the handle to a new ESC_GUI or the handle to
%      the existing singleton*.
%
%      ESC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ESC_GUI.M with the given input arguments.
%
%      ESC_GUI('Property','Value',...) creates a new ESC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before esc_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to esc_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help esc_gui

% Last Modified by GUIDE v2.5 07-Jan-2012 22:28:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @esc_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @esc_gui_OutputFcn, ...
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


% --- Executes just before esc_gui is made visible.
function esc_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to esc_gui (see VARARGIN)

% Choose default command line output for esc_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using esc_gui.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(1));
	 cla;
end

%Find available serial ports and display them
pushbuttonRefreshSerialPorts_Callback(hObject, eventdata, handles);

% UIWAIT makes esc_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%Add paths for backend functions
rootPath=pwd;
path(path,fullfile(rootPath, '../'));

%Create esc serial object
global esc
if isempty(esc)
	esc = EscSerial;
end

if ispc
	set(handles.editLogFileName, 'String', 'C:\foo\bar\')
else
	set(handles.editLogFileName, 'String', '/foo/bar/')
end


% --- Outputs from this function are returned to the command line.
function varargout = esc_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end


% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)


% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenuSerialConnection.
function popupmenuSerialConnection_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSerialConnection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuSerialConnection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuSerialConnection


% --- Executes during object creation, after setting all properties.
function popupmenuSerialConnection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSerialConnection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxTransmitCommand.
function checkboxTransmitCommand_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTransmitCommand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTransmitCommand


% --- Executes on button press in pushbuttonRefreshSerialPorts.
function pushbuttonRefreshSerialPorts_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRefreshSerialPorts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

serialInfo = instrhwinfo('serial');
if ~isempty(serialInfo.SerialPorts)
	set(handles.popupmenuSerialConnection, 'String', serialInfo.AvailableSerialPorts)
else
	warndlg('No serial ports found')
end

% --- Executes on button press in buttonSerialConnection.
function buttonSerialConnection_Callback(hObject, eventdata, handles)
% hObject    handle to buttonSerialConnection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Find all serial ports
serialInfo = instrhwinfo('serial');

%Get GUI information
connectionState=get(handles.buttonSerialConnection, 'String');
popup_sel_index = get(handles.popupmenuSerialConnection, 'Value');
serialPortNames=get(handles.popupmenuSerialConnection, 'String');

global esc

if strcmpi(connectionState, 'Connect')
	%Open the selected serial port
	try
		esc = openPort(esc, serialPortNames{popup_sel_index});
	catch SER_ERR
		if strfind(SER_ERR.message, [serialPortNames{popup_sel_index} ' is not available'])
			%%
			choice = questdlg('The selected serial port is already open inside Matlab. Would you like to force close it? (Currently, this forces closed ALL serial ports.)', 'Serial error', 'Yes', 'No', 'Yes');
			if strcmpi(choice, 'yes')
				fclose(instrfind);
				esc = openPort(esc, serialPortNames{popup_sel_index});		
			elseif strcmpi(choice, 'no')
				error('Serial port is already open.')
			end
		end
	end
		

%  	s=serial(serialPortNames{popup_sel_index});
%  	fopen(s);
	
	if strcmp(esc.ser.Status, 'open')
		%Port opened, change button text...
		set(handles.buttonSerialConnection, 'String', 'Disconnect')
	else
		%...it did not open, warn user
		warndlg('Serial port failed to open');
	end
elseif strcmpi(connectionState, 'Disconnect')
	if ~strcmpi(get(handles.buttonRun, 'String'), 'Run')
		errordlg('ESC must be stopped first');		
	else
		esc=closePort(esc);
		
		%Close *all* serial ports. NOT THE BEST WAY TO DO THIS
		% 	fclose(instrfind);
		% 	fclose(serialInfo.AvailableSerialPorts{popup_sel_index});
		
		%Change button text
		set(handles.buttonSerialConnection, 'String', 'Connect')
	end
end


% --- Executes on slider movement.
function sliderAccel_Callback(hObject, eventdata, handles)
% hObject    handle to sliderAccel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=get(handles.sliderAccel, 'Value');
set(handles.editAccel, 'String', [num2str(val, '%0.2f') ' ']);



% --- Executes during object creation, after setting all properties.
function sliderAccel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderAccel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function editAccel_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

val=str2num(get(handles.editAccel, 'String'));
set(handles.sliderAccel, 'Value', val);


% --- Executes during object creation, after setting all properties.
function editAccel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function sliderPhase_Callback(hObject, eventdata, handles)
% hObject    handle to sliderPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val=get(handles.sliderPhase, 'Value');
set(handles.editPhase, 'String', [num2str(val, '%0.2f') ' ']);


% --- Executes during object creation, after setting all properties.
function sliderPhase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function editPhase_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

val=str2num(get(handles.editPhase, 'String'));
set(handles.sliderPhase, 'Value', val);



% --- Executes during object creation, after setting all properties.
function editPhase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuRunType.
function popupmenuRunType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuRunType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuRunType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuRunType


% --- Executes during object creation, after setting all properties.
function popupmenuRunType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuRunType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttonRun.
function buttonRun_Callback(hObject, eventdata, handles)
% hObject    handle to buttonRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global esc

% if portIsOpen
portState=get(handles.buttonSerialConnection, 'String');

if strcmpi(portState, 'Disconnect') %Serial connection says "Disconnect" only when program has successfully connected to serial port
	runState=get(handles.buttonRun, 'String');
	popup_sel_index = get(handles.popupmenuRunType, 'Value');
	% runTypeNames=get(handles.popupmenuRunType, 'String');
	
	
	if strcmpi(runState, 'Run')
		set(handles.buttonRun, 'Value', 1)
		set(handles.buttonRun, 'String', 'Stop')
		
		switch popup_sel_index
			case 1				
				% flush old data and enable logging
				esc = disableLogging(esc);
				if(get(esc.ser,'BytesAvailable') > 0)
					fread(esc.ser, get(esc.ser,'BytesAvailable'), 'uint8');
				end
				
				esc = enableLogging(esc);
				esc.packet = [];
				
				%Basic data gathering
				esc_gather_data(handles)
			case 2
				computeKff
			case 3
				StepResponse
			case 4
				SweepKp
			case 5
				visualizeAdc
		end
		
		% 	eval([runTypeNames{popup_sel_index} ';'])
	elseif strcmpi(runState, 'Stop')
		set(handles.buttonRun, 'Value', 0)
		set(handles.buttonRun, 'String', 'Run')
		esc = setSerialSpeed(esc, 0);
		esc = disableLogging(esc);
	end
	
else
	errordlg('Serial connection must be open first');
end



function editLogFileName_Callback(hObject, eventdata, handles)
% hObject    handle to editLogFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLogFileName as text
%        str2double(get(hObject,'String')) returns contents of editLogFileName as a double


% --- Executes during object creation, after setting all properties.
function editLogFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLogFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxEnableLogging.
function checkboxEnableLogging_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxEnableLogging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxEnableLogging

logOnOff=get(handles.checkboxEnableLogging, 'Value');
if(logOnOff)
	[pathstr, name, ext]=fileparts(get(handles.editLogFileName, 'String'))
	if ~isempty(name)
		warndlg('Only file directories can be set here. The file name is generated as a function of the time and date.')
		set(handles.editLogFileName, 'String', [pathstr '/']); %Remove filename from path
	end
else
end
