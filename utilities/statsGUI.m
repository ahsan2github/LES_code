function varargout = statsGUI(varargin)
% STATSGUI M-file for statsGUI.fig
%      STATSGUI, by itself, creates a new STATSGUI or raises the existing
%      singleton*.
%
%      H = STATSGUI returns the handle to a new STATSGUI or the handle to
%      the existing singleton*.
%
%      STATSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STATSGUI.M with the given input arguments.
%
%      STATSGUI('Property','Value',...) creates a new STATSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before statsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to statsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help statsGUI

% Last Modified by GUIDE v2.5 16-Sep-2014 10:38:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @statsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @statsGUI_OutputFcn, ...
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


% --- Executes just before statsGUI is made visible.
function statsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to statsGUI (see VARARGIN)

% Choose default command line output for statsGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

inputs.Ts=1;
set(handles.inputsPanel,'Userdata',inputs)

updateInputs(handles);

updatePlot(handles)


% UIWAIT makes statsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = statsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Browse.
function Browse_Callback(hObject, eventdata, handles)
% hObject    handle to Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

directoryname=get(handles.simDirectory,'String');

directoryname = uigetdir(directoryname, 'Pick a Directory');

%check that directory was input with trailing slash
if(strcmp(directoryname(end),'/')~=1)
    directoryname=[directoryname,'/']; %if not, fix it
end

%check that directory exists
if(exist(directoryname,'dir')~=7)
    error('%s%s%s','Error: Directory ',directoryname,' does not exist.')
end

%check that LESinputs.txt file exists
if(exist([directoryname,'input/LESinputs.txt'],'file')~=2)
    error('%s%s','Error: LESinputs.txt file does not exist in directory ',...
        directoryname)
end

set(handles.simDirectory,'String',directoryname);

updateInputs(handles);


function simDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to simDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of simDirectory as text
%        str2double(get(hObject,'String')) returns contents of simDirectory as a double

directoryname = get(hObject,'String');

%check that directory was input with trailing slash
if(strcmp(directoryname(end),'/')~=1)
    directoryname=[directoryname,'/']; %if not, fix it
end

%check that directory exists
if(exist(directoryname,'dir')~=7)
    error('%s%s%s','Error: Directory ',directoryname,' does not exist.')
end

%check that LESinputs.txt file exists
if(exist([directoryname,'LESinputs.txt'],'file')~=2)
    error('%s%s','Error: LESinputs.txt file does not exist in directory ',...
        directoryname)
end

rmpath([get(handles.simDirectory,'String'),'utilities'])

set(handles.simDirectory,'String',directoryname);

addpath([directoryname,'utilities'])

inputs=get(handles.inputsPanel,'UserData');
inputs.Ts=1;
set(handles.inputsPanel,'UserData',inputs);
set(handles.Ts_box,'String','1');

updatePlot(handles)


% --- Executes during object creation, after setting all properties.
function simDirectory_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

if(exist('../utilities','dir'))
    addpath('../utilities')
end


function Nx_box_Callback(hObject, eventdata, handles)
% hObject    handle to Nx_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nx_box as text
%        str2double(get(hObject,'String')) returns contents of Nx_box as a double

inputs.Nx=str2num(get(handles.Nx_box,'String'));

set(handles.inputsPanel,'Userdata',inputs)


% --- Executes during object creation, after setting all properties.
function Nx_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nx_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ny_box_Callback(hObject, eventdata, handles)
% hObject    handle to Ny_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ny_box as text
%        str2double(get(hObject,'String')) returns contents of Ny_box as a double

inputs.Ny=str2num(get(handles.Ny_box,'String'));

set(handles.inputsPanel,'Userdata',inputs)


% --- Executes during object creation, after setting all properties.
function Ny_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ny_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nz_box_Callback(hObject, eventdata, handles)
% hObject    handle to Nz_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nz_box as text
%        str2double(get(hObject,'String')) returns contents of Nz_box as a double

inputs.Nz=str2num(get(handles.Nz_box,'String'));

set(handles.inputsPanel,'Userdata',inputs)


% --- Executes during object creation, after setting all properties.
function Nz_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nz_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function zi_box_Callback(hObject, eventdata, handles)
% hObject    handle to zi_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zi_box as text
%        str2double(get(hObject,'String')) returns contents of zi_box as a double


% --- Executes during object creation, after setting all properties.
function zi_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zi_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lz_box_Callback(hObject, eventdata, handles)
% hObject    handle to Lz_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lz_box as text
%        str2double(get(hObject,'String')) returns contents of Lz_box as a double


% --- Executes during object creation, after setting all properties.
function Lz_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lz_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pcount_box_Callback(hObject, eventdata, handles)
% hObject    handle to pcount_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pcount_box as text
%        str2double(get(hObject,'String')) returns contents of pcount_box as a double


inputs.p_count=str2num(get(handles.pcount_box,'String'));

set(handles.inputsPanel,'Userdata',inputs)

% --- Executes during object creation, after setting all properties.
function pcount_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pcount_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ustar_box_Callback(hObject, eventdata, handles)
% hObject    handle to ustar_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ustar_box as text
%        str2double(get(hObject,'String')) returns contents of ustar_box as a double

inputs.u_star=str2num(get(handles.ustar_box,'String'));

set(handles.inputsPanel,'Userdata',inputs)


% --- Executes during object creation, after setting all properties.
function ustar_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ustar_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Updates 'Input Parameters' panel
function updateInputs(handles)

inputs=get(handles.inputsPanel,'Userdata');
directoryname=get(handles.simDirectory,'String');

if(isempty(directoryname))
    error('??? You must slect a directory.')
end

Nx=readFile('Nx',directoryname);
set(handles.Nx_box,'string',Nx)
inputs.Nx=Nx;

Ny=readFile('Ny',directoryname);
set(handles.Ny_box,'string',Ny)
inputs.Ny=Ny;

Nz=readFile('Nz',directoryname);
set(handles.Nz_box,'string',Nz)
inputs.Nz=Nz;

z_i=readFile('z_i',directoryname);
set(handles.zi_box,'string',z_i)
inputs.z_i=z_i;

l_z=readFile('l_z',directoryname);
set(handles.Lz_box,'string',l_z)
inputs.l_z=l_z;

p_count=readFile('p_count',directoryname);
set(handles.pcount_box,'string',p_count)
inputs.p_count=p_count;

u_star=readFile('u_star',directoryname);
set(handles.ustar_box,'string',u_star)
inputs.u_star=u_star;

mom_nodes=readFile('mom_nodes',directoryname);
inputs.mom_nodes=mom_nodes;

inputs.dz=l_z/(Nz-1);

null=loadbin([directoryname,'output/au.bin'],inputs.Nz,'l');
inputs.Te=size(null,1)*inputs.p_count;
set(handles.Te_box,'string',inputs.Te);

set(handles.inputsPanel,'Userdata',inputs);

function result=readFile(name,SA)
str=['grep ''= ',name,''' ',SA,'input/LESinputs.txt | cut -f1 -d''='''];
[status,result]=system(str);
if(status==0)
    %check if result is a double (i.e., contains 'd0')
    if(isempty(strfind(result,'d'))==0)
        result(strfind(result,'d'))='0';
    end
    result=str2double(result);
else
    sprintf('%s%s%s','Warning: Variable "',name,'" was not found')
    result=0;
end

function Ts_box_Callback(hObject, eventdata, handles)
% hObject    handle to Ts_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ts_box as text
%        str2double(get(hObject,'String')) returns contents of Ts_box as a double

inputs=get(handles.inputsPanel,'Userdata');
inputs.Ts=str2double(get(hObject,'String'));
set(handles.inputsPanel,'Userdata',inputs);
updatePlot(handles)

% --- Executes during object creation, after setting all properties.
function Ts_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ts_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Ts_whatsthis.
function Ts_whatsthis_Callback(hObject, eventdata, handles)
% hObject    handle to Ts_whatsthis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox('Ts is the starting timestep for averaging of stats.')



function Te_box_Callback(hObject, eventdata, handles)
% hObject    handle to Te_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Te_box as text
%        str2double(get(hObject,'String')) returns contents of Te_box as a double

inputs=get(handles.inputsPanel,'Userdata');
inputs.Te=str2double(get(hObject,'String'));
set(handles.inputsPanel,'Userdata',inputs);
updatePlot(handles)

% --- Executes during object creation, after setting all properties.
function Te_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Te_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Te_whatsthis.
function Te_whatsthis_Callback(hObject, eventdata, handles)
% hObject    handle to Te_whatsthis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox('Te is the ending timestep for averaging of stats. Default is nsteps.')

function checkinputs(inputs)

if(isfield(inputs,'Nx')==0)
    errordlg('ERROR: Input "Nx" has not been specified')
elseif(isfield(inputs,'Ny')==0)
    errordlg('ERROR: Input "Ny" has not been specified')
elseif(isfield(inputs,'Nz')==0)
    errordlg('ERROR: Input "Nz" has not been specified')
end
    
% --- Executes on selection change in statType.
function statType_Callback(hObject, eventdata, handles)
% hObject    handle to statType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns statType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from statType

updatePlot(handles)


function updatePlot(handles)

stat=get(handles.statType,'Value');

directoryname=get(handles.simDirectory,'String');

inputs=get(handles.inputsPanel,'Userdata');

switch stat
    
    case 1 %check convergence
        data1=loadbin([directoryname,'output/auw.bin'],inputs.Nz,'l')*inputs.u_star^2;
        data2=loadbin([directoryname,'output/atxz.bin'],inputs.Nz,'l')*inputs.u_star^2;
        data=data1+data2;
        int=mean(data,2);
        axis(handles.axes1);
        t=[0:length(int)-1]*inputs.p_count;
        plot(t,int/min(int),'-k')
        xlabel('timesteps','FontSize',14)
        return
    case 2 %au.bin
        data=loadbin([directoryname,'output/au.bin'],inputs.Nz,'l')*inputs.u_star;
        nodes=1;
        xlab='$\langle u \rangle$';
        ylab='$z$';
    case 3 %av.bin
        data=loadbin([directoryname,'output/av.bin'],inputs.Nz,'l')*inputs.u_star;
        nodes=1;
        xlab='$\langle v \rangle$';
        ylab='$z$';
    case 4 %aw.bin
        data=loadbin([directoryname,'output/aw.bin'],inputs.Nz,'l')*inputs.u_star;
        nodes=0;
        xlab='$\langle w \rangle$';
        ylab='$z$';
    case 5 %u2.bin
        data=loadbin([directoryname,'output/u2.bin'],inputs.Nz,'l')*inputs.u_star^2;
        nodes=1;
        xlab='$\langle u^{\prime 2} \rangle$';
        ylab='$z$';
    case 6 %v2.bin
        data=loadbin([directoryname,'output/v2.bin'],inputs.Nz,'l')*inputs.u_star^2;
        nodes=1;
        xlab='$\langle v^{\prime 2} \rangle$';
        ylab='$z$';
    case 7 %w2.bin
        data=loadbin([directoryname,'output/w2.bin'],inputs.Nz,'l')*inputs.u_star^2;
        nodes=0;
        xlab='$\langle w^{\prime 2} \rangle$';
        ylab='$z$';
    case 8 %w3.bin
        data=loadbin([directoryname,'output/w3.bin'],inputs.Nz,'l')*inputs.u_star^3;
        nodes=0;
        xlab='$\langle w^{\prime 3} \rangle$';
        ylab='$z$';
    case 9 %auw.bin+atxz.bin
        data1=loadbin([directoryname,'output/auw.bin'],inputs.Nz,'l')*inputs.u_star^2;
        data2=loadbin([directoryname,'output/atxz.bin'],inputs.Nz,'l')*inputs.u_star^2;
        data=data1+data2;
        nodes=0;
        xlab='$\langle u^\prime w^\prime \rangle$';
        ylab='$z$';
    case 10 %Cs_ALL.bin
        data=loadbin([directoryname,'output/Cs_ALL.bin'],inputs.Nz,'l');
        nodes=inputs.mom_nodes;
        xlab='$\langle C_s \rangle$';
        ylab='$z$';
    case 11 %beta1.bin
        data=loadbin([directoryname,'output/beta1.bin'],inputs.Nz,'l');
        nodes=inputs.mom_nodes;
        xlab='$\langle \beta \rangle$';
        ylab='$z$';
        
    case 12 % spectra
        data = loadbin(directoryname,'output/spectru.bin',inputs.Nz,'l');
        nodes=inputs.mom_nodes;
        xlab='$\langle \beta \rangle$';
        ylab='$z$';
        r = Nx/2.0 + 1; c = size(data, 0) / r;
        data = reshape(data, [r c]);
end

if(nodes==0) %w-nodes
    z=linspace(0,inputs.l_z,inputs.Nz);
else   %u-nodes
    z=linspace(inputs.dz/2,inputs.l_z-inputs.dz/2,inputs.Nz);
end

axes(handles.axes1)

if(inputs.Te/inputs.p_count>size(data,1))
    error('Te is out of range.')
end

data_ave=squeeze(mean(data(ceil(inputs.Ts/inputs.p_count):floor(inputs.Te/inputs.p_count),:)));

plot(data_ave,z,'-k','linewidth',1.5)
xlabel(xlab,'FontSize',14,'Interpreter','Latex');
ylabel(ylab,'FontSize',14,'Interpreter','Latex');
    

% --- Executes during object creation, after setting all properties.
function statType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
