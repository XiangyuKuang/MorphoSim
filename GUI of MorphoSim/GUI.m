function varargout = MorphoSim(varargin)
         gui_Singleton = 1;
         gui_State = struct('gui_Name'      ,mfilename,       ...
                            'gui_Singleton' ,gui_Singleton,   ...
                            'gui_OpeningFcn',@GUI_OpeningFcn, ...
                            'gui_OutputFcn' ,@GUI_OutputFcn,  ...
                            'gui_LayoutFcn' ,[],              ...
                            'gui_Callback'  ,[])   ;
if  nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1}) ;
end

if  nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State,varargin{:}) ;
else
    gui_mainfcn(gui_State,varargin{:})             ;
end
function GUI_OpeningFcn(hObject,eventdata,handles,varargin)
         handles.output = hObject ; guidata(hObject,handles)    ;
function varargout = GUI_OutputFcn(hObject,eventdata,handles) 
         varargout{1} = handles.output             ;



function Text1_CreateFcn(hObject,eventdata,handles)
function Edit1_Callback( hObject,eventdata,handles)
function Edit1_CreateFcn(hObject,eventdata,handles)
if  ispc && isequal(get(hObject,'BackgroundColor' ),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white') ;
end
function Edit2_CreateFcn(hObject,eventdata,handles)
if  ispc && isequal(get(hObject,'BackgroundColor' ),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white') ;
end
function Edit2_Callback( hObject,eventdata,handles)
function Edit3_Callback( hObject,eventdata,handles)
function Edit3_CreateFcn(hObject,eventdata,handles)
if  ispc && isequal(get(hObject,'BackgroundColor' ),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white') ;
end
function Edit4_Callback( hObject,eventdata,handles)
function Edit4_CreateFcn(hObject,eventdata,handles)
if  ispc && isequal(get(hObject,'BackgroundColor' ),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white') ;
end
function Edit5_Callback( hObject,eventdata,handles)
function Edit5_CreateFcn(hObject,eventdata,handles)
if  ispc && isequal(get(hObject,'BackgroundColor' ),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white') ;
end
function Edit6_Callback( hObject,eventdata,handles)
function Edit6_CreateFcn(hObject,eventdata,handles)
if  ispc && isequal(get(hObject,'BackgroundColor' ),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white') ;
end
function Edit7_Callback( hObject,eventdata,handles)
function Edit7_CreateFcn(hObject,eventdata,handles)
if  ispc && isequal(get(hObject,'BackgroundColor' ),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white') ;
end
function Edit8_Callback( hObject,eventdata,handles)
function Edit8_CreateFcn(hObject,eventdata,handles)
if  ispc && isequal(get(hObject,'BackgroundColor' ),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white') ;
end
function Edit9_Callback( hObject,eventdata,handles)
function Edit9_CreateFcn(hObject,eventdata,handles)
if  ispc && isequal(get(hObject,'BackgroundColor' ),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white') ;
end
function Edit10_Callback( hObject,eventdata,handles)
function Edit10_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'   ), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white') ;
end



function PushButton1_Callback(hObject,eventdata,handles)
    data=importdata(get(handles.Edit1,'string')) ; sigmaIn=data.data ;
    CellName=data.textdata              ; CellName=CellName(2:end,1) ;
    load(get(handles.Edit2,'string'))   ;
    TimeNum=str2num(get(handles.Edit3,'string'))        ;
    dt=str2num(get(handles.Edit10,'string'))            ;
    SavingInterval=str2num(get(handles.Edit4,'string')) ;
    SimName=get(handles.Edit8,'string') ;
if  strcmp(get(handles.Edit7,'string'),'Yes')==1
    EnableEggshell=1  ;
end
if  strcmp(get(handles.Edit7,'string'),'No' )==1
    EnableEggshell=0  ;
end
if  strcmp(get(handles.Edit9,'string'),'GPU')==1
    EnableGPU=1       ;
end
if  strcmp(get(handles.Edit9,'string'),'CPU')==1
    EnableGPU=0       ;
end
    GuiSim(CellName,sigmaIn,Mask,dt,TimeNum,SavingInterval,SimName,EnableEggshell,EnableGPU) ;
    File=dir(SimName) ; fclose('all')   ;
for I=3:length(File)
    Name=File(I).name ;
end
function PushButton2_Callback(hObject, eventdata, handles)
    CellName=importdata(get(handles.Edit5,'string')) ;
    load(get(handles.Edit6,'string'))   ; color = hsv(length(phi_save)) ;
    CellParameters.Name=CellName        ; figure     ;
    MultiphaseDisp(phi_save,CellIdx,0.5,CellParameters,'color',color)