function varargout = mark(varargin)
% MARK MATLAB code for mark.fig
%      MARK, by itself, creates a new MARK or raises the existing
%      singleton*.
%
%      H = MARK returns the handle to a new MARK or the handle to
%      the existing singleton*.
%
%      MARK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARK.M with the given input arguments.
%
%      MARK('Property','Value',...) creates a new MARK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mark_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mark_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mark

% Last Modified by GUIDE v2.5 17-Jul-2017 16:48:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mark_OpeningFcn, ...
                   'gui_OutputFcn',  @mark_OutputFcn, ...
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


% --- Executes just before mark is made visible.
function mark_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mark (see VARARGIN)

% Choose default command line output for mark
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mark wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mark_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_cr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cr as text
%        str2double(get(hObject,'String')) returns contents of edit_cr as a double
global cr
cr = str2double(get(handles.edit_cr, 'String'));
set(handles.slid_cr, 'Value', cr);
draw_images(handles);


% --- Executes during object creation, after setting all properties.
function edit_cr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slid_cr_Callback(hObject, eventdata, handles)
% hObject    handle to slid_cr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global cr
cr = round(get(handles.slid_cr, 'Value'));
if cr < 1, cr = 1; end
set(handles.edit_cr, 'String', num2str(cr));
draw_images(handles);


% --- Executes during object creation, after setting all properties.
function slid_cr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slid_cr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_cc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cc as text
%        str2double(get(hObject,'String')) returns contents of edit_cc as a double
global cc
cc = str2double(get(handles.edit_cc, 'String'));
set(handles.slid_cc, 'Value', cc);
draw_images(handles);


% --- Executes during object creation, after setting all properties.
function edit_cc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slid_cc_Callback(hObject, eventdata, handles)
% hObject    handle to slid_cc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global cc
cc = round(get(handles.slid_cc, 'Value'));
if cc < 1, cc = 1; end
set(handles.edit_cc, 'String', num2str(cc));
draw_images(handles);


% --- Executes during object creation, after setting all properties.
function slid_cc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slid_cc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_rh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rh as text
%        str2double(get(hObject,'String')) returns contents of edit_rh as a double
global rh
rh = str2double(get(handles.edit_rh, 'String'));
set(handles.slid_rh, 'Value', rh);
draw_images(handles);


% --- Executes during object creation, after setting all properties.
function edit_rh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slid_rh_Callback(hObject, eventdata, handles)
% hObject    handle to slid_rh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global rh
rh = round(get(handles.slid_rh, 'Value'));
if rh < 1, rh = 1; end
set(handles.edit_rh, 'String', num2str(rh));
draw_images(handles);


% --- Executes during object creation, after setting all properties.
function slid_rh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slid_rh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_rw_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rw as text
%        str2double(get(hObject,'String')) returns contents of edit_rw as a double
global rw
rw = str2double(get(handles.edit_rw, 'String'));
set(handles.slid_rw, 'Value', rw);
draw_images(handles);


% --- Executes during object creation, after setting all properties.
function edit_rw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slid_rw_Callback(hObject, eventdata, handles)
% hObject    handle to slid_rw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global rw
rw = round(get(handles.slid_rw, 'Value'));
if rw < 1, rw = 1; end
set(handles.edit_rw, 'String', num2str(rw));
draw_images(handles);


% --- Executes during object creation, after setting all properties.
function slid_rw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slid_rw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in button_save.
function button_save_Callback(hObject, eventdata, handles)
% hObject    handle to button_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filename pathname cr cc rh rw
fid = fopen([pathname,'..\groundtruth.txt'], 'r');
if fid < 0
    gt = [];
else
    gt = fscanf(fid, '%d %d %d %d', [4 Inf]);
    fclose(fid);
end
str = strsplit(filename,'.');
file_num = str2num(str{1});
if size(gt,2) >= file_num
    choice = questdlg('此ground truth已存在，是否要取代？','儲存錯誤','是','否','否');
    if strcmp(choice,'否'), return; end
    gt(:,file_num) = [cr;cc;rh;rw];
    fid = fopen([pathname,'..\groundtruth.txt'], 'w');
    for i = 1:size(gt,2)
        fprintf(fid, '%d %d %d %d\n', gt(1,i), gt(2,i), gt(3,i), gt(4,i));
    end
    fclose(fid);
    msgbox('取代儲存成功！', '儲存成功');
else
    fid = fopen([pathname,'..\groundtruth.txt'], 'a');
    for i = size(gt,2)+1:file_num
        fprintf(fid, '%d %d %d %d\n', cr, cc, rh, rw);
    end
    fclose(fid);
    msgbox('儲存成功！', '儲存成功');
end


% --- Executes on button press in button_load.
function button_load_Callback(hObject, eventdata, handles)
% hObject    handle to button_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filename pathname
fn = filename; pn = pathname;
[filename, pathname] = uigetfile({'*.png','*.jpg'},'開啟圖片');
if isequal(filename,0), filename = fn; pathname = pn; return; end
load_image(handles);

% --- Executes on button press in button_save.
function button_last_Callback(hObject, eventdata, handles)
% hObject    handle to button_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filename
str = strsplit(filename,'.');
file_num = str2num(str{1});
if file_num <= 1, return; end
file_num = file_num-1;
filename = [num2str(file_num,'%04d'),'.',str{2}];
load_image(handles);

% --- Executes on button press in button_save.
function button_next_Callback(hObject, eventdata, handles)
% hObject    handle to button_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filename pathname
str = strsplit(filename,'.');
file_num = str2num(str{1});
all_file = dir([pathname,'*.',str{2}]);
if file_num >= length(all_file), return; end
file_num = file_num+1;
filename = [num2str(file_num,'%04d'),'.',str{2}];
load_image(handles);


function [] = load_image(handles)
global filename pathname im ih iw cr cc rh rw
im = imread([pathname,filename]);
[ih,iw,~] = size(im);
fid = fopen([pathname,'..\groundtruth.txt'], 'r');
if fid < 0 % gt.txt not found
    cr = ih/2; cc = iw/2;
    rh = 100; rw = 100;
else
    gt = fscanf(fid, '%d %d %d %d', [4 Inf]);
    fclose(fid);
    str = strsplit(filename,'.');
    file_num = str2num(str{1});
    if size(gt,2) < file_num, file_num = size(gt,2); end
    cr = gt(1,file_num);
    cc = gt(2,file_num);
    rh = gt(3,file_num);
    rw = gt(4,file_num);
end
set(handles.slid_cr, 'Max', ih, 'sliderstep', [1/ih,10/ih]);
set(handles.slid_cr, 'Value', cr);
set(handles.edit_cr, 'String', num2str(cr));
set(handles.slid_cc, 'Max', iw, 'sliderstep', [1/iw,10/iw]);
set(handles.slid_cc, 'Value', cc);
set(handles.edit_cc, 'String', num2str(cc));
set(handles.slid_rh, 'Max', ih, 'sliderstep', [1/ih,10/ih]);
set(handles.slid_rh, 'Value', rh);
set(handles.edit_rh, 'String', num2str(rh));
set(handles.slid_rw, 'Max', iw, 'sliderstep', [1/iw,10/iw]);
set(handles.slid_rw, 'Value', rw);
set(handles.edit_rw, 'String', num2str(rw));
draw_images(handles);
set(handles.text_image, 'String', [pathname,filename]);
set(handles.axes_image.Parent, 'WindowButtonDownFcn', {@mouse_down, handles});

function [] = mouse_down(object, eventdata, handles)
set(handles.axes_image.Parent, 'WindowButtonMotionFcn', {@mouse_move, handles});
set(handles.axes_image.Parent, 'WindowButtonUpFcn', {@mouse_up, handles});
global cr cc ih iw
currPt = get(handles.axes_image, 'CurrentPoint');
if currPt(1,1) < 1 || currPt(1,2) < 1 || currPt(1,1) > iw || currPt(1,2) > ih, return; end
cc = round(currPt(1,1));
cr = round(currPt(1,2));
set(handles.slid_cc, 'Value', cc);
set(handles.edit_cc, 'String', num2str(cc));
set(handles.slid_cr, 'Value', cr);
set(handles.edit_cr, 'String', num2str(cr));
draw_images(handles);

function [] = mouse_move(object, eventdata, handles)
global cr cc ih iw
currPt = get(handles.axes_image, 'CurrentPoint');
if currPt(1,1) < 1 || currPt(1,2) < 1 || currPt(1,1) > iw || currPt(1,2) > ih, return; end
cc = round(currPt(1,1));
cr = round(currPt(1,2));
set(handles.slid_cc, 'Value', cc);
set(handles.edit_cc, 'String', num2str(cc));
set(handles.slid_cr, 'Value', cr);
set(handles.edit_cr, 'String', num2str(cr));
draw_images(handles);

function [] = mouse_up(object, eventdata, handles)
set(handles.axes_image.Parent, 'WindowButtonMotionFcn', '');
set(handles.axes_image.Parent, 'WindowButtonUpFcn', '');

function [] = draw_images(handles)
global im cr cc rh rw ih iw
lines = get_box(im, cr, cc, rh, rw);
axes(handles.axes_image);
imshow(im); axis image
hold on
plot(handles.axes_image, lines(2,:),lines(1,:),'r.')
plot(handles.axes_image, cc,cr,'g.');
hold off
if max(rh/ih,rw/iw) <= 1/3, scale = 1;
elseif max(rh/ih,rw/iw) > 2/3, scale = 4;
else scale = 2;
end
result = get_sceen(im, cr, cc, rh, rw, scale);
axes(handles.axes_result);
imshow(result); axis image
