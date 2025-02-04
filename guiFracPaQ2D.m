function varargout = guiFracPaQ2D(varargin)
% GUIFRACPAQ2D MATLAB code for guiFracPaQ2D.fig
%      GUIFRACPAQ2D, by itself, creates a new GUIFRACPAQ2D or raises the existing
%      singleton*.
%
%      H = GUIFRACPAQ2D returns the handle to a new GUIFRACPAQ2D or the handle to
%      the existing singleton*.
%
%      GUIFRACPAQ2D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIFRACPAQ2D.M with the given input arguments.
%
%      GUIFRACPAQ2D('Property','Value',...) creates a new GUIFRACPAQ2D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiFracPaQ2D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiFracPaQ2D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%% Copyright
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.


% Edit the above text to modify the response to help guiFracPaQ2D

% Last Modified by GUIDE v2.5 21-Apr-2021 15:55:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiFracPaQ2D_OpeningFcn, ...
                   'gui_OutputFcn',  @guiFracPaQ2D_OutputFcn, ...
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


%   my initialisation can go here... 

% --- Executes just before guiFracPaQ2D is made visible.
function guiFracPaQ2D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiFracPaQ2D (see VARARGIN)

% Choose default command line output for guiFracPaQ2D
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiFracPaQ2D wait for user response (see UIRESUME)
% uiwait(handles.figure_main);
axes(handles.axes_logo) ;
imLogo = imread('FracPaQlogo.jpeg') ; 
imshow(imLogo) ; 

axes(handles.axes_icon) ;
imIcon = imread('FracPaQicon.jpeg') ; 
imshow(imIcon) ; 

axes(handles.axes_tracemap) ; 

set(handles.popupmenu_histolengthbins,'Value',2) ;
set(handles.popupmenu_histoanglebins,'Value',2) ;
set(handles.popupmenu_rosebins,'Value',2) ;

handles.FracPaQversion = '2.8' ; 
disp(['FracPaQ version ', handles.FracPaQversion]) ; 
handles.graphx = [0, 0] ; 
handles.graphy = [0, 0] ; 
% Update handles structure
guidata(hObject, handles);

v = ver('MATLAB') ; 
iRelease = regexp(v.Release, 'R.....') ; 
nReleaseYear = str2num(v.Release(iRelease+1:iRelease+4)) ; 
if nReleaseYear < 2015 
    h = msgbox('Warning: FracPaQ needs MATLAB release of R2015a, or later.', 'Important Warning!') ;  
    uiwait(h) ; 
end ; 

%Reposition each panel to same location as panel 1
set(handles.uibgTabLengths,'position',get(handles.uibgTabMaps,'position'));
set(handles.uibgTabAngles,'position',get(handles.uibgTabMaps,'position'));
set(handles.uibgTabFluid,'position',get(handles.uibgTabMaps,'position'));
set(handles.uibgTabWavelets,'position',get(handles.uibgTabMaps,'position'));
set(handles.uibgTabGraphs,'position',get(handles.uibgTabMaps,'position'));

set(handles.uibgTabMaps, 'Visible', 'off') ; 
set(handles.uibgTabLengths, 'Visible', 'off') ; 
set(handles.uibgTabAngles, 'Visible', 'off') ; 
set(handles.uibgTabFluid, 'Visible', 'off') ; 
set(handles.uibgTabWavelets, 'Visible', 'off') ; 
set(handles.uibgTabGraphs, 'Visible', 'off') ; 

set(handles.edit_filename, 'Enable', 'on') ; 

% --- Outputs from this function are returned to the command line.
function varargout = guiFracPaQ2D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox_triangle.
function checkbox_triangle_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_triangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_triangle
if get(hObject, 'Value')
    set(handles.label_numblocksmap, 'Enable', 'on') ; 
    set(handles.edit_numblocksmap, 'Enable', 'on') ; 
    set(handles.label_numpixelsItoY, 'Enable', 'on') ; 
    set(handles.edit_numpixelsItoY, 'Enable', 'on') ; 
else 
    set(handles.label_numblocksmap, 'Enable', 'off') ; 
    set(handles.edit_numblocksmap, 'Enable', 'off') ; 
    set(handles.label_numpixelsItoY, 'Enable', 'off') ; 
    set(handles.edit_numpixelsItoY, 'Enable', 'off') ; 
end ; 

% --- Executes on button press in checkbox_permellipse.
function checkbox_permellipse_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_permellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_permellipse
if get(hObject, 'Value')
    set(handles.label_aperturefactor, 'Enable', 'on') ; 
    set(handles.edit_aperturefactor, 'Enable', 'on') ; 
    set(handles.label_apertureexponent, 'Enable', 'on') ; 
    set(handles.edit_apertureexponent, 'Enable', 'on') ; 
    set(handles.edit_lambda, 'Enable', 'on') ; 
    set(handles.label_lambda, 'Enable', 'on') ; 
    set(handles.edit_fixedaperture, 'Enable', 'on') ; 
    set(handles.label_fixedaperture, 'Enable', 'on') ; 
    set(handles.radiobutton_fixedaperture, 'Enable', 'on') ; 
    set(handles.radiobutton_scaledaperture, 'Enable', 'on') ; 
else
    set(handles.label_aperturefactor, 'Enable', 'off') ; 
    set(handles.edit_aperturefactor, 'Enable', 'off') ; 
    set(handles.label_apertureexponent, 'Enable', 'off') ; 
    set(handles.edit_apertureexponent, 'Enable', 'off') ; 
    set(handles.edit_lambda, 'Enable', 'off') ; 
    set(handles.label_lambda, 'Enable', 'off') ; 
    set(handles.edit_fixedaperture, 'Enable', 'off') ; 
    set(handles.label_fixedaperture, 'Enable', 'off') ; 
    set(handles.radiobutton_fixedaperture, 'Enable', 'off') ; 
    set(handles.radiobutton_scaledaperture, 'Enable', 'off') ; 
end ; 

% --- Executes on button press in checkbox_intensitymap.
function checkbox_intensitymap_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_intensitymap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_intensitymap
if get(hObject, 'Value') || get(handles.checkbox_densitymap, 'Value')
    set(handles.checkbox_showcircles, 'Enable', 'on') ; 
    set(handles.edit_numscancircles, 'Enable', 'on') ; 
    set(handles.label_numscancircles, 'Enable', 'on') ; 
    set(handles.edit_circlespacing, 'Enable', 'on') ; 
    set(handles.label_circlespacing, 'Enable', 'on') ; 
    set(handles.checkbox_rotategrid, 'Enable', 'on') ; 
    if get(handles.checkbox_rotategrid, 'Value')
        set(handles.label_gridrotation, 'Enable', 'on') ;     
        set(handles.edit_gridrotation, 'Enable', 'on') ; 
    end;    
else 
    set(handles.checkbox_showcircles, 'Enable', 'off') ; 
    set(handles.checkbox_showcircles, 'Value', 0) ; 
    set(handles.edit_numscancircles, 'Enable', 'off') ; 
    set(handles.edit_numscancircles, 'Value', 0) ; 
    set(handles.label_numscancircles, 'Enable', 'off') ; 
    set(handles.edit_circlespacing, 'Enable', 'off') ; 
    set(handles.edit_circlespacing, 'Value', 0) ; 
    set(handles.label_circlespacing, 'Enable', 'off') ; 
    set(handles.checkbox_rotategrid, 'Enable', 'off') ; 
%    set(handles.checkbox_rotategrid, 'Value', 0) ; 
    set(handles.edit_gridrotation, 'Enable', 'off') ;
%    set(handles.edit_gridrotation, 'Value', 0) ; 
    set(handles.label_gridrotation, 'Enable', 'off') ; 
end ; 


% --- Executes on button press in checkbox_densitymap.
function checkbox_densitymap_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_densitymap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_densitymap
if get(hObject, 'Value') || get(handles.checkbox_intensitymap, 'Value')
    set(handles.checkbox_showcircles, 'Enable', 'on') ; 
    set(handles.edit_numscancircles, 'Enable', 'on') ; 
    set(handles.label_numscancircles, 'Enable', 'on') ; 
    set(handles.edit_circlespacing, 'Enable', 'on') ; 
    set(handles.label_circlespacing, 'Enable', 'on') ; 
    set(handles.checkbox_rotategrid, 'Enable', 'on') ; 
    if get(handles.checkbox_rotategrid, 'Value')
        set(handles.label_gridrotation, 'Enable', 'on') ;     
        set(handles.edit_gridrotation, 'Enable', 'on') ; 
    end;
else 
    set(handles.checkbox_showcircles, 'Enable', 'off') ; 
    set(handles.checkbox_showcircles, 'Value', 0) ; 
    set(handles.edit_numscancircles, 'Enable', 'off') ; 
    set(handles.edit_numscancircles, 'Value', 0) ; 
    set(handles.label_numscancircles, 'Enable', 'off') ; 
    set(handles.edit_circlespacing, 'Enable', 'off') ; 
    set(handles.edit_circlespacing, 'Value', 0) ; 
    set(handles.label_circlespacing, 'Enable', 'off') ; 
    set(handles.checkbox_rotategrid, 'Enable', 'off') ; 
%    set(handles.checkbox_rotategrid, 'Value', 0) ; 
    set(handles.edit_gridrotation, 'Enable', 'off') ; 
%    set(handles.edit_gridrotation, 'Value', 0) ; 
    set(handles.label_gridrotation, 'Enable', 'off') ; 
end ; 


% --- Executes on button press in checkbox_tracemap.
function checkbox_tracemap_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_tracemap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_tracemap
if get(hObject, 'Value') 
    set(handles.checkbox_shownodes, 'Enable', 'on') ; 
else 
    set(handles.checkbox_shownodes, 'Enable', 'off') ; 
    set(handles.checkbox_shownodes, 'Value', 0) ; 
end ; 
    
% --- Executes on button press in checkbox_shownodes.
function checkbox_shownodes_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_shownodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_shownodes


% --- Executes on button press in checkbox_histolength.
function checkbox_histolength_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_histolength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_histolength
if get(hObject, 'Value')
    set(handles.text_histolengthbins, 'Enable', 'on') ; 
    set(handles.popupmenu_histolengthbins, 'Enable', 'on') ; 
else 
    set(handles.text_histolengthbins, 'Enable', 'off') ; 
    set(handles.popupmenu_histolengthbins, 'Enable', 'off') ; 
end ; 


% --- Executes on button press in checkbox_logloglength.
function checkbox_logloglength_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_logloglength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_logloglength
if get(hObject, 'Value')
    set(handles.checkbox_censor, 'Enable', 'on') ; 
    set(handles.checkbox_mle, 'Enable', 'on') ; 
else    
    set(handles.checkbox_censor, 'Enable', 'off') ; 
    set(handles.checkbox_mle, 'Enable', 'off') ; 
    set(handles.checkbox_censor, 'Value', 0) ; 
    set(handles.checkbox_mle, 'Value', 0) ; 
    set(handles.edit_lc, 'Enable', 'off') ; 
    set(handles.edit_uc, 'Enable', 'off') ; 
    set(handles.label_lc, 'Enable', 'off') ; 
    set(handles.label_uc, 'Enable', 'off') ; 
end ; 

% --- Executes on button press in checkbox_mle.
function checkbox_mle_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_mle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_mle
if get(hObject, 'Value')
    set(handles.edit_lc, 'Enable', 'on') ; 
    set(handles.edit_uc, 'Enable', 'on') ; 
    set(handles.label_lc, 'Enable', 'on') ; 
    set(handles.label_uc, 'Enable', 'on') ; 
else
    set(handles.edit_lc, 'Enable', 'off') ; 
    set(handles.edit_uc, 'Enable', 'off') ; 
    set(handles.label_lc, 'Enable', 'off') ; 
    set(handles.label_uc, 'Enable', 'off') ; 
end ; 

% --- Executes on button press in checkbox_histoangle.
function checkbox_histoangle_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_histoangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_histoangle
if get(hObject, 'Value')
    set(handles.text_histoanglebins, 'Enable', 'on') ; 
    set(handles.popupmenu_histoanglebins, 'Enable', 'on') ; 
else 
    set(handles.text_histoanglebins, 'Enable', 'off') ; 
    set(handles.popupmenu_histoanglebins, 'Enable', 'off') ; 
end ; 


% --- Executes on button press in checkbox_rose.
function checkbox_rose_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_rose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_rose
if get(hObject, 'Value')
    set(handles.text_rosebins, 'Enable', 'on') ; 
    set(handles.popupmenu_rosebins, 'Enable', 'on') ; 
    set(handles.edit_degfromnorth, 'Enable', 'on') ; 
    set(handles.label_degfromnorth, 'Enable', 'on') ; 
    set(handles.checkbox_roselengthweighted, 'Enable', 'on') ; 
    set(handles.checkbox_showrosemean, 'Enable', 'on') ; 
    set(handles.checkbox_rosecolour, 'Enable', 'on') ; 
else 
    set(handles.text_rosebins, 'Enable', 'off') ; 
    set(handles.popupmenu_rosebins, 'Enable', 'off') ; 
    set(handles.edit_degfromnorth, 'Enable', 'off') ; 
    set(handles.label_degfromnorth, 'Enable', 'off') ; 
    set(handles.checkbox_roselengthweighted, 'Enable', 'off') ; 
    set(handles.checkbox_showrosemean, 'Enable', 'off') ; 
    set(handles.checkbox_rosecolour, 'Enable', 'off') ; 
    set(handles.checkbox_roselengthweighted, 'Value', false) ; 
    set(handles.checkbox_showrosemean, 'Value', false) ; 
    set(handles.checkbox_rosecolour, 'Value', false) ; 
end ; 


% --- Executes on button press in checkbox_crossplot.
function checkbox_crossplot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_crossplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_crossplot


% --- Executes on button press in checkbox_showcircles.
function checkbox_showcircles_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showcircles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showcircles


% --- Executes on button press in checkbox_censor.
function checkbox_censor_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_censor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_censor



function edit_houghpeaks_Callback(hObject, eventdata, handles)
% hObject    handle to edit_houghpeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_houghpeaks as text
%        str2double(get(hObject,'String')) returns contents of edit_houghpeaks as a double


% --- Executes during object creation, after setting all properties.
function edit_houghpeaks_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_houghpeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_houghthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_houghthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_houghthreshold as text
%        str2double(get(hObject,'String')) returns contents of edit_houghthreshold as a double


% --- Executes during object creation, after setting all properties.
function edit_houghthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_houghthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fillgap_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fillgap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fillgap as text
%        str2double(get(hObject,'String')) returns contents of edit_fillgap as a double


% --- Executes during object creation, after setting all properties.
function edit_fillgap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fillgap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minlength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minlength as text
%        str2double(get(hObject,'String')) returns contents of edit_minlength as a double


% --- Executes during object creation, after setting all properties.
function edit_minlength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_degfromnorth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_degfromnorth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_degfromnorth as text
%        str2double(get(hObject,'String')) returns contents of edit_degfromnorth as a double


% --- Executes during object creation, after setting all properties.
function edit_degfromnorth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_degfromnorth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_scaling_Callback(hObject, eventdata, handles)
% hObject    handle to edit_scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_scaling as text
%        str2double(get(hObject,'String')) returns contents of edit_scaling as a double


% --- Executes during object creation, after setting all properties.
function edit_scaling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_numscancircles_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numscancircles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numscancircles as text
%        str2double(get(hObject,'String')) returns contents of edit_numscancircles as a double
 if length(handles.traces) > 0
        nCount = str2double(get(hObject,'String'));
        if get(handles.checkbox_rotategrid, 'Value')
            rotation = round(str2double(get(handles.edit_gridrotation, 'String'))) ; 
        else
            rotation = 0;
        end
        [dist, iCount, jCount] = circleSpacingFromNumber(handles.traces, rotation, nCount);
        set(handles.edit_circlespacing,'String', num2str(dist,'%.1f'));
    end

% --- Executes during object creation, after setting all properties.
function edit_numscancircles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numscancircles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lambda as text
%        str2double(get(hObject,'String')) returns contents of edit_lambda as a double


% --- Executes during object creation, after setting all properties.
function edit_lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_exit.
function pushbutton_exit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all ; 
close(gcf) ; 


function edit_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filename as text
%        str2double(get(hObject,'String')) returns contents of edit_filename as a double
cla(handles.axes_tracemap) ; 
handles.axes_tracemap.YDir = 'normal' ; 

handles.selfile = get(hObject, 'String') ;
handles.selpath = pwd ; 
if ispc 
    handles.selpath = [handles.selpath, '\'] ;
else
    handles.selpath = [handles.selpath, '/'] ;
end 

if ~isempty(handles.selfile)
    
    sFileName = char(handles.selfile) ; 

    if ~isempty(regexp(sFileName, 'colour', 'ONCE')) 
        iHex = regexp(sFileName, '[A-F,0-9]{6}', 'ONCE') ; 
        %   check each file has 'HHHHHH' in the name i.e. a valid hex string for colour 
        if isempty(iHex) 
            hError = errordlg('For a colour file, all filenames must contain a 6 character hexadecimal colour code', 'Input error', 'modal') ; 
            uicontrol(handles.edit_filename) ; 
            flagError = true ; 
            return ; 
        else 
            iRGB = strfind(sFileName, 'colour') ; 
            sRGB = sFileName(iRGB+6:iRGB+11) ; 
            sColour = hexRGB(sRGB) ; 
            handles.cstrColours = cellstr(sColour) ; 
        end ; 

    else
        %   MATLAB default blue 
        sColour = '[0, 0, 1]' ; 
        handles.cstrColours = cellstr(sColour) ; 
    end ; 

    handles.iCurrentColour = 1 ; 

    set(handles.pushbutton_preview, 'Enable', 'on') ; 

    if ~isempty(strfind(upper(handles.selfile), '.TIF')) || ...
       ~isempty(strfind(upper(handles.selfile), '.TIFF')) || ...
       ~isempty(strfind(upper(handles.selfile), '.JPEG')) || ...
       ~isempty(strfind(upper(handles.selfile), '.JPG'))  
        set(handles.radiobutton_image, 'Value', 1) ;
        set(handles.radiobutton_node, 'Value', 0) ;
        set(handles.edit_houghpeaks, 'Enable', 'on') ; 
        set(handles.edit_houghthreshold, 'Enable', 'on') ; 
        set(handles.edit_fillgap, 'Enable', 'on') ; 
        set(handles.edit_minlength, 'Enable', 'on') ; 
    else 
        set(handles.radiobutton_node, 'Value', 1) ;
        set(handles.radiobutton_image, 'Value', 0) ;
        set(handles.edit_houghpeaks, 'Enable', 'off') ; 
        set(handles.edit_houghthreshold, 'Enable', 'off') ; 
        set(handles.edit_fillgap, 'Enable', 'off') ; 
        set(handles.edit_minlength, 'Enable', 'off') ; 
    end ; 

end 

set(handles.pushbutton_flipx, 'Enable', 'off') ; 
set(handles.pushbutton_flipy, 'Enable', 'off') ; 
set(handles.pushbutton_run, 'Enable', 'off') ; 
set(handles.text_message, 'String', 'Click Preview to view the file contents.') ;                                 

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white') ;
end


% --- Executes on button press in pushbutton_selectgraphendpoints.
function pushbutton_selectgraphendpoints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectgraphendpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[x, y] = ginput(2) ; 
handles.graphx = x ; 
handles.graphy = y ; 

cla(handles.axes_tracemap) ; 

%   allow for multiple files, one for each colour 
if ~isempty(handles.selmultifile)

    hold on ; 
    for i = 1:length(handles.selmultifile) 

        sFilename = [ char(handles.selpath), char(handles.selmultifile(i)) ] ;
        sColour = char(handles.cstrColours(i)) ; 
        [ ~, ~, ~, ~, ~ ] = guiFracPaQ2Dlist(sFilename, handles.nPixels, handles.axes_tracemap, sColour) ; 

    end ; 
    plot(handles.axes_tracemap, x, y, '-sr', 'LineWidth', 1) ;
    hold off ; 

else
    hold on ; 
    [ ~, ~, ~, ~, ~ ] = guiFracPaQ2Dlist(handles.sFilename, handles.nPixels, handles.axes_tracemap, handles.sColour) ; 
    plot(handles.axes_tracemap, x, y, '-sr', 'LineWidth', 1) ;
    hold off ; 
end ; 

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_browse.
function pushbutton_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.axes_tracemap) ; 
handles.axes_tracemap.YDir = 'normal' ; 

%   from v1.8 - all platforms get access to .svg input files 
sFileTypes = {'*.txt'; '*.svg'; '*.jpeg'; '*.jpg'; '*.tiff'; '*.tif'} ; 

[ handles.selfile, handles.selpath ] = uigetfile(sFileTypes, 'Select input file(s) for FracPaQ2D') ; 

if ~isempty(handles.selfile)
    
    sFileName = char(handles.selfile) ; 

    if ~isempty(regexp(sFileName, 'colour', 'ONCE')) 
        iHex = regexp(sFileName, '[A-F,0-9]{6}', 'ONCE') ; 
        %   check each file has 'HHHHHH' in the name i.e. a valid hex string for colour 
        if isempty(iHex) 
            hError = errordlg('For a colour file, all filenames must contain a 6 character hexadecimal colour code', 'Input error', 'modal') ; 
            uicontrol(handles.edit_filename) ; 
            flagError = true ; 
            return ; 
        else 
            iRGB = strfind(sFileName, 'colour') ; 
            sRGB = sFileName(iRGB+6:iRGB+11) ; 
            sColour = hexRGB(sRGB) ; 
            handles.cstrColours = cellstr(sColour) ; 
        end ; 

    else
        %   MATLAB default blue 
        sColour = '[0, 0, 1]' ; 
        handles.cstrColours = cellstr(sColour) ; 
    end ; 

    handles.iCurrentColour = 1 ; 

    set(handles.edit_filename, 'String', sFileName) ; 
    set(handles.pushbutton_preview, 'Enable', 'on') ; 

    if ~isempty(strfind(upper(sFileName), '.TIF')) || ...
       ~isempty(strfind(upper(sFileName), '.TIFF')) || ...
       ~isempty(strfind(upper(sFileName), '.JPEG')) || ...
       ~isempty(strfind(upper(sFileName), '.JPG'))  
        set(handles.radiobutton_image, 'Value', 1) ;
        set(handles.radiobutton_node, 'Value', 0) ;
        set(handles.edit_houghpeaks, 'Enable', 'on') ; 
        set(handles.edit_houghthreshold, 'Enable', 'on') ; 
        set(handles.edit_fillgap, 'Enable', 'on') ; 
        set(handles.edit_minlength, 'Enable', 'on') ; 
    else 
        set(handles.radiobutton_node, 'Value', 1) ;
        set(handles.radiobutton_image, 'Value', 0) ;
        set(handles.edit_houghpeaks, 'Enable', 'off') ; 
        set(handles.edit_houghthreshold, 'Enable', 'off') ; 
        set(handles.edit_fillgap, 'Enable', 'off') ; 
        set(handles.edit_minlength, 'Enable', 'off') ; 
    end ; 
        
else
    
    set(handles.edit_filename, 'String', '(no file selected)') ; 
    set(handles.pushbutton_preview, 'Enable', 'off') ; 
    
end ; 
    
set(handles.pushbutton_flipx, 'Enable', 'off') ; 
set(handles.pushbutton_flipy, 'Enable', 'off') ; 
set(handles.pushbutton_run, 'Enable', 'off') ; 
set(handles.text_message, 'String', 'Click Preview to view the file contents.') ;                                 

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_preview.
function pushbutton_preview_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global sTag ; 

cla(handles.axes_tracemap) ; 

flagError = false ; 

%   image file validation 
if get(handles.radiobutton_image, 'Value')
    
    sValue = get(handles.edit_houghpeaks, 'String') ; 
    if isnan(str2double(sValue)) || str2double(sValue) < 100 || str2double(sValue) > 5000
        hError = errordlg('Hough peaks must be a whole number, 100-5000', 'Input error', 'modal') ; 
        uicontrol(handles.edit_houghpeaks) ; 
        flagError = true ; 
        return ; 
    end ; 

    sValue = get(handles.edit_houghthreshold, 'String') ; 
    if isnan(str2double(sValue)) || str2double(sValue) <= 0.0 || str2double(sValue) > 1.0  
        hError = errordlg('Hough threshold must be a number > 0.0, and <= 1.0', 'Input error', 'modal') ; 
        uicontrol(handles.edit_houghthreshold) ; 
        flagError = true ; 
        return ; 
    end ; 

    sValue = get(handles.edit_fillgap, 'String') ; 
    if isnan(str2double(sValue)) || str2double(sValue) <= 1.0    
        hError = errordlg('Merge gap must be a number of pixels, > 1 and < max. image size', 'Input error', 'modal') ; 
        uicontrol(handles.edit_fillgap) ; 
        flagError = true ; 
        return ; 
    end ; 

    sValue = get(handles.edit_minlength, 'String') ; 
    if isnan(str2double(sValue)) || str2double(sValue) <= 1   
        hError = errordlg('Discard length must be a number of pixels, > 1 and < max. image size', 'Input error', 'modal') ; 
        uicontrol(handles.edit_minlength) ; 
        flagError = true ; 
        return ; 
    end ; 
    
    Iinfo = imfinfo(get(handles.edit_filename, 'String')) ; 
    if Iinfo.BitDepth ~= 8 || ~strcmp('grayscale', Iinfo.ColorType)
        disp('Image file format error - check BitDepth (8) and ColorType(''grayscale'')') ; 
        disp(Iinfo) ; 
        hError = errordlg({'Selected image file is not 8-bit & grayscale';'Check Command Window for file info, BitDepth and ColorType'}, ... 
                                'Input error', 'modal') ; 
        uicontrol(handles.edit_filename) ; 
        flagError = true ; 
        return ; 
    end ; 

    handles.selmultifile = { } ; 

%   node file validation 
else 
    
    sFilename = [ char(handles.selpath), char(handles.selfile) ] ; 
    disp(sFilename);
    %   convert .svg file to .txt file 
    if ~isempty(strfind(sFilename, '.svg'))  
        
        nConvert = convertSVG2txt_colour2(sFilename) ; 
        
        if nConvert > 0 
            
            %   now handle possible multiple files 
            if nConvert > 1 
                
                sFilenameHead = sFilename(1:strfind(sFilename, '.svg')-1) ; 
                infoConvertedFiles = dir([sFilenameHead, '_colour*_converted.txt']) ;
                handles.selmultifile = { infoConvertedFiles.name } ;
                sColour = '' ;
                for i = 1:nConvert 
                    
                    sThisFileName = char(handles.selmultifile(i)) ; 
                    iHex = regexp(sThisFileName, '[A-F,0-9]{6}', 'ONCE') ; 
                    
%                     disp(sThisFileName) ; 
%                     disp(iHex) ; 
                     
                    %   check each file has 'HHHHHH' in the name i.e. a valid hex string for colour 
                    if isempty(iHex) 
                        hError = errordlg('For a colour file, all filenames must contain a 6 character hexadecimal colour code', 'Input error', 'modal') ; 
                        uicontrol(handles.edit_filename) ; 
                        flagError = true ; 
                        return ; 
                    else 
                        iRGB = strfind(sThisFileName, 'colour') ; 
                        sRGB = sThisFileName(iRGB+6:iRGB+11) ; 
                        sColour = [ sColour ; hexRGB(sRGB) ]  ; 
                    end ; 

                end ; 
                
                handles.cstrColours = cellstr(sColour) ; 
                    
            else
                
                handles.selmultifile = { } ; 
                handles.selfile = strrep(handles.selfile, '.svg', '_converted.txt') ; 
                
            end ; 
            
        else
            
            hError = errordlg('Error converting .svg file to .txt; check .svg file is version 1.1', 'Input error', 'modal') ; 
            uicontrol(handles.edit_filename) ; 
            flagError = true ; 
            
        end ; 
        
    else
        
        handles.selmultifile = { } ; 

    end ; 
    
end ; 

sValue = get(handles.edit_scaling, 'String') ; 
if ~isempty(sValue)
    if isnan(str2double(sValue)) || str2double(sValue) <= 0.0 
        hError = errordlg('Scaling can be blank (i.e. scale as is) or a positive number', 'Input error', 'modal') ; 
        uicontrol(handles.edit_scaling) ; 
        flagError = true ; 
        return ; 
    end ; 
end ; 

if ~flagError 

    sFilename = [ char(handles.selpath), char(handles.selfile) ] ; 
    %   get the current RGB colour string 
    sColour = char(handles.cstrColours(handles.iCurrentColour)) ; 

    if isempty(get(handles.edit_scaling, 'String'))
        nPixels = 0 ; 
    else
        nPixels = str2double(get(handles.edit_scaling, 'String')) ; 
    end ; 
    
    handles.nPixels = nPixels ; 
    handles.sColour = sColour ; 
    handles.sFilename = sFilename ; 
    
    if get(handles.radiobutton_image, 'Value') 
        
        nHoughPeaks = str2double(get(handles.edit_houghpeaks, 'String')) ;
        dHoughThreshold = str2double(get(handles.edit_houghthreshold, 'String')) ;
        dfillgap = str2double(get(handles.edit_fillgap, 'String')) ;
        dminlength = str2double(get(handles.edit_minlength, 'String')) ;

        set(handles.text_message, 'String', 'Running Hough transform on image...') ;  
        [ handles.traces, handles.xmin, handles.ymin, handles.xmax, handles.ymax ] = guiFracPaQ2Dimage(sFilename, nPixels, ...
                                                            nHoughPeaks, dHoughThreshold, dfillgap, dminlength, handles.axes_tracemap) ; 
    else
        
        set(handles.text_message, 'String', 'Reading node file...') ;                                 

        %   allow for multiple files, one for each colour 
        if ~isempty(handles.selmultifile)

            handles.traces = struct([]) ; 
            handles.nTraces = 0 ; 
            handles.nSegments = 0 ; 
            handles.nNodes = 0 ;
            handles.xmin = 3e38 ; 
            handles.ymin = 3e38 ; 
            handles.xmax = -3e38 ;
            handles.ymax = -3e38 ; 
            hold on ; 
            for i = 1:length(handles.selmultifile) 
                
                sFilename = [ char(handles.selpath), char(handles.selmultifile(i)) ] ;
                sColour = char(handles.cstrColours(i)) ; 
                
                [ traces, xmin, ymin, xmax, ymax ] = guiFracPaQ2Dlist(sFilename, nPixels, handles.axes_tracemap, sColour) ; 
                
                handles.traces = [ handles.traces, traces ] ; 

                if xmin < handles.xmin
                    handles.xmin = xmin ; 
                end ; 
                if ymin < handles.ymin 
                    handles.ymin = ymin ; 
                end ; 
                if xmax > handles.xmax
                    handles.xmax = xmax ; 
                end ; 
                if ymax > handles.ymax 
                    handles.ymax = ymax ; 
                end ; 
                
            end ; 
            hold off ; 
           
        else
            
            hold on ; 
            [ handles.traces, handles.xmin, handles.ymin, handles.xmax, handles.ymax ] = guiFracPaQ2Dlist(sFilename, nPixels, handles.axes_tracemap, sColour) ; 
            hold off ; 
            
        end ; 
        
    end ; 
    
    set(handles.text_message, 'String', 'Collecting fracture trace data...') ;     

    setappdata(handles.pushbutton_run, 'Traces', handles.traces) ;                                                 
    setappdata(handles.pushbutton_run, 'xMax', handles.xmax) ;                                                 
    setappdata(handles.pushbutton_run, 'yMax', handles.ymax) ;                                                 
    setappdata(handles.pushbutton_run, 'xMin', handles.xmin) ;                                                 
    setappdata(handles.pushbutton_run, 'yMin', handles.ymin) ;                                                 
    
    sTag = ['_', get(handles.edit_FilenameTag, 'String')] ; 
    
    handles.nTraces = length(handles.traces) ; 
    handles.nSegments = sum([handles.traces(:).nSegments]) ; 
    handles.nNodes = sum([handles.traces(:).nNodes]) ; 

    set(handles.label_stats, 'String', {['Min. X coordinate: ', num2str(handles.xmin)], ...
                                        ['Min. Y coordinate: ', num2str(handles.ymin)], ...
                                        ['Max. X coordinate: ', num2str(handles.xmax)], ...
                                        ['Max. Y coordinate: ', num2str(handles.ymax)], ...
                                        ['Number of traces: ', num2str(handles.nTraces)], ...
                                        ['Number of segments: ', num2str(handles.nSegments)], ...
                                        ['Number of nodes: ', num2str(handles.nNodes)]}) ;

    xyMaxRange = max([handles.xmax-handles.xmin, handles.ymax-handles.ymin]) ;
    set(handles.edit_numpixelsItoY, 'String', num2str(xyMaxRange/1000, '%2.2f')) ; 
    
    set(handles.text_message, 'String', 'Ready. Click Run to generate maps and graphs.') ;     
  
    set(handles.pushbutton_run, 'Enable', 'on') ; 
    set(handles.pbMapsTab, 'Enable', 'on') ; 
    set(handles.pbLengthsTab, 'Enable', 'on') ; 
    set(handles.pbAnglesTab, 'Enable', 'on') ; 
    set(handles.pbFluidTab, 'Enable', 'on') ; 
    set(handles.pbWaveletsTab, 'Enable', 'on') ; 
    set(handles.pbGraphsTab, 'Enable', 'on') ; 
    set(handles.pushbutton_flipx, 'Enable', 'on') ; 
    set(handles.pushbutton_flipy, 'Enable', 'on') ; 
%    uicontrol(handles.pushbutton_run) ; 
    uicontrol(handles.pbMapsTab) ; 
    set(handles.uibgTabMaps, 'Visible', 'on') ; 
    
    % Update handles structure
    guidata(hObject, handles);
    
end ; 

% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

flagError = false ; 
myEps = 1e-8 ; 

sValue = get(handles.edit_sigma1, 'String') ; 
if isnan(str2double(sValue)) || abs(str2double(sValue)) < myEps 
    hError = errordlg('Sigma 1 must be a number, and not 0', 'Input error', 'modal') ; 
    uicontrol(handles.edit_sigma1) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_sigma2, 'String') ; 
if isnan(str2double(sValue)) || abs(str2double(sValue)) < myEps 
    hError = errordlg('Sigma 2 must be a number, and not 0', 'Input error', 'modal') ; 
    uicontrol(handles.edit_sigma2) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_thetasigma1, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) > 180.0 || str2double(sValue) < 0.0
    hError = errordlg('Angle must be a number and between 0-180', 'Input error', 'modal') ; 
    uicontrol(handles.edit_thetasigma1) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_C0, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 0.0
    hError = errordlg('Cohesion (C0) must be a positive number or 0.0', 'Input error', 'modal') ; 
    uicontrol(handles.edit_C0) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_Pf, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 0.0
    hError = errordlg('Pore pressure must be a positive number or 0.0', 'Input error', 'modal') ; 
    uicontrol(handles.edit_Pf) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_mu, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 0.0
    hError = errordlg('Friction coefficient must be a positive number or 0.0', 'Input error', 'modal') ; 
    uicontrol(handles.edit_mu) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_numscancircles, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 2
    hError = errordlg('Number of scan circles must be a positive whole number (e.g. 10-20)', 'Input error', 'modal') ; 
    uicontrol(handles.edit_numscancircles) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_degfromnorth, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) > 360.0 || str2double(sValue) < -360.0
    hError = errordlg('Degrees from North must be a number <= 360 (use negative for CCW)', 'Input error', 'modal') ; 
    uicontrol(handles.edit_degfromnorth) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_lambda, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 0.0 || str2double(sValue) > 1.0
    hError = errordlg('Lamda factor must be a number, >= 0, <= 1.0', 'Input error', 'modal') ; 
    uicontrol(handles.edit_lambda) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_fixedaperture, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 0.0
    hError = errordlg('Fixed aperture must be a number, >= 0', 'Input error', 'modal') ; 
    uicontrol(handles.edit_fixedaperture) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_aperturefactor, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 0.0 || str2double(sValue) > 0.1
    hError = errordlg('Aperture scaling factor must be a number, >= 0, <= 0.1', 'Input error', 'modal') ; 
    uicontrol(handles.edit_aperturefactor) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_apertureexponent, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 0.0 || str2double(sValue) > 1.0 || str2double(sValue) < 0.5 
    hError = errordlg('Aperture scaling exponent must be a number, >= 0.5, <= 1.0', 'Input error', 'modal') ; 
    uicontrol(handles.edit_apertureexponent) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_numblocksmap, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) <= 1.0 
    hError = errordlg('Number of blocks in heat maps must be an integer, >= 2', 'Input error', 'modal') ; 
    uicontrol(handles.edit_numblocksmap) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_numscancirclesgraph, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) <= 1.0 
    hError = errordlg('Number of scan circles along line must be an integer, >= 2', 'Input error', 'modal') ; 
    uicontrol(handles.edit_numscancirclesgraph) ; 
    flagError = true ; 
    return ; 
end ; 

flag_tracelengthgraphalong = get(handles.checkbox_tracelengthgraphalong, 'Value') ;
flag_segmentlengthgraphalong = get(handles.checkbox_segmentlengthgraphalong, 'Value') ;
flag_segmentanglegraphalong = get(handles.checkbox_segmentstrikegraphalong, 'Value') ;
flag_intensitygraphalong = get(handles.checkbox_intensitygraphalong, 'Value') ;
flag_densitygraphalong = get(handles.checkbox_densitygraphalong, 'Value') ;
flag_tracelengthgraphfrom = get(handles.checkbox_tracelengthgraphfrom, 'Value') ;
flag_segmentlengthgraphfrom = get(handles.checkbox_segmentlengthgraphfrom, 'Value') ;
flag_segmentanglegraphfrom = get(handles.checkbox_segmentanglegraphfrom, 'Value') ;

flag_graph_any = false ; 

if sum([flag_tracelengthgraphalong, flag_segmentlengthgraphalong, flag_segmentanglegraphalong, ...
        flag_intensitygraphalong, flag_densitygraphalong, ... 
        flag_tracelengthgraphfrom, flag_segmentlengthgraphfrom, flag_segmentanglegraphfrom]) > 0 
    
    flag_graph_any = true ; 

    linex = handles.graphx ; 
    liney = handles.graphy ; 
    if sum(linex) == 0 && sum(liney) == 0
        hError = errordlg('Please select a line with two end points', 'Input error', 'modal') ; 
        uicontrol(handles.pushbutton_selectgraphendpoints) ; 
        flagError = true ; 
        return ; 
    end 

end ; 

if ~flagError 
    
    traces = getappdata(hObject, 'Traces') ; 

    %   get the current RGB colour string 
    if ~isempty(handles.selmultifile)
        sColour = '[ 0 0 1 ]' ; 
    else
        sColour = char(handles.cstrColours(handles.iCurrentColour)) ; 
    end ; 
    
    if length(traces) > 0 

        %   get values from text boxes 
        if isempty(get(handles.edit_scaling, 'String'))
            nPixels = 0 ; 
        else
            nPixels = str2double(get(handles.edit_scaling, 'String')) ; 
        end ; 
        
        if strcmp(handles.axes_tracemap.YDir, 'reverse')
            flag_reverseY = true ; 
        else
            flag_reverseY = false ; 
        end ;
        
        if strcmp(handles.axes_tracemap.XDir, 'reverse')
            flag_reverseX = true ; 
        else
            flag_reverseX = false ; 
        end ;
        
        nNorth = str2double(get(handles.edit_degfromnorth, 'String')) ;
        nMLE_lc = str2double(get(handles.edit_lc, 'String')) ;  
        nMLE_uc = str2double(get(handles.edit_uc, 'String')) ;  
        
        nSigma1 = str2double(get(handles.edit_sigma1, 'String')) ; 
        nSigma2 = str2double(get(handles.edit_sigma2, 'String')) ; 
        nThetaSigma1 = str2double(get(handles.edit_thetasigma1, 'String')) ; 
        C0 = str2double(get(handles.edit_C0, 'String')) ; 
        Pf = str2double(get(handles.edit_Pf, 'String')) ; 
        mu = str2double(get(handles.edit_mu, 'String')) ; 
        
        nBlocks = str2double(get(handles.edit_numblocksmap, 'String')) ; 
        nPixelsItoY = str2double(get(handles.edit_numpixelsItoY, 'String')) ; 
        nScanCircles = str2double(get(handles.edit_numscancirclesgraph, 'String')) ; 

        %   get values from drop down lists 
        list = get(handles.popupmenu_histolengthbins, 'String') ; 
        nHistoLengthBins = str2double(list{get(handles.popupmenu_histolengthbins, 'Value')}) ;
        
        list = get(handles.popupmenu_histoanglebins, 'String') ; 
        nHistoAngleBins = str2double(list{get(handles.popupmenu_histoanglebins, 'Value')}) ;
        
        list = get(handles.popupmenu_rosebins, 'String') ; 
        nRoseBins = str2double(list{get(handles.popupmenu_rosebins, 'Value')}) ; 
    
        %   call the various maps & graphs 
        flag_tracemap = get(handles.checkbox_tracemap, 'Value') ; 
        flag_shownodes = get(handles.checkbox_shownodes, 'Value') ; 
        flag_sliptendency = get(handles.checkbox_sliptendency, 'Value') ; 
        flag_dilationtendency = get(handles.checkbox_dilationtendency, 'Value') ; 
        flag_fracsuscep = get(handles.checkbox_fracsuscep, 'Value') ; 
        flag_CSF = get(handles.checkbox_CSF, 'Value') ; 
        flag_tracesbylength = get(handles.checkbox_tracemaplength, 'Value') ; 
        flag_segmentsbylength = get(handles.checkbox_segmentmaplength, 'Value') ; 
        flag_segmentsbystrike = get(handles.checkbox_segmentmapstrike, 'Value') ; 
        
        flag_tracemap_any = false ; 
        if sum([flag_tracemap, flag_sliptendency, flag_dilationtendency, ... 
                flag_tracesbylength, flag_segmentsbylength, flag_segmentsbystrike, ...
                flag_fracsuscep, flag_CSF]) > 0 
            flag_tracemap_any = true ; 
        end ; 
        
        flag_intensitymap = get(handles.checkbox_intensitymap, 'Value') ; 
        flag_densitymap = get(handles.checkbox_densitymap, 'Value') ; 
        flag_showcircles = get(handles.checkbox_showcircles, 'Value') ; 
        nCircles = round(str2double(get(handles.edit_numscancircles, 'String'))) ; 
        if  get(handles.checkbox_rotategrid, 'Value')
            gridrotation = round(str2double(get(handles.edit_gridrotation, 'String'))) ; 
        else
            gridrotation = 0;
        end;
         
        flag_histolength = get(handles.checkbox_histolength, 'Value') ; 
        flag_logloglength = get(handles.checkbox_logloglength, 'Value') ;
        flag_crossplot = get(handles.checkbox_crossplot, 'Value') ;
        flag_censor = get(handles.checkbox_censor, 'Value') ; 
        flag_mle = get(handles.checkbox_mle, 'Value') ;
        flag_blocksize = get(handles.checkbox_blocksize, 'Value') ; 
        flag_seglenvario = get(handles.checkbox_seglenvariogram, 'Value') ; 

        flag_histoangle = get(handles.checkbox_histoangle, 'Value') ;
        flag_roseangle = get(handles.checkbox_rose, 'Value') ;
        flag_rosemean = get(handles.checkbox_showrosemean, 'Value') ; 
        flag_roselengthweighted = get(handles.checkbox_roselengthweighted, 'Value') ; 
        flag_rosecolour = get(handles.checkbox_rosecolour, 'Value') ; 
        
        flag_cracktensor = get(handles.checkbox_cracktensor, 'Value') ;

        flag_triangle = get(handles.checkbox_triangle, 'Value') ;
        flag_permellipse = get(handles.checkbox_permellipse, 'Value') ;
        
        if length(handles.cstrColours) > 1 
            flag_multicolour = 1 ; 
        else 
            flag_multicolour = 0 ; 
        end ;

        set(handles.text_message, 'String', 'Running...') ;   
        if flag_tracemap_any 
            guiFracPaQ2Dtracemap(traces, nPixels, nNorth, ...
                    flag_shownodes, flag_reverseY, flag_reverseX, sColour, ...
                    flag_multicolour, flag_tracemap, flag_sliptendency, flag_dilationtendency, ...
                    flag_tracesbylength, flag_segmentsbylength, flag_segmentsbystrike, ...
                    nSigma1, nSigma2, nThetaSigma1, flag_fracsuscep, flag_CSF, C0, Pf, mu) ; 
        end ; 

        flag_length = sum([flag_histolength, flag_logloglength, flag_blocksize, flag_seglenvario]) ; 
        if flag_length
            guiFracPaQ2Dlength(traces, nPixels, nNorth, ...
                               nHistoLengthBins, flag_histolength, flag_logloglength, ...
                               flag_crossplot, flag_mle, flag_censor, ...
                               flag_blocksize, flag_reverseY, flag_reverseX, ...
                               sColour, nMLE_lc, nMLE_uc, flag_seglenvario) ; 
        end ; 

        flag_angle = sum([flag_histoangle, flag_roseangle, flag_cracktensor]) ; 
        if flag_angle
            guiFracPaQ2Dangle(traces, nNorth, nHistoAngleBins, nRoseBins, ...
                                    flag_histoangle, flag_roseangle, ...
                                    flag_reverseY, flag_reverseX, ...
                                    flag_cracktensor, ...
                                    flag_roselengthweighted, flag_rosemean, ...
                                    flag_rosecolour, sColour) ; 
        end ; 

        flag_pattern = sum([flag_intensitymap, flag_densitymap, flag_triangle]) ; 
        if flag_pattern
            guiFracPaQ2Dpattern(traces, nPixels, nBlocks, ... 
                                flag_intensitymap, flag_densitymap, ...
                                flag_triangle, flag_showcircles, ...
                                nCircles, flag_reverseY, flag_reverseX, ...
                                sColour, nPixelsItoY, gridrotation) ; 
        end ; 

        if flag_permellipse
            nLambda = str2double(get(handles.edit_lambda, 'String')) ;
            if get(handles.radiobutton_fixedaperture, 'Value') 
                nFixedAperture = str2double(get(handles.edit_fixedaperture, 'String')) ;
                nApertureFactor = 0 ;
                nApertureExponent = 0 ;
            else 
                nFixedAperture = 0 ; 
                nApertureFactor = str2double(get(handles.edit_aperturefactor, 'String')) ;
                nApertureExponent = str2double(get(handles.edit_apertureexponent, 'String')) ;
            end ; 
            guiFracPaQ2Dtensor(traces, nNorth, nLambda, ...
                                    flag_reverseY, flag_reverseX, ...
                                        nFixedAperture, nApertureFactor, nApertureExponent, ...
                                            nPixels, sColour) ; 
        end ; 

        if flag_graph_any 
            guiFracPaQ2Dgraphs(traces, nPixels, ...
                                handles.graphx, handles.graphy, sColour, ...
                                nScanCircles, ... 
                                flag_reverseX, flag_reverseY, ... 
                                flag_tracelengthgraphalong, ...
                                flag_segmentlengthgraphalong, ...
                                flag_segmentanglegraphalong, ...
                                flag_intensitygraphalong, ...
                                flag_densitygraphalong, ...
                                flag_tracelengthgraphfrom, ...
                                flag_segmentlengthgraphfrom, ...
                                flag_segmentanglegraphfrom) ; 
        end 
        
        set(handles.text_message, 'String', 'All done. Check current folder for plot figure .tif files.') ;   
        uicontrol(handles.pushbutton_run) ; 

    else 
        
        errordlg('No traces found; check the input file is OK.') ; 
        uicontrol(handles.pushbutton_browse) ; 

    end ; 
   
end ; 

%   only validate wavelet input if wavelet analysis selected...
if get(handles.checkbox_wavelet, 'Value') 

    flagWaveError = false ; 

    sValue = get(handles.edit_aString, 'String') ; 
    if isempty(str2double(sValue)) 
        hError = errordlg('a values must all be whole numbers, e.g. 2 4 8', 'Input error', 'modal') ; 
        uicontrol(handles.edit_aString) ; 
        flagWaveError = true ; 
        return ; 
    else 
        handles.a = str2num(sValue) ; 
    end ; 

    sValue = get(handles.edit_LString, 'String') ; 
    if isempty(str2double(sValue)) 
        hError = errordlg('L values must all be whole numbers, e.g. 2 4 8', 'Input error', 'modal') ; 
        uicontrol(handles.edit_LString) ; 
        flagWaveError = true ; 
        return ; 
    else 
        handles.L = 1 ./ str2num(sValue) ; 
    end ; 

    sValue = get(handles.edit_DegString, 'String') ; 
    if isnan(str2double(sValue)) || str2double(sValue) < 1 || str2double(sValue) > 360
        hError = errordlg('Angular increment must be a whole number, 1-360', 'Input error', 'modal') ; 
        uicontrol(handles.edit_DegString) ; 
        flagWaveError = true ; 
        return ; 
    else 
        handles.nTheta = str2num(sValue) ; 
    end ; 

    if get(handles.radiobutton_Morlet, 'Value')
        handles.fMorlet = true ; 
    else 
        handles.fMorlet = false ; 
    end ; 

    set(handles.text_message, 'String', 'Running...') ;   
    if ~flagWaveError
        guiFracPaQ2Dwavelet(traces, handles.a, handles.L, ... 
                            handles.nTheta, handles.fMorlet, '[0 0 1]', ...
                            flag_reverseX, flag_reverseY) ; 
    end ; 
    set(handles.text_message, 'String', 'All done. Check current folder for plot figure .tif files.') ;   
    uicontrol(handles.pushbutton_run) ; 

end ; 

% --- Executes on button press in radiobutton_image.
function radiobutton_image_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_image
set(handles.edit_houghpeaks, 'Enable', 'on') ; 
set(handles.edit_houghthreshold, 'Enable', 'on') ; 
set(handles.edit_fillgap, 'Enable', 'on') ; 
set(handles.edit_minlength, 'Enable', 'on') ; 


% --- Executes on button press in radiobutton_node.
function radiobutton_node_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_node (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_node
set(handles.edit_houghpeaks, 'Enable', 'off') ; 
set(handles.edit_houghthreshold, 'Enable', 'off') ; 
set(handles.edit_fillgap, 'Enable', 'off') ; 
set(handles.edit_minlength, 'Enable', 'off') ; 


% --- Executes during object creation, after setting all properties.
function figure_main_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in popupmenu_histolengthbins.
function popupmenu_histolengthbins_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_histolengthbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_histolengthbins contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_histolengthbins


% --- Executes during object creation, after setting all properties.
function popupmenu_histolengthbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_histolengthbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_histoanglebins.
function popupmenu_histoanglebins_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_histoanglebins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_histoanglebins contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_histoanglebins


% --- Executes during object creation, after setting all properties.
function popupmenu_histoanglebins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_histoanglebins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_rosebins.
function popupmenu_rosebins_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_rosebins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_rosebins contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_rosebins


% --- Executes during object creation, after setting all properties.
function popupmenu_rosebins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_rosebins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_flipy.
function pushbutton_flipy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_flipy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.axes_tracemap.YDir, 'reverse')
    handles.axes_tracemap.YDir = 'normal' ;
else
    handles.axes_tracemap.YDir = 'reverse' ;
end ;     
set(handles.text_message, 'String', 'Y-axis flipped.') ;   

% --- Executes on button press in pushbutton_flipx.
function pushbutton_flipx_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_flipx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.axes_tracemap.XDir, 'reverse')
    handles.axes_tracemap.XDir = 'normal' ;
else
    handles.axes_tracemap.XDir = 'reverse' ;
end ;     
set(handles.text_message, 'String', 'X-axis flipped.') ;   


% --- Executes on button press in checkbox_cracktensor.
function checkbox_cracktensor_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cracktensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_cracktensor


% --- Executes on button press in checkbox_roselengthweighted.
function checkbox_roselengthweighted_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_roselengthweighted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_roselengthweighted



function edit_apertureexponent_Callback(hObject, eventdata, handles)
% hObject    handle to edit_apertureexponent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_apertureexponent as text
%        str2double(get(hObject,'String')) returns contents of edit_apertureexponent as a double


% --- Executes during object creation, after setting all properties.
function edit_apertureexponent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_apertureexponent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_aperturefactor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_aperturefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_aperturefactor as text
%        str2double(get(hObject,'String')) returns contents of edit_aperturefactor as a double


% --- Executes during object creation, after setting all properties.
function edit_aperturefactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_aperturefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fixedaperture_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fixedaperture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fixedaperture as text
%        str2double(get(hObject,'String')) returns contents of edit_fixedaperture as a double


% --- Executes during object creation, after setting all properties.
function edit_fixedaperture_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fixedaperture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_showrosemean.
function checkbox_showrosemean_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showrosemean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showrosemean


% --- Executes on button press in checkbox_blocksize.
function checkbox_blocksize_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_blocksize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_blocksize



function edit_lc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lc as text
%        str2double(get(hObject,'String')) returns contents of edit_lc as a double


% --- Executes during object creation, after setting all properties.
function edit_lc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_uc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_uc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_uc as text
%        str2double(get(hObject,'String')) returns contents of edit_uc as a double


% --- Executes during object creation, after setting all properties.
function edit_uc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_uc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_rosecolour.
function checkbox_rosecolour_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_rosecolour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_rosecolour


% --- Executes on button press in checkbox_seglenvariogram.
function checkbox_seglenvariogram_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_seglenvariogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_seglenvariogram


% --- Executes on button press in pbMapsTab.
function pbMapsTab_Callback(hObject, eventdata, handles)
% hObject    handle to pbMapsTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uibgTabMaps, 'Visible', 'on') ; 
set(handles.uibgTabLengths, 'Visible', 'off') ; 
set(handles.uibgTabAngles, 'Visible', 'off') ; 
set(handles.uibgTabFluid, 'Visible', 'off') ; 
set(handles.uibgTabGraphs, 'Visible', 'off') ; 
set(handles.uibgTabWavelets, 'Visible', 'off') ; 


% --- Executes on button press in pbLengthsTab.
function pbLengthsTab_Callback(hObject, eventdata, handles)
% hObject    handle to pbLengthsTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uibgTabMaps, 'Visible', 'off') ; 
set(handles.uibgTabLengths, 'Visible', 'on') ; 
set(handles.uibgTabAngles, 'Visible', 'off') ; 
set(handles.uibgTabFluid, 'Visible', 'off') ; 
set(handles.uibgTabGraphs, 'Visible', 'off') ; 
set(handles.uibgTabWavelets, 'Visible', 'off') ; 


% --- Executes on button press in pbAnglesTab.
function pbAnglesTab_Callback(hObject, eventdata, handles)
% hObject    handle to pbAnglesTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uibgTabMaps, 'Visible', 'off') ; 
set(handles.uibgTabLengths, 'Visible', 'off') ; 
set(handles.uibgTabAngles, 'Visible', 'on') ; 
set(handles.uibgTabFluid, 'Visible', 'off') ; 
set(handles.uibgTabGraphs, 'Visible', 'off') ; 
set(handles.uibgTabWavelets, 'Visible', 'off') ; 


% --- Executes on button press in pbFluidTab.
function pbFluidTab_Callback(hObject, eventdata, handles)
% hObject    handle to pbFluidTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uibgTabMaps, 'Visible', 'off') ; 
set(handles.uibgTabLengths, 'Visible', 'off') ; 
set(handles.uibgTabAngles, 'Visible', 'off') ; 
set(handles.uibgTabFluid, 'Visible', 'on') ; 
set(handles.uibgTabGraphs, 'Visible', 'off') ; 
set(handles.uibgTabWavelets, 'Visible', 'off') ; 


% --- Executes on button press in pbWaveletsTab.
function pbWaveletsTab_Callback(hObject, eventdata, handles)
% hObject    handle to pbWaveletsTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uibgTabMaps, 'Visible', 'off') ; 
set(handles.uibgTabLengths, 'Visible', 'off') ; 
set(handles.uibgTabAngles, 'Visible', 'off') ; 
set(handles.uibgTabFluid, 'Visible', 'off') ; 
set(handles.uibgTabGraphs, 'Visible', 'off') ; 
set(handles.uibgTabWavelets, 'Visible', 'on') ; 


% --- Executes on button press in checkbox_tracemaplength.
function checkbox_tracemaplength_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_tracemaplength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_tracemaplength


% --- Executes on button press in checkbox_segmentmaplength.
function checkbox_segmentmaplength_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_segmentmaplength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_segmentmaplength


% --- Executes on button press in checkbox_segmentmapstrike.
function checkbox_segmentmapstrike_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_segmentmapstrike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_segmentmapstrike



function edit_sigma1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sigma1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sigma1 as text
%        str2double(get(hObject,'String')) returns contents of edit_sigma1 as a double


% --- Executes during object creation, after setting all properties.
function edit_sigma1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sigma1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_dilationtendency.
function checkbox_dilationtendency_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_dilationtendency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_dilationtendency
if get(hObject, 'Value') || ...
   get(handles.checkbox_sliptendency, 'Value') || ... 
   get(handles.checkbox_fracsuscep, 'Value') || ... 
   get(handles.checkbox_CSF, 'Value')
    set(handles.edit_sigma1, 'Enable', 'on') ; 
    set(handles.edit_sigma2, 'Enable', 'on') ; 
    set(handles.edit_thetasigma1, 'Enable', 'on') ; 
    set(handles.edit_C0, 'Enable', 'on') ; 
    set(handles.edit_Pf, 'Enable', 'on') ; 
    set(handles.edit_mu, 'Enable', 'on') ; 
    set(handles.label_sigma1, 'Enable', 'on') ; 
    set(handles.label_sigma2, 'Enable', 'on') ; 
    set(handles.label_thetasigma1, 'Enable', 'on') ; 
    set(handles.label_C0, 'Enable', 'on') ; 
    set(handles.label_Pf, 'Enable', 'on') ; 
    set(handles.label_mu, 'Enable', 'on') ; 
else 
    set(handles.edit_sigma1, 'Enable', 'off') ; 
    set(handles.edit_sigma2, 'Enable', 'off') ; 
    set(handles.edit_thetasigma1, 'Enable', 'off') ; 
    set(handles.edit_C0, 'Enable', 'off') ; 
    set(handles.edit_Pf, 'Enable', 'off') ; 
    set(handles.edit_mu, 'Enable', 'off') ; 
    set(handles.label_sigma1, 'Enable', 'off') ; 
    set(handles.label_sigma2, 'Enable', 'off') ; 
    set(handles.label_thetasigma1, 'Enable', 'off') ; 
    set(handles.label_C0, 'Enable', 'off') ; 
    set(handles.label_Pf, 'Enable', 'off') ; 
    set(handles.label_mu, 'Enable', 'off') ; 
end ; 


% --- Executes on button press in checkbox_sliptendency.
function checkbox_sliptendency_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_sliptendency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_sliptendency
if get(hObject, 'Value') || ...
   get(handles.checkbox_dilationtendency, 'Value') || ... 
   get(handles.checkbox_fracsuscep, 'Value') || ... 
   get(handles.checkbox_CSF, 'Value')
    set(handles.edit_sigma1, 'Enable', 'on') ; 
    set(handles.edit_sigma2, 'Enable', 'on') ; 
    set(handles.edit_thetasigma1, 'Enable', 'on') ; 
    set(handles.edit_C0, 'Enable', 'on') ; 
    set(handles.edit_Pf, 'Enable', 'on') ; 
    set(handles.edit_mu, 'Enable', 'on') ; 
    set(handles.label_sigma1, 'Enable', 'on') ; 
    set(handles.label_sigma2, 'Enable', 'on') ; 
    set(handles.label_thetasigma1, 'Enable', 'on') ; 
    set(handles.label_C0, 'Enable', 'on') ; 
    set(handles.label_Pf, 'Enable', 'on') ; 
    set(handles.label_mu, 'Enable', 'on') ; 
else 
    set(handles.edit_sigma1, 'Enable', 'off') ; 
    set(handles.edit_sigma2, 'Enable', 'off') ; 
    set(handles.edit_thetasigma1, 'Enable', 'off') ; 
    set(handles.edit_C0, 'Enable', 'off') ; 
    set(handles.edit_Pf, 'Enable', 'off') ; 
    set(handles.edit_mu, 'Enable', 'off') ; 
    set(handles.label_sigma1, 'Enable', 'off') ; 
    set(handles.label_sigma2, 'Enable', 'off') ; 
    set(handles.label_thetasigma1, 'Enable', 'off') ; 
    set(handles.label_C0, 'Enable', 'off') ; 
    set(handles.label_Pf, 'Enable', 'off') ; 
    set(handles.label_mu, 'Enable', 'off') ; 
end ; 


function edit_sigma2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sigma2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sigma2 as text
%        str2double(get(hObject,'String')) returns contents of edit_sigma2 as a double


% --- Executes during object creation, after setting all properties.
function edit_sigma2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sigma2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_thetasigma1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thetasigma1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thetasigma1 as text
%        str2double(get(hObject,'String')) returns contents of edit_thetasigma1 as a double


% --- Executes during object creation, after setting all properties.
function edit_thetasigma1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thetasigma1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_aString_Callback(hObject, eventdata, handles)
% hObject    handle to edit_aString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_aString as text
%        str2double(get(hObject,'String')) returns contents of edit_aString as a double


% --- Executes during object creation, after setting all properties.
function edit_aString_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_aString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_LString_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LString as text
%        str2double(get(hObject,'String')) returns contents of edit_LString as a double


% --- Executes during object creation, after setting all properties.
function edit_LString_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_DegString_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DegString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DegString as text
%        str2double(get(hObject,'String')) returns contents of edit_DegString as a double


% --- Executes during object creation, after setting all properties.
function edit_DegString_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DegString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_Morlet.
function radiobutton_Morlet_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Morlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_Morlet
if get(hObject, 'Value') 
    handles.fMorlet = true ; 
end ; 

% --- Executes on button press in radiobutton_Mexican.
function radiobutton_Mexican_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Mexican (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_Mexican
if get(hObject, 'Value') 
    handles.fMorlet = false ; 
end ; 

% --- Executes on button press in checkbox_wavelet.
function checkbox_wavelet_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_wavelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_wavelet
if get(hObject, 'Value') 
    set(handles.edit_aString, 'Enable', 'on') ; 
    set(handles.edit_LString, 'Enable', 'on') ; 
    set(handles.edit_DegString, 'Enable', 'on') ; 
    set(handles.radiobutton_Morlet, 'Enable', 'on') ; 
    set(handles.radiobutton_Mexican, 'Enable', 'on') ; 
    set(handles.label_aString, 'Enable', 'on') ; 
    set(handles.label_LString, 'Enable', 'on') ; 
    set(handles.label_DegString, 'Enable', 'on') ; 
else 
    set(handles.edit_aString, 'Enable', 'off') ; 
    set(handles.edit_LString, 'Enable', 'off') ; 
    set(handles.edit_DegString, 'Enable', 'off') ; 
    set(handles.radiobutton_Morlet, 'Enable', 'off') ; 
    set(handles.radiobutton_Mexican, 'Enable', 'off') ; 
    set(handles.label_aString, 'Enable', 'off') ; 
    set(handles.label_LString, 'Enable', 'off') ; 
    set(handles.label_DegString, 'Enable', 'off') ; 
end ; 


% --- Executes on button press in pbGraphsTab.
function pbGraphsTab_Callback(hObject, eventdata, handles)
% hObject    handle to pbGraphsTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uibgTabMaps, 'Visible', 'off') ; 
set(handles.uibgTabLengths, 'Visible', 'off') ; 
set(handles.uibgTabAngles, 'Visible', 'off') ; 
set(handles.uibgTabFluid, 'Visible', 'off') ; 
set(handles.uibgTabGraphs, 'Visible', 'on') ; 
set(handles.uibgTabWavelets, 'Visible', 'off') ; 


% --- Executes on button press in checkbox_segmentlengthgraphalong.
function checkbox_segmentlengthgraphalong_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_segmentlengthgraphalong (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_segmentlengthgraphalong


% --- Executes on button press in checkbox_segmentstrikegraphalong.
function checkbox_segmentstrikegraphalong_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_segmentstrikegraphalong (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_segmentstrikegraphalong


% --- Executes on button press in checkbox_intensitygraphalong.
function checkbox_intensitygraphalong_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_intensitygraphalong (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_intensitygraphalong
if get(hObject, 'Value') || get(handles.checkbox_densitygraphalong, 'Value')
    set(handles.label_numscancirclesgraph, 'Enable', 'on') ; 
    set(handles.edit_numscancirclesgraph, 'Enable', 'on') ; 
else 
    set(handles.label_numscancirclesgraph, 'Enable', 'off') ; 
    set(handles.edit_numscancirclesgraph, 'Enable', 'off') ; 
end ; 

% --- Executes on button press in checkbox_densitygraphalong.
function checkbox_densitygraphalong_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_densitygraphalong (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_densitygraphalong
if get(hObject, 'Value') || get(handles.checkbox_intensitygraphalong, 'Value')
    set(handles.label_numscancirclesgraph, 'Enable', 'on') ; 
    set(handles.edit_numscancirclesgraph, 'Enable', 'on') ; 
else 
    set(handles.label_numscancirclesgraph, 'Enable', 'off') ; 
    set(handles.edit_numscancirclesgraph, 'Enable', 'off') ; 
end ; 

% --- Executes on button press in checkbox_tracelengthgraphalong.
function checkbox_tracelengthgraphalong_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_tracelengthgraphalong (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_tracelengthgraphalong


% --- Executes on button press in checkbox_tracelengthgraphfrom.
function checkbox_tracelengthgraphfrom_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_tracelengthgraphfrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_tracelengthgraphfrom


% --- Executes on button press in checkbox_segmentlengthgraphfrom.
function checkbox_segmentlengthgraphfrom_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_segmentlengthgraphfrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_segmentlengthgraphfrom


% --- Executes on button press in checkbox_segmentanglegraphfrom.
function checkbox_segmentanglegraphfrom_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_segmentanglegraphfrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_segmentanglegraphfrom



function edit_numblocksmap_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numblocksmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numblocksmap as text
%        str2double(get(hObject,'String')) returns contents of edit_numblocksmap as a double


% --- Executes during object creation, after setting all properties.
function edit_numblocksmap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numblocksmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_numscancirclesgraph_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numscancirclesgraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numscancirclesgraph as text
%        str2double(get(hObject,'String')) returns contents of edit_numscancirclesgraph as a double


% --- Executes during object creation, after setting all properties.
function edit_numscancirclesgraph_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numscancirclesgraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_FilenameTag_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FilenameTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FilenameTag as text
%        str2double(get(hObject,'String')) returns contents of edit_FilenameTag as a double
global sTag ; 

sTag = ['_', get(hObject,'String')] ; 

% --- Executes during object creation, after setting all properties.
function edit_FilenameTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FilenameTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_numpixelsItoY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numpixelsItoY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numpixelsItoY as text
%        str2double(get(hObject,'String')) returns contents of edit_numpixelsItoY as a double


% --- Executes during object creation, after setting all properties.
function edit_numpixelsItoY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numpixelsItoY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_fracsuscep.
function checkbox_fracsuscep_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_fracsuscep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_fracsuscep
if get(hObject, 'Value') || ...
   get(handles.checkbox_dilationtendency, 'Value') || ...
   get(handles.checkbox_sliptendency, 'Value') || ... 
   get(handles.checkbox_CSF, 'Value')
    set(handles.edit_sigma1, 'Enable', 'on') ; 
    set(handles.edit_sigma2, 'Enable', 'on') ; 
    set(handles.edit_thetasigma1, 'Enable', 'on') ; 
    set(handles.edit_C0, 'Enable', 'on') ; 
    set(handles.edit_Pf, 'Enable', 'on') ; 
    set(handles.edit_mu, 'Enable', 'on') ; 
    set(handles.label_sigma1, 'Enable', 'on') ; 
    set(handles.label_sigma2, 'Enable', 'on') ; 
    set(handles.label_thetasigma1, 'Enable', 'on') ; 
    set(handles.label_C0, 'Enable', 'on') ; 
    set(handles.label_Pf, 'Enable', 'on') ; 
    set(handles.label_mu, 'Enable', 'on') ; 
else 
    set(handles.edit_sigma1, 'Enable', 'off') ; 
    set(handles.edit_sigma2, 'Enable', 'off') ; 
    set(handles.edit_thetasigma1, 'Enable', 'off') ; 
    set(handles.edit_C0, 'Enable', 'off') ; 
    set(handles.edit_Pf, 'Enable', 'off') ; 
    set(handles.edit_mu, 'Enable', 'off') ; 
    set(handles.label_sigma1, 'Enable', 'off') ; 
    set(handles.label_sigma2, 'Enable', 'off') ; 
    set(handles.label_thetasigma1, 'Enable', 'off') ; 
    set(handles.label_C0, 'Enable', 'off') ; 
    set(handles.label_Pf, 'Enable', 'off') ; 
    set(handles.label_mu, 'Enable', 'off') ; 
end ; 


% --- Executes on button press in checkbox_CSF.
function checkbox_CSF_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_CSF
if get(hObject, 'Value') || ...
   get(handles.checkbox_dilationtendency, 'Value') || ...
   get(handles.checkbox_sliptendency, 'Value') || ...
   get(handles.checkbox_fracsuscep, 'Value')
    set(handles.edit_sigma1, 'Enable', 'on') ; 
    set(handles.edit_sigma2, 'Enable', 'on') ; 
    set(handles.edit_thetasigma1, 'Enable', 'on') ; 
    set(handles.edit_C0, 'Enable', 'on') ; 
    set(handles.edit_Pf, 'Enable', 'on') ; 
    set(handles.edit_mu, 'Enable', 'on') ; 
    set(handles.label_sigma1, 'Enable', 'on') ; 
    set(handles.label_sigma2, 'Enable', 'on') ; 
    set(handles.label_thetasigma1, 'Enable', 'on') ; 
    set(handles.label_C0, 'Enable', 'on') ; 
    set(handles.label_Pf, 'Enable', 'on') ; 
    set(handles.label_mu, 'Enable', 'on') ; 
else 
    set(handles.edit_sigma1, 'Enable', 'off') ; 
    set(handles.edit_sigma2, 'Enable', 'off') ; 
    set(handles.edit_thetasigma1, 'Enable', 'off') ; 
    set(handles.edit_C0, 'Enable', 'off') ; 
    set(handles.edit_Pf, 'Enable', 'off') ; 
    set(handles.edit_mu, 'Enable', 'off') ; 
    set(handles.label_sigma1, 'Enable', 'off') ; 
    set(handles.label_sigma2, 'Enable', 'off') ; 
    set(handles.label_thetasigma1, 'Enable', 'off') ; 
    set(handles.label_C0, 'Enable', 'off') ; 
    set(handles.label_Pf, 'Enable', 'off') ; 
    set(handles.label_mu, 'Enable', 'off') ; 
end ; 


function edit_C0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_C0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_C0 as text
%        str2double(get(hObject,'String')) returns contents of edit_C0 as a double


% --- Executes during object creation, after setting all properties.
function edit_C0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_C0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Pf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Pf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Pf as text
%        str2double(get(hObject,'String')) returns contents of edit_Pf as a double


% --- Executes during object creation, after setting all properties.
function edit_Pf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Pf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_mu_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mu as text
%        str2double(get(hObject,'String')) returns contents of edit_mu as a double


% --- Executes during object creation, after setting all properties.
function edit_mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_rotategrid.
function checkbox_rotategrid_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_rotategrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_rotategrid
if get(hObject,'Value') 
    set(handles.edit_gridrotation, 'Enable', 'on') ; 
    set(handles.label_gridrotation, 'Enable', 'on') ; 
else
    set(handles.edit_gridrotation, 'Enable', 'off') ; 
    set(handles.label_gridrotation, 'Enable', 'off') ;     
end;



function edit_gridrotation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gridrotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gridrotation as text
%        str2double(get(hObject,'String')) returns contents of edit_gridrotation as a double


% --- Executes during object creation, after setting all properties.
function edit_gridrotation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gridrotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_circlespacing_Callback(hObject, eventdata, handles)
% hObject    handle to edit_circlespacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_circlespacing as text
%        str2double(get(hObject,'String')) returns contents of edit_circlespacing as a double
   
    if length(handles.traces) > 0
        spacing = str2double(get(hObject,'String'));
        if get(handles.checkbox_rotategrid, 'Value')
            rotation = round(str2double(get(handles.edit_gridrotation, 'String'))) ; 
        else
            rotation = 0;
        end
        [dist, iCount, jCount] = circleNumberFromSpacing(handles.traces, rotation, spacing);
        disp(min(iCount, jCount));
        set(handles.edit_numscancircles,'String', num2str(min(iCount, jCount)));
    end


% --- Executes during object creation, after setting all properties.
function edit_circlespacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_circlespacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_saveconfig.
function pushbutton_saveconfig_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_saveconfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sDefaultName = [get(handles.edit_FilenameTag,'String'),'.fpq'];
[sFile, sPath] = uiputfile(sDefaultName);

if ~isequal(sFile,0) && ~isequal(sPath,0)
    hMain = get(hObject,'parent');
    hAll = findall(hMain, 'Type', 'UIControl');
    S.Control=struct('Tag',{},'Style',{},'Enable',{},'Visible',{},'Value',{}, 'String',{});
    
    if isfield(handles, 'selfile')
        S.File=struct('Name',char(handles.selfile),'Path',char(handles.selpath));
    end
    for i = 1:length(hAll)
        S.Control(i).Tag = get(hAll(i),'Tag');
        S.Control(i).Style = get(hAll(i),'Style');
        S.Control(i).Enable = get(hAll(i),'Enable');
        S.Control(i).Visible = get(hAll(i),'Visible');
        if get(hAll(i),'Style') ~= "pushbutton" && get(hAll(i),'Style') ~= "text" && get(hAll(i),'Style') ~= "popupmenu"            
            if get(hAll(i),'Style') == "edit"
                S.Control(i).String = get(hAll(i),'String');
            else
                S.Control(i).Value = get(hAll(i),'Value');
            end        
        end
    end

    writestruct(S,[sPath,sFile],"FileType","xml");
end


% --- Executes on button press in pushbutton_loadconfig.
function pushbutton_loadconfig_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadconfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sFileTypes = {'*.fpq'} ; 

[sFile, sPath] = uigetfile(sFileTypes);

if ~isequal(sFile,0) && ~isequal(sPath,0)
    S = readstruct([sPath,sFile],"FileType","xml");
    drawnow limitrate nocallbacks;
    f = waitbar(0, ['Loading Config ',sFile]);
 
    for i = 1:length(S.Control)
        waitbar(i / length(S.Control), f);
        h = findobj('Tag', S.Control(i).Tag);
        if length(h) > 0
            set(h,'Visible', S.Control(i).Visible);
            set(h,'Enable', S.Control(i).Enable);
            if get(h,'Style') ~= "pushbutton" && get(h,'Style') ~= "text" && get(h,'Style') ~= "popupmenu"            
                if get(h,'Style') == "edit"
                    set(h, 'String', S.Control(i).String);
                else
                    set(h, 'Value', S.Control(i).Value);
                end
            end
        end
    end
    close(f);
    drawnow;

    if isfield(S, 'File')
        handles.selfile = S.File.Name;
        handles.selpath = S.File.Path;
        handles.iCurrentColour = 1;
        handles.cstrColours = cellstr('[0, 0, 1]') ; 
        guidata(handles.pushbutton_browse, handles);
        
        cla(handles.axes_tracemap);
    end
end
