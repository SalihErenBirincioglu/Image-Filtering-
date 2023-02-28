function varargout = assignment(varargin)
% ASSIGNMENT MATLAB code for assignment.fig
%      ASSIGNMENT, by itself, creates a new ASSIGNMENT or raises the existing
%      singleton*.
%
%      H = ASSIGNMENT returns the handle to a new ASSIGNMENT or the handle to
%      the existing singleton*.
%
%      ASSIGNMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ASSIGNMENT.M with the given input arguments.
%
%      ASSIGNMENT('Property','Value',...) creates a new ASSIGNMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before assignment_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to assignment_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help assignment

% Last Modified by GUIDE v2.5 21-May-2022 20:51:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @assignment_OpeningFcn, ...
                   'gui_OutputFcn',  @assignment_OutputFcn, ...
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


% --- Executes just before assignment is made visible.
function assignment_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to assignment (see VARARGIN)
set(handles.slider1,'value',20)
data = get(handles.slider1,'Value');
data1 = num2str(data);
set(handles.edit1,'String',data1);
set(handles.ideal, 'Value', 1);
set(handles.low, 'Value', 1);
% Choose default command line output for assignment
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes assignment wait for user response (see UIRESUME)
% uiwait(handles.figure1);

I = im2double(imread('lena.jpg'));
axes(handles.axes1);
imshow(I);

[m n]=size (I);
H = zeros(m,n);
duv= zeros(m,n);
D0 = get(handles.slider1,'Value');
if (handles.low.Value == 1 & handles.ideal.Value == 1)
    %IDEAL FILTER
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            if (Duv <D0)
            H(i,j) = 1;
            end
        end
    end
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';

end

% --- Outputs from this function are returned to the command line.
function varargout = assignment_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
data = get(handles.slider1,'Value');
data1 = num2str(data);
set(handles.edit1,'String',data1);

I = im2double(imread('lena.jpg'));
[m n]=size (I);
H = zeros(m,n);
D0 = get(handles.slider1,'Value');
if (handles.low.Value == 1 & handles.ideal.Value == 1)
    %IDEAL FILTER
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            if (Duv <D0)
            H(i,j) = 1;
            end
        end
    end
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';

end

if (handles.high.Value == 1 & handles.ideal.Value == 1) %%high pass ideal
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            if (Duv < D0)
            H(i,j) = 1;
            end
        end
    end
    H=1-H;
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end

if (handles.low.Value == 1 & handles.btw.Value == 1)  
    n_order = 1;
    %Butterworth FILTER
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            H(i,j)=1/(1+Duv/D0)^(2*n_order);
        end
    end
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end
if (handles.high.Value == 1 & handles.btw.Value == 1) %%high pass btw
    n_order = 1;
    %Butterworth FILTER
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            H(i,j)=1/(1+Duv/D0)^(2*n_order);
        end
    end
    H=1-H;
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end

if (handles.low.Value == 1 & handles.gaussian.Value == 1)
        %Gaussian FILTER
        for i = 1:m
            for j = 1:n
                Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
                duv(i,j)=Duv;
                H(i,j)=exp(-Duv^2/(2*D0^2));
            end
        end
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)
    
    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end
if (handles.high.Value == 1 & handles.gaussian.Value == 1)
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            H(i,j)=exp(-Duv^2/(2*D0^2));
        end
    end
    H=1-H;
    axes(handles.axes2);
    imshow(H)
    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end



% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



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


% --------------------------------------------------------------------
function uipanel1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in uibuttongroup2.
function uibuttongroup2_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup2 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = im2double(imread('lena.jpg'));
[m n]=size (I);
H = zeros(m,n);
duv= zeros(m,n);
D0 = get(handles.slider1,'Value');
if (handles.low.Value == 1 & handles.ideal.Value == 1)
    %IDEAL FILTER
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            if (Duv <D0)
            H(i,j) = 1;
            end
        end
    end
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';

end

if (handles.high.Value == 1 & handles.ideal.Value == 1) %%high pass ideal
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            if (Duv < D0)
            H(i,j) = 1;
            end
        end
    end
    H=1-H;
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end

if (handles.low.Value == 1 & handles.btw.Value == 1)  
    n_order = 1;
    %Butterworth FILTER
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            H(i,j)=1/(1+Duv/D0)^(2*n_order);
        end
    end
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end
if (handles.high.Value == 1 & handles.btw.Value == 1) %%high pass btw
    n_order = 1;
    %Butterworth FILTER
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            H(i,j)=1/(1+Duv/D0)^(2*n_order);
        end
    end
    H=1-H;
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end

if (handles.low.Value == 1 & handles.gaussian.Value == 1)
        %Gaussian FILTER
        for i = 1:m
            for j = 1:n
                Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
                duv(i,j)=Duv;
                H(i,j)=exp(-Duv^2/(2*D0^2));
            end
        end
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)
    
    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end
if (handles.high.Value == 1 & handles.gaussian.Value == 1)
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            H(i,j)=exp(-Duv^2/(2*D0^2));
        end
    end
    H=1-H;
    axes(handles.axes2);
    imshow(H)
    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end


% --- Executes when selected object is changed in uibuttongroup3.
function uibuttongroup3_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup3 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = im2double(imread('lena.jpg'));
[m n]=size (I);
H = zeros(m,n);
D0 = get(handles.slider1,'Value');
if (handles.low.Value == 1 & handles.ideal.Value == 1)
    %IDEAL FILTER
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            if (Duv <D0)
            H(i,j) = 1;
            end
        end
    end
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';

end

if (handles.high.Value == 1 & handles.ideal.Value == 1) %%high pass ideal
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            if (Duv < D0)
            H(i,j) = 1;
            end
        end
    end
    H=1-H;
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end

if (handles.low.Value == 1 & handles.btw.Value == 1)  
    n_order = 1;
    %Butterworth FILTER
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            H(i,j)=1/(1+Duv/D0)^(2*n_order);
        end
    end
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end
if (handles.high.Value == 1 & handles.btw.Value == 1) %%high pass btw
    n_order = 1;
    %Butterworth FILTER
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            H(i,j)=1/(1+Duv/D0)^(2*n_order);
        end
    end
    H=1-H;
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end

if (handles.low.Value == 1 & handles.gaussian.Value == 1)
        %Gaussian FILTER
        for i = 1:m
            for j = 1:n
                Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
                duv(i,j)=Duv;
                H(i,j)=exp(-Duv^2/(2*D0^2));
            end
        end
    axes(handles.axes2);
    imshow(H)

    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)
    
    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end
if (handles.high.Value == 1 & handles.gaussian.Value == 1)
    for i = 1:m
        for j = 1:n
            Duv = sqrt(power(i-m/2,2)+power(j-n/2,2));
            duv(i,j)=Duv;
            H(i,j)=exp(-Duv^2/(2*D0^2));
        end
    end
    H=1-H;
    axes(handles.axes2);
    imshow(H)
    F = fft2(I);
    Fshift = fftshift(F);
    F1_shift = Fshift.*H;
    F1 = ifftshift(F1_shift);
    G = real(ifft2(F1));
    axes(handles.axes4);
    imshow(G)

    axes(handles.axes3);
    plot(duv,H)
    xlabel 'D(U,V)';
    ylabel 'H(U,V)';
    title 'Graphic of filter';
end
