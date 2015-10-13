function varargout = pEM_gui(varargin)
%-------------------------------------------------------------------------- 
% GUI to employ pEM.  Please see the README.txt for more information
%
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

% load path
addpath('pEM');

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pEM_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @pEM_gui_OutputFcn, ...
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


% --- Executes just before pEM_gui is made visible.
function pEM_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pEM_gui (see VARARGIN)

% Choose default command line output for pEM_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pEM_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pEM_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function dt_Callback(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt as text
%        str2double(get(hObject,'String')) returns contents of dt as a double


% --- Executes during object creation, after setting all properties.
function dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dE_Callback(hObject, eventdata, handles)
% hObject    handle to dE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dE as text
%        str2double(get(hObject,'String')) returns contents of dE as a double


% --- Executes during object creation, after setting all properties.
function dE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% RUN pEM
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = guidata(hObject);
X = data.X;
numTracks = length(X);

% load movie parameters
dt = get(findobj('Tag', 'dt'), 'String');
dE = get(findobj('Tag', 'dE'), 'String');
set(findobj('Tag', 'dt'), 'Value',str2double(dt));
set(findobj('Tag', 'dE'), 'Value',str2double(dE));
dt = get(findobj('Tag', 'dt'), 'Value');
dE = get(findobj('Tag', 'dE'), 'Value');
dim = size(X{1},2);
trackInfo.numberOfTracks = numTracks;        % number of tracks
trackInfo.dimensions = dim;                    % particle track dimensions
trackInfo.dt = dt;                         % frame duration
trackInfo.R = 1/6*dE/dt;                           % motion blur coefficient

% load rEM parameters
numReinitialize = get(findobj('Tag', 'numTrials'), 'String');

% load pEM parameters
minStates = get(findobj('Tag', 'minSize'), 'String');
maxStates = get(findobj('Tag', 'maxSize'), 'String');
numPerturb = get(findobj('Tag', 'numPerturb'), 'String');
maxiter = get(findobj('Tag', 'maxiter'), 'String');
convergence = get(findobj('Tag', 'convergence'),'String');

set(findobj('Tag', 'numTrials'), 'Value',str2double(numReinitialize));
set(findobj('Tag', 'minSize'), 'Value',str2double(minStates));
set(findobj('Tag', 'maxSize'), 'Value',str2double(maxStates));
set(findobj('Tag', 'numPerturb'), 'Value',str2double(numPerturb));
set(findobj('Tag', 'maxiter'), 'Value',str2double(maxiter));
set(findobj('Tag', 'convergence'),'Value',str2double(convergence));

numReinitialize = get(findobj('Tag', 'numTrials'), 'Value');
minStates = get(findobj('Tag', 'minSize'), 'Value');
maxStates = get(findobj('Tag', 'maxSize'), 'Value');
numPerturb = get(findobj('Tag', 'numPerturb'), 'Value');
maxiter = get(findobj('Tag', 'maxiter'), 'Value');
convergence = get(findobj('Tag', 'convergence'),'Value');
showplot = get(findobj('Tag', 'showplot'), 'Value');
params.numPerturbation = numPerturb;   % number of perturbations trials
params.converged = convergence;       % convergence condition for EM
params.maxiter = maxiter;         % maximum number of iterations for EM
params.showplot = showplot;            % displays progress of parameter estimates (0,1)
params.verbose = 1;             % display progress on command window (0,1)

% calculate the displacements for each particle track
deltaX = cell(trackInfo.numberOfTracks,1);
for i = 1:trackInfo.numberOfTracks
    deltaX{i} = diff(X{i});
end

% calculate relevant properties to enhance compuatational time
[trackInfo.trackLength trackInfo.uniqueLength] = TrackLengthParameters(deltaX);
[trackInfo.diagonals trackInfo.correlations trackInfo.C] = CovarianceProperties(deltaX);

% diffusivity and static localization estimate from covariance-based estimator
trackInfo.D_cve = mean((trackInfo.diagonals+2*trackInfo.correlations)/(2*trackInfo.dt),2);
trackInfo.sigma_cve = mean(trackInfo.diagonals,2)/2 - trackInfo.D_cve*trackInfo.dt*(1-2*trackInfo.R); 


% BIC Model Selection Loop
results = struct;
BIC = zeros(maxStates,1); 
for numStates = minStates:maxStates
    startTime = tic;
    disp([num2str(numStates) ' state model']);

    % random initialization
    [D0 P0 S0] = RandomInitialization(numStates,trackInfo.D_cve,trackInfo.sigma_cve);

    % run rEM
    [baseD baseS baseP Lmax] = rEM(deltaX,D0,P0,S0,params,trackInfo,numReinitialize);
    
    % run pEM
    [baseD baseS baseP Lmax posteriorProb] = pEM(deltaX,baseD,baseP,baseS,Lmax,params,trackInfo);

    % calculate BIC
    BIC(numStates) = Lmax - numStates/2*log(trackInfo.numberOfTracks);

    % display results    
    disp('-------------------------------------------------------');
    disp([num2str(numStates) ' state model results:']);
    disp(['D_k = ' num2str(baseD) ' um^2/s']);
    disp(['sigma_k = ' num2str(baseS) ' um']);
    disp(['pi_k = ' num2str(baseP) ]);
    disp(['L = ' num2str(Lmax)]);
    disp(['BIC = ' num2str(BIC(numStates))]);
    disp('-------------------------------------------------------');
    
    % store results
    results(numStates).numberOfStates = numStates;
    results(numStates).BIC = BIC(numStates);
    results(numStates).optimalD = baseD;
    results(numStates).optimalS = baseS;
    results(numStates).optimalP = baseP;
    results(numStates).optimalL = Lmax;
    results(numStates).posteriorProb = posteriorProb;
    results(numStates).elapsedTime = toc(startTime);
    
    % check BIC model selection
    if numStates > 1 
        if  BIC(numStates-1) ~= 0
            if BIC(numStates) < BIC(numStates-1)
                display(['Lower BIC found.  Optimal number of State:' num2str(numStates-1)]);
                break;
            end
            if numStates == maxStates
                display(['Optimal number of states not found is larger than ' num2str(maxStates)]);
                break;
            end
        end
    end
end
[MAX numStates] = max(BIC);

% store results
data.X = X;
data.params = params;
data.trackInfo = trackInfo;
data.results = results;
data.BIC = BIC(numStates);
data.optimalD = results(numStates).optimalD;
data.optimalP = results(numStates).optimalP;
data.optimalS = results(numStates).optimalS;
data.optimalL = results(numStates).optimalL;
data.posteriorProb = results(numStates).posteriorProb;
guidata(hObject,data);

disp('Finished analysis');
disp('-------------------------------------------------------');
disp(['Optimal size: ' num2str(numStates) ' states']);
disp(['D_k = ' num2str(data.optimalD) ' um^2/s']);
disp(['sigma_k = ' num2str(data.optimalS) ' um']);
disp(['pi_k = ' num2str(data.optimalP) ]);
disp(['L = ' num2str(data.optimalL(end))]);
disp(['BIC = ' num2str(BIC(numStates))]);
disp('-------------------------------------------------------');




% --- Executes on button press in loaddata.
function loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
[filename dirpath] = uigetfile('Select data file');
[tmp name ext] = fileparts(filename);

disp(['Loaded ' filename]);
set(findobj('Tag', 'fileloaded'), 'String',filename);

status = 1;
switch ext
    case '.mat'
        load(fullfile(dirpath,filename));
        duration = zeros(length(X),1);
        for i = 1:length(X)
             duration(i) = size(X{i},1);
        end
        
    case '.txt'
        data = dlmread(fullfile(dirpath,filename));
        id = data(:,end);
        duration = zeros(id(end),1);
        X = cell(id(end),1);
        for i = 1:id(end)
            index = find(id == i);
            duration(i) = length(x);
            X{i} = data(index,1:2);
        end
    otherwise
        disp('Data is not in the right format!');
        status = 0;
end

handles.X = X; 
guidata(hObject,handles);


% 

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function numTracks_Callback(hObject, eventdata, handles)
% hObject    handle to numTracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numTracks as text
%        str2double(get(hObject,'String')) returns contents of numTracks as a double


% --- Executes during object creation, after setting all properties.
function numTracks_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numTracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function duration_Callback(hObject, eventdata, handles)
% hObject    handle to duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of duration as text
%        str2double(get(hObject,'String')) returns contents of duration as a double


% --- Executes during object creation, after setting all properties.
function duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function minSize_Callback(hObject, eventdata, handles)
% hObject    handle to minSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minSize as text
%        str2double(get(hObject,'String')) returns contents of minSize as a double


% --- Executes during object creation, after setting all properties.
function minSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxSize_Callback(hObject, eventdata, handles)
% hObject    handle to maxSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxSize as text
%        str2double(get(hObject,'String')) returns contents of maxSize as a double


% --- Executes during object creation, after setting all properties.
function maxSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showplot.
function showplot_Callback(hObject, eventdata, handles)
% hObject    handle to showplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showplot



function numPerturb_Callback(hObject, eventdata, handles)
% hObject    handle to numPerturb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numPerturb as text
%        str2double(get(hObject,'String')) returns contents of numPerturb as a double


% --- Executes during object creation, after setting all properties.
function numPerturb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numPerturb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxiter_Callback(hObject, eventdata, handles)
% hObject    handle to maxiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxiter as text
%        str2double(get(hObject,'String')) returns contents of maxiter as a double


% --- Executes during object creation, after setting all properties.
function maxiter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function convergence_Callback(hObject, eventdata, handles)
% hObject    handle to convergence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of convergence as text
%        str2double(get(hObject,'String')) returns contents of convergence as a double


% --- Executes during object creation, after setting all properties.
function convergence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to convergence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot2.
function plot2_Callback(hObject, eventdata, handles)
% hObject    handle to plot2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp(['Plotting Posterior Tracks']);

data = guidata(hObject);
X = data.X;
posteriorProb = data.posteriorProb;
DisplayPosteriorTracks(X,posteriorProb);


function savefile_Callback(hObject, eventdata, handles)
% hObject    handle to savefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savefile as text
%        str2double(get(hObject,'String')) returns contents of savefile as a double


% --- Executes during object creation, after setting all properties.
function savefile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveresults.
function saveresults_Callback(hObject, eventdata, handles)
% hObject    handle to saveresults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hint: get(hObject,'Value') returns toggle state of saveresults
data = guidata(hObject);
filename = get(findobj('Tag', 'fileloaded'));
[tmp name] = fileparts(filename);
savename = [name '_Results'];

savefolder = 'Results';
if exist(savefolder) ~= 7
    mkdir(savefolder);
end
results = data.results;
save(fullfile(savefolder,savename),'results');
disp(['Saving Results to: ' savefolder filesep savename]);


% --- Executes on button press in display.
function display_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
data = guidata(hObject);

% display results    
disp('-------------------------------------------------------');
disp(['Final Results:']);
disp(['D_k = ' num2str(data.optimalD) ' um^2/s']);
disp(['sigma_k = ' num2str(data.optimalS) ' um']);
disp(['pi_k = ' num2str(data.optimalP) ]);
disp(['L = ' num2str(data.optimalL)]);
disp(['BIC = ' num2str(data.BIC)]);
disp('-------------------------------------------------------');


% --- Executes on button press in plot1.
function plot1_Callback(hObject, eventdata, handles)
% hObject    handle to plot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Plotting MSD for each state');

data = guidata(hObject);
X = data.X;
posteriorProb = data.posteriorProb;
trackInfo = data.trackInfo;
numLags = 10;
DisplayWeightedMSD(X,posteriorProb,numLags,trackInfo.dt);

function dim_Callback(hObject, eventdata, handles)
% hObject    handle to dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dim as text
%        str2double(get(hObject,'String')) returns contents of dim as a double


% --- Executes during object creation, after setting all properties.
function dim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numTrials_Callback(hObject, eventdata, handles)
% hObject    handle to numTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numTrials as text
%        str2double(get(hObject,'String')) returns contents of numTrials as a double


% --- Executes during object creation, after setting all properties.
function numTrials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
