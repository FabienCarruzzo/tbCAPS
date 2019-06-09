%% This is the script of what is supposed to be a friendly user interface
% for application of co-activation pattern analysis. 

function varargout = CAP_TBsmall(varargin)
% CAP_TBsmall MATLAB code for CAP_TBsmall.fig
%      CAP_TBsmall, by itself, creates a new CAP_TBsmall or raises the existing
%      singleton*.
%
%      H = CAP_TBsmall returns the handle to a new CAP_TBsmall or the handle to
%      the existing singleton*.
%
%      CAP_TBsmall('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAP_TBsmall.M with the given input arguments.
%
%      CAP_TBsmall('Property','Value',...) creates a new CAP_TBsmall or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CAP_TBsmall_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CAP_TBsmall_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CAP_TBsmall

% Last Modified by GUIDE v2.5 03-Mar-2017 13:41:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CAP_TBsmall_OpeningFcn, ...
                   'gui_OutputFcn',  @CAP_TBsmall_OutputFcn, ...
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

% Executes when the window opens; all handle variables are created at this
% time
function CAP_TBsmall_OpeningFcn(hObject, eventdata, handles, varargin)

%%%%%%%%%%%%%%%%%%%%
% Path and other miscellaneous settings

% Adds the paths to the subfolders of the toolbox that will be important
% for the plotting and the analysis
addpath(genpath('./Plotting'));
addpath(genpath('./Analysis'));
addpath(genpath('./DefaultData'));

% Sets warnings off
warning('off');

% Choose default command line output for CAP_TBsmall
handles.output = hObject;



%%%%%%%%%%%%%%%%%%%%
% Data loading

% TC will contain the time courses of the subjects from the different
% populations
handles.TC = {};

% FD contains the traces of framewise displacement for the subjects (n_TP x
% n_subj per cell of the array, one cell per dataset)
handles.FD = {};

% Information on the NIFTI files from which the data originate
handles.brain_info = {};

% Mask used prior to CAP analysis
handles.mask = {};

% Number of datasets added to the interface. A dataset is defined as a
% population of subjects from the same experimental group (e.g. an ensemble 
% of subjects suffering from the same disorder)
handles.n_datasets = 0;

% Stores the number of subjects that have been loaded
handles.n_subjects = {};

% SubjNames contains the names of the files from which subject data have
% been sampled (full paths)
handles.SubjNames = {};

% MotName contains the name(s) of the file(s) loaded as motion ones
handles.MotName = {};

% TP and VOX contain the number of time points (of frames) and of brain
% voxels that are present in the loaded datasets. Those values are
% initialized at -inf, and then take the values of the first file that is
% being loaded if that file looks reasonable dimensionally speaking. In the
% scripts below, it is assumed that all subject populations loaded have a
% similar number of time points and of voxels
handles.SubjSize.TP = -inf;
handles.SubjSize.VOX = -inf;

% By default, the reference population from which CAPs will be extracted
% will be the first uploaded one
handles.ReferencePopulation = 1;

% Loads and sets the brain underlay used for plotting purposes
Underlay = load_nii('Underlay.nii');
Underlay_mat = [Underlay.hdr.hist.srow_x; Underlay.hdr.hist.srow_y; Underlay.hdr.hist.srow_z; 0 0 0 1];
Underlay_dim = Underlay.hdr.dime.dim;
Underlay_dim = Underlay_dim(2:4);
handles.Underlay_info.dim = Underlay_dim;
handles.Underlay_info.mat = Underlay_mat;
clear Underlay
clear Underlay_dim
clear Underlay_mat
load('brain.mat');
assignin('base','brain', brain);
handles.brain = brain;
clear brain

% Handles for the TR
handles.TR = -inf;
handles.isTROK = false;

%%%%%%%%%%%%%%%%%%%%
% Seed selection and seed maps

% Seed used for the analysis
handles.seed = [];
handles.seed2 = [];

% One map per subject
handles.SeedMaps = {};

% One average map throughout subjects
handles.AvgSeedMap = [];

%%%%%%%%%%%%%%%%%%%%
% Time points selection

% Motion threshold for scrubbing
handles.Tmot = 0.5;

% Threshold for frame selection in the analysis
handles.T = 0.5;

% Sets the right text header in front of the frame selection threshold box
% (threshold or retention percentage)
if get(handles.TRadio,'Value')
    set(handles.TText,'String','T [-]');
    handles.SelMode = 'Threshold';
else
    set(handles.TText,'String','P [%]');
    handles.SelMode = 'Percentage';
end



% Activation and deactivation frames kept for all datasets
handles.Xonp = {};
handles.Xonn = {};

% Percentage of frames retained for CAP analysis (discarding both the
% baseline time points and the scrubbed time points)
handles.RetainedPercentage = {};

% Indices of the frames that have been retained (i.e. when do they occur in
% the full time course), of baseline frames, and of scrubbed frames
handles.FrameIndices = {};

%%%%%%%%%%%%%%%%%%%%
% CAP analysis

% Number of times that clustering is run
handles.n_rep = 20;

% Percentage voxels to keep for clustering (positive - Pp - and negative -
% Pn - ones)
handles.Pp = 100;
handles.Pn = 100;

% Number of clusters to use in the analysis
handles.K = 5;

% Type of CAPs computed: can take a value of 'Act','Deact' or 'Both'
handles.CAPType = '';

% Indices of the CAP to which frames from the reference population and from
% the other populations are assigned
handles.idx = {};

% Value of correlation of the control group frame that is the Tper-th least
% close to its CAP
handles.CorrDist = [];

% Contains the CAPs
handles.CAP = [];

% Contains the standard deviation for the CAPs
handles.STDCAP = [];

% Parameters for the GMM
handles.Gamma_GMM = [];
handles.Priors_GMM = [];
handles.Mu_GMM = [];
% This one must be put into a structure to have several sparse matrices
% filling it
handles.Sigma_GMM = {};

% Percentile threshold used in frame assignment
handles.percentile = 5;

%%%%%%%%%%%%%%%%%%%%
% Metrics

% Will contain the metrics
% State matrix (n_subjects x n_time points)
handles.TPM = {};

% State counts (raw and frac)
handles.Counts = {};

% Number of times entering a state
handles.Number = {};

% Average duration within a state
handles.Avg_Duration = {};

% Duration of all the excursions within a state
handles.Duration = {};

% Transition probabilities 
handles.TM = {};

% Cumulative sum of states
handles.TPMCum = {};

%%%%%%%%%%%%%%%%%%%%
% General utilities

% Log containing the different events summoned from the toolbox
handles.Log = {};

% Colors used in plotting of all populations
handles.PopColor{1} = [1,0.9,0.4; 0.8, 1.0, 1.0; 0.2, 1, 0.2; 0.9, 0.9, 0.9];
handles.PopColor{2} = [1, 0, 0.2; 0.2, 0.2, 0.8; 0, 0.4, 0; 0, 0, 0];

% Project title, by default 'Untitled'
handles.project_title = 'Untitled';

% Directory to which data is to be saved (initially loaded as ./SavedData)
handles.savedir = fullfile(pwd,'SavedData');
set(handles.SaveFolderText,'String',handles.savedir);

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = CAP_TBsmall_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% SECTION 1: LOADING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data Button Click

% Executes when adding a subject population (clicking on 'A. Load data')
function DataButton_Callback(hObject, eventdata, handles)
% Opens up a menu to choose the required files for the analysis; the user
% must select four files:
% 1. Data file
% 2. Mask file
% 3. Info file (header of NIFTI)
% 4. Motion file
%
% He can select them in the order he likes
[filename1,pathname1]=uigetfile({'*.*','All Files'},...
  'Select data, motion, mask and brain info files...','MultiSelect','on');

% If the user has indeed entered files
if ~isequal(filename1,0) || ~isequal(pathname1,0)
    % There should be four selected files. In this switch, we test
    % for the amount of entered files
    if length(filename1) == 4

        % The files are loaded sequentially
        for i = 1:length(filename1)
            File{i} = fullfile(pathname1, filename1{i});
            tmp = load(File{i});
            assignin('base','tmp', tmp);
            tmp = struct2array(tmp);
            
            % Finds what type of file DataType is between the four
            % possibilities
            DataType = CAP_FindDataType(tmp);
            
            % Accordingly, fill the right handle with the information
            switch DataType
                case 'Data'
                    % We store the data into handles.TC and the file name that goes
                    % with it
                    handles.TC{handles.n_datasets+1} = tmp;
                    
                    % Takes only the last two parts of the file name and
                    % puts them in tmp_file
                    [tmp_file,n_delim] = strsplit(File{i},'/');
                    
                    if isempty(n_delim)
                        tmp_file = strsplit(File{i},'\');
                    end
                    
                    tmp_file = tmp_file(end-1:end);
                    
                    % This is what is saved and displayed in the main
                    % window then
                    handles.SubjNames{handles.n_datasets+1} = fullfile(tmp_file{1},tmp_file{2});
                    handles.n_subjects{handles.n_datasets+1} = size(handles.TC{handles.n_datasets+1},2);
                    
                    % Some commands are run only for the first dataset that we add
                    if handles.n_datasets == 0
                            % We compute and store the number of voxels and the number of time
                            % points, as well as the number of subjects
                            handles.SubjSize.VOX = size(handles.TC{1}{1},2);
                            handles.SubjSize.TP = size(handles.TC{1}{1},1);
                    end
                    
                    % Sets the text label about data dimensions
                    set(handles.Dimensionality_Text, 'String', [num2str(handles.SubjSize.TP),...
                        ' frames x ',num2str(handles.SubjSize.VOX),' voxels (',...
                        strjoin(arrayfun(@(x) num2str(x),cell2mat(handles.n_subjects),...
                        'UniformOutput',false),'+'),')']);
        
                case 'Motion'  
                    % We store the path of the motion file added
                    handles.MotName{handles.n_datasets+1} = File{i};
                    % If the dimensions hold, we store the file into the FD variable
                    % and then plot the FD ribbon graph
                    handles.FD{handles.n_datasets+1} = tmp; 
                    
                case 'Mask'
                    handles.mask{handles.n_datasets+1} = tmp;
                case 'Info'
                    % If so, we store the value and we validate the choice
                    handles.brain_info{handles.n_datasets+1} = tmp;

                % If the data file is unknown, then we return an error and
                % the user must enter files again
                case 'Unknown'
                    errordlg('At least one of the selected files is not recognized; please try again !');
                    handles = ClearDataButton_Callback(handles.ClearDataButton, eventdata, handles);
            end
            
        end
        
        % Check if the dimensionality of the entered data holds between
        % the file types. It may be that the user entered four files of
        % the same type (e.g. four data files), rather than one of each
        % type as required
        [is_DataOK,Data_problems] = CAP_IsDataOK(handles.TC{handles.n_datasets+1},handles.FD{handles.n_datasets+1},...
                handles.mask{handles.n_datasets+1},handles.brain_info{handles.n_datasets+1});
        if is_DataOK

            % We increment handles.n_datasets
            handles.n_datasets = handles.n_datasets + 1;

            % We fill the list of loaded populations, and make it visible
            handles = FillPopulationList(handles.RefPop,handles);
            set(handles.RefPop,'Visible','on');

            % We also fill the list of subjects from the seed menu, but do
            % not make it visible yet
            handles = FillSubjectList(handles.SubjectMenu,handles);

            % We can now enable the seed selection
            set(handles.SeedButton,'Enable','on');

            % Also, we can now color the button in green
            set(hObject,'BackgroundColor', [0.4 0.6 0.4]);
            
            % If we are loading the first dataset, we convert the underlay
            % to the resolution of the functional data for plotting
            if handles.n_datasets == 1

                % The brain variable now contains a good resolution
                % underlay that can directly be overlapped with the
                % functional data
                handles.brain = CAP_V2V(handles.brain,handles.Underlay_info.dim,...
                    handles.Underlay_info.mat,handles.brain_info{1}.dim,handles.brain_info{1}.mat);

            elseif handles.n_datasets > 1 && handles.n_datasets < 5 && ~isempty(handles.CAP)
                set(handles.AssignButton,'Enable','on');
            elseif handles.n_datasets > 4
                errordlg('Please enter at most four different populations in the interface !');
                handles = ClearDataButton_Callback(handles.ClearDataButton, eventdata, handles);
            end
            
            handles.Log = CAP_AddToLog(handles.Log,'Data correctly loaded');

        % If it doesn't hold, then we return an error
        else
            errordlg(['There is a dimensionality problem in your data: ',Data_problems]);
            handles = ClearDataButton_Callback(handles.ClearDataButton, eventdata, handles);
        end 
         
    % If a different number of files is entered, then there is a problem,
    % and everything is reset
    else
            errordlg('You did not enter the correct number of files !');
            handles = ClearDataButton_Callback(handles.ClearDataButton, eventdata, handles);
    end
% Else, an error is displayed and the user is prompted to enter files
else
    errordlg('Cancelling data entry will not solve your problems !');
    handles = ClearDataButton_Callback(handles.ClearDataButton, eventdata, handles);
end


% Update handles structure
guidata(hObject, handles);

%% TR Textbox Interaction

% Executes when we go to the TR field to add the TR of the experiment
function TR_Entry_Callback(hObject, eventdata, handles)

% If the TR takes a reasonable value, then we validate it
if (~isempty(str2double(get(hObject,'String')))) && ...
        (str2double(get(hObject,'String')) > 0.5) && ...
        (str2double(get(hObject,'String')) <= 5)  
    
    handles.TR = str2double(get(hObject,'String'));    
    set(hObject,'BackgroundColor', [0.4 0.6 0.4]);
    handles.isTROK = true;
    
    handles.Log = CAP_AddToLog(handles.Log,'Correct value of TR entered',{handles.TR},{'TR'});

% Else, the TR value is not accepted
else 
    set(hObject,'BackgroundColor', [0.93 0.84 0.84]);
    handles.isTROK = false;  
end

guidata(hObject, handles); 

% Executes during creation of the TR textbox
function handles = TR_Entry_CreateFcn(hObject, eventdata, handles)

set(hObject,'Enable','off');
set(hObject,'String','Click to enter...');
set(hObject,'FontAngle','italic');

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','r');
end

guidata(hObject, handles); 

% Executes when clicking on the TR text space
function TR_Entry_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to TMotEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Enable','on');
set(hObject,'String','');
set(hObject,'FontAngle','normal');
uicontrol(hObject);
guidata(hObject, handles); 

%% Seed sliders interactions
% For the below functions, the goal is to change the value of the slider
% textboxes when the sliders are moved, and to update the graph display
% accordingly. For this purpose, cla is used to clear graph content prior
% to a new display

% --- Executes on slider movement.
function SliderX_Callback(hObject, eventdata, handles)

cla(handles.SeedGraphX);
set(handles.XCoordText,'String',['X: ',sprintf('%.2f',get(hObject,'Value'))]);
handles.SeedGraphX = plot_slice(handles.seed,get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'X',get(hObject,'Value'),handles.SeedGraphX);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SliderX_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
guidata(hObject, handles); 

% --- Executes on slider movement.
function SliderY_Callback(hObject, eventdata, handles)

cla(handles.SeedGraphY);
set(handles.YCoordText,'String',['Y: ',sprintf('%.2f',get(hObject,'Value'))]);
handles.SeedGraphY = plot_slice(handles.seed,get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Y',get(hObject,'Value'),handles.SeedGraphY);
guidata(hObject, handles); 

% --- Executes during object creation, after setting all properties.
function SliderY_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
guidata(hObject, handles); 

% --- Executes on slider movement.
function SliderZ_Callback(hObject, eventdata, handles)

cla(handles.SeedGraphZ);
set(handles.ZCoordText,'String',['Z: ',sprintf('%.2f',get(hObject,'Value'))]);
handles.SeedGraphZ = plot_slice(handles.seed,get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Z',get(hObject,'Value'),handles.SeedGraphZ);
guidata(hObject, handles); 

% --- Executes during object creation, after setting all properties.
function SliderZ_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
guidata(hObject, handles); 


%% Seed Button Controls
% We want to define what happens when loading data (SeedButton) or when
% attempting to plot them (PlotSeedButton)

% Executes when clicking on 'B. Select a seed'
function SeedButton_Callback(hObject, eventdata, handles)

% No seed union asked
if ~get(handles.Union_Checkbox,'Value')
    
    [filename_seed,pathname_seed]=uigetfile({'*.*','All Files'},...
      'Select Seed File...');

    File_seed = fullfile(pathname_seed, filename_seed);
                tmp = load(File_seed);
                assignin('base','tmp', tmp);
                tmp = struct2array(tmp);
% Seed union asked
else
    
    [filename_seed,pathname_seed]=uigetfile({'*.*','All Files'},...
      'Select Seed File...','MultiSelect','on');

    File_seed = fullfile(pathname_seed, filename_seed{1});
    tmp = load(File_seed);
    assignin('base','tmp', tmp);
    tmp = struct2array(tmp);
    handles.seed = tmp;
    
    File_seed = fullfile(pathname_seed, filename_seed{2});
    tmp2 = load(File_seed);
    assignin('base','tmp2', tmp2);
    tmp2 = struct2array(tmp2);
    handles.seed2 = tmp2;
end

% If the user has indeed entered files
if ~isequal(filename_seed,0) || ~isequal(pathname_seed,0) 
    
    % If the file is of suitable dimensions
    if islogical(tmp) && size(tmp,2) == 1 && size(tmp,1) == sum(handles.mask{1})
        
        % Then we put it in the handles, enable the plotting button, and
        % make the seed selection button green
        handles.seed = tmp;
        set(handles.PlotSeedButton,'Enable','on');
        set(handles.SeedButton,'BackgroundColor', [0.4 0.6 0.4]);
        
        % We can now go through the next parts of the analysis, so we
        % enable the related buttons
        set(handles.TPSelectionButton,'Enable','on');
        set(handles.SeedMapPushButton,'Enable','on');
        
        handles.Log = CAP_AddToLog(handles.Log,'Seed chosen',{File_seed},{'Seed file'});
    else
        errordlg('The file you entered appears to be of wrong dimensions...');
    end
    
else
    errordlg('Please enter a seed file !');
end
    
guidata(hObject, handles);

% Change in the status of the checkbox for seed union
function Union_Checkbox_Callback(hObject, eventdata, handles)

handles.Log = CAP_AddToLog(handles.Log,'Seed union status changed',{get(hObject,'Value')},{'Status'});

guidata(hObject, handles);

% Executes when clicking on 'Plot Seed'
function PlotSeedButton_Callback(hObject, eventdata, handles)

% Clears the present graph content
cla(handles.SeedGraphX);
cla(handles.SeedGraphX);
cla(handles.SeedGraphX);

% Plots the slices within the graph windows
handles.SeedGraphX = plot_slice(handles.seed,get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'X',get(handles.SliderX,'Value'),handles.SeedGraphX);
handles.SeedGraphY = plot_slice(handles.seed,get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Y',get(handles.SliderY,'Value'),handles.SeedGraphY);
handles.SeedGraphZ = plot_slice(handles.seed,get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Z',get(handles.SliderZ,'Value'),handles.SeedGraphZ);

% Sets the sliders to visible
set(handles.SliderX,'Visible','on');
set(handles.SliderY,'Visible','on');
set(handles.SliderZ,'Visible','on');

% Sets the text values at the ones of the sliders
set(handles.XCoordText,'String',['X: ',sprintf('%.2f',get(handles.SliderX,'Value'))]);
set(handles.YCoordText,'String',['Y: ',sprintf('%.2f',get(handles.SliderY,'Value'))]);
set(handles.ZCoordText,'String',['Z: ',sprintf('%.2f',get(handles.SliderZ,'Value'))]);

% Sets the visibility of the slider texts to on
set(handles.XCoordText,'Visible','on');
set(handles.YCoordText,'Visible','on');
set(handles.ZCoordText,'Visible','on');

handles.Log = CAP_AddToLog(handles.Log,'Seed plots activated');

guidata(hObject, handles);

%% Reference group List Control
% If we change the entry from this list, we change the reference group on
% which seed maps will be computed, and CAPs calculated

% Executes on selection change in RefPop
function RefPop_Callback(hObject, eventdata, handles)

% We take the actual value of the population as the reference group
handles.ReferencePopulation = get(hObject,'Value');

% We want to fill the subject lists again, according to how many subjects
% now lie in the considered reference population
handles = FillSubjectList(handles.SubjectMenu,handles);
handles = FillSubjectList(handles.SubjectMenuMetrics,handles);

% We want to reset the graph displays and only enable what can be computed

% Resetting the seed map section
set(handles.SubjectMenu,'Visible','off');
set(handles.SubjectMenu,'Value',1);

handles = ResetGraphDisplay(handles.SeedMapX,handles);
handles = ResetGraphDisplay(handles.SeedMapY,handles);
handles = ResetGraphDisplay(handles.SeedMapZ,handles);

handles = ResetGraphDisplay(handles.SubjSeedMapX,handles);
handles = ResetGraphDisplay(handles.SubjSeedMapY,handles);
handles = ResetGraphDisplay(handles.SubjSeedMapZ,handles);

set(handles.SeedMap_SliderX,'Visible','off');
set(handles.SeedMapSliderX,'Visible','off');
set(handles.SeedMap_SliderY,'Visible','off');
set(handles.SeedMapSliderY,'Visible','off');
set(handles.SeedMap_SliderZ,'Visible','off');
set(handles.SeedMapSliderZ,'Visible','off');

handles = ResetGraphDisplay(handles.ColorbarSeed,handles);
set(handles.TSeed,'Visible','off');
set(handles.TSeed_Slider,'Visible','off');

% Resetting the time points selection section
handles = ResetGraphDisplay(handles.TPViolin,handles);

% Resetting the CAP analysis section
set(handles.ClusterButton,'Enable','off');
set(handles.GMM_Button,'Enable','off');

set(handles.TVIS,'Visible','off');
set(handles.TVIS_Slider,'Visible','off');

handles = ResetGraphDisplay(handles.ColorbarCAP,handles);

set(handles.CAP1_Frames,'Visible','off');
set(handles.CAP2_Frames,'Visible','off');
set(handles.CAP3_Frames,'Visible','off');
set(handles.CAP4_Frames,'Visible','off');
set(handles.CAP5_Frames,'Visible','off');
set(handles.CAP6_Frames,'Visible','off');

tmp_CAPX = {handles.CAP1X,handles.CAP2X,handles.CAP3X,handles.CAP4X,handles.CAP5X,handles.CAP6X};
tmp_CAPY = {handles.CAP1Y,handles.CAP2Y,handles.CAP3Y,handles.CAP4Y,handles.CAP5Y,handles.CAP6Y};
tmp_CAPZ = {handles.CAP1Z,handles.CAP2Z,handles.CAP3Z,handles.CAP4Z,handles.CAP5Z,handles.CAP6Z};

for i = 1:6
    handles = ResetGraphDisplay(tmp_CAPX{i},handles);
    handles = ResetGraphDisplay(tmp_CAPY{i},handles);
    handles = ResetGraphDisplay(tmp_CAPZ{i},handles);
end

set(handles.CAP_XC,'Visible','off');
set(handles.CAP_YC,'Visible','off');
set(handles.CAP_ZC,'Visible','off');

set(handles.CAP_SliderX,'Visible','off');
set(handles.CAP_SliderY,'Visible','off');
set(handles.CAP_SliderZ,'Visible','off');

% Resetting the metrics section
set(handles.MetricsButton,'Enable','off');

handles = ResetGraphDisplay(handles.CAP_Mat,handles);
handles = ResetGraphDisplay(handles.TMGraph,handles);
handles = ResetGraphDisplay(handles.TM_Subject,handles);
handles = ResetGraphDisplay(handles.ColorbarSimMat,handles);
handles = ResetGraphDisplay(handles.ColorbarTransMat,handles);
handles = ResetGraphDisplay(handles.DynStates,handles);
handles = ResetGraphDisplay(handles.CumStates,handles);
handles = ResetGraphDisplay(handles.ViolinCounts,handles);
handles = ResetGraphDisplay(handles.ViolinCountsFrac,handles);
handles = ResetGraphDisplay(handles.ViolinNumber,handles);
handles = ResetGraphDisplay(handles.ViolinDuration,handles);

set(handles.SubjectMenuMetrics,'Visible','off');
set(handles.SubjectMenuMetrics,'Value',1);

set(handles.StateMenu,'Visible','off');

handles.Log = CAP_AddToLog(handles.Log,'Reference population change',{handles.ReferencePopulation},{'New reference population index'});

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function RefPop_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%SECTION 2: SEED MAPS COMPUTATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Seed Map Computation Button
% When pressing on this button, classical seed maps are computed for the
% population of subjects chosen as the reference one in the loading part.
% The last entered seed is used

function SeedMapPushButton_Callback(hObject, eventdata, handles)

% Computes seed maps for each subject and for the population, using the
% data from the chosen reference population
[handles.SeedMaps,handles.AvgSeedMap] = CAP_Compute_SeedMap(handles.TC{handles.ReferencePopulation},handles.seed);

% Graphical displays

% Making the plots, texts and sliders visible
set(handles.SubjectMenu,'Visible','on');
set(handles.SeedMapX,'Visible','on');
set(handles.SeedMapY,'Visible','on');
set(handles.SeedMapZ,'Visible','on');
set(handles.SubjSeedMapX,'Visible','on');
set(handles.SubjSeedMapY,'Visible','on');
set(handles.SubjSeedMapZ,'Visible','on');
set(handles.SeedMap_SliderX,'Visible','on');
set(handles.SeedMap_SliderY,'Visible','on');
set(handles.SeedMap_SliderZ,'Visible','on');
set(handles.SeedMapSliderX,'Visible','on');
set(handles.SeedMapSliderY,'Visible','on');
set(handles.SeedMapSliderZ,'Visible','on');
set(handles.TSeed_Slider,'Visible','on');
set(handles.TSeed,'Visible','on');
set(handles.ColorbarSeed,'Visible','on');

% Writing down the text with current MNI coordinates
set(handles.SeedMap_SliderX,'String',['X: ',sprintf('%.2f',get(handles.SeedMapSliderX,'Value'))]);
set(handles.SeedMap_SliderY,'String',['Y: ',sprintf('%.2f',get(handles.SeedMapSliderY,'Value'))]);
set(handles.SeedMap_SliderZ,'String',['Z: ',sprintf('%.2f',get(handles.SeedMapSliderZ,'Value'))]);

% Clears previous plot contents (in case we want to re-plot after changing
% the seed)
cla(handles.SeedMapX);
cla(handles.SeedMapY);
cla(handles.SeedMapZ);
cla(handles.SubjSeedMapX);
cla(handles.SubjSeedMapY);
cla(handles.SubjSeedMapZ);

% Plots new slices
handles.SeedMapX = plot_slice(handles.AvgSeedMap,0.25,1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'X',get(handles.SeedMapSliderX,'Value'),handles.SeedMapX);
handles.SeedMapY = plot_slice(handles.AvgSeedMap,0.25,1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Y',get(handles.SeedMapSliderY,'Value'),handles.SeedMapY);
handles.SeedMapZ = plot_slice(handles.AvgSeedMap,0.25,1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Z',get(handles.SeedMapSliderZ,'Value'),handles.SeedMapZ);

% Plots subject specific slices after having selected the subject
handles.SubjSeedMapX = plot_slice(handles.SeedMaps{get(handles.SubjectMenu,'Value')},get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'X',get(handles.SeedMapSliderX,'Value'),handles.SubjSeedMapX);
handles.SubjSeedMapY = plot_slice(handles.SeedMaps{get(handles.SubjectMenu,'Value')},get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Y',get(handles.SeedMapSliderY,'Value'),handles.SubjSeedMapY);
handles.SubjSeedMapZ = plot_slice(handles.SeedMaps{get(handles.SubjectMenu,'Value')},get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Z',get(handles.SeedMapSliderZ,'Value'),handles.SubjSeedMapZ);

% Adds the colorbar for the seed maps (between -1 and 1)
handles.ColorbarSeed = Create_CAP_colorbar(-1,1,0.5,get(handles.TSeed_Slider,'Value'),'',handles.ColorbarSeed,'Horizontal');

handles.Log = CAP_AddToLog(handles.Log,'Seed maps displayed');

guidata(hObject,handles);

%% Slider Controls (MNI coordinates)
% We want to reload the seed images with the new parameters when changing a
% slider, so we clear the previous display, change the text summarizing the
% MNI coordinate where we stand, and plot the new image

function SeedMapSliderX_Callback(hObject, eventdata, handles)

% Clears graphs
cla(handles.SeedMapX);
cla(handles.SubjSeedMapX);

% Changes slider texts
set(handles.SeedMap_SliderX,'String',['X: ',sprintf('%.2f',get(hObject,'Value'))]);

% Plots new slices
handles.SeedMapX = plot_slice(handles.AvgSeedMap,get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'X',get(hObject,'Value'),handles.SeedMapX);
handles.SubjSeedMapX = plot_slice(handles.SeedMaps{get(handles.SubjectMenu,'Value')},get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'X',get(hObject,'Value'),handles.SubjSeedMapX);
guidata(hObject, handles);

function SeedMapSliderX_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function SeedMapSliderY_Callback(hObject, eventdata, handles)

% Clears graphs
cla(handles.SeedMapY);
cla(handles.SubjSeedMapY);

% Changes slider texts
set(handles.SeedMap_SliderY,'String',['Y: ',sprintf('%.2f',get(hObject,'Value'))]);

% Plots new slices
handles.SeedMapY = plot_slice(handles.AvgSeedMap,get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Y',get(hObject,'Value'),handles.SeedMapY);
handles.SubjSeedMapY = plot_slice(handles.SeedMaps{get(handles.SubjectMenu,'Value')},get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Y',get(hObject,'Value'),handles.SubjSeedMapY);
guidata(hObject, handles);

function SeedMapSliderY_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function SeedMapSliderZ_Callback(hObject, eventdata, handles)

cla(handles.SeedMapZ);
cla(handles.SubjSeedMapZ);

% Changes slider texts
set(handles.SeedMap_SliderZ,'String',['Z: ',sprintf('%.2f',get(hObject,'Value'))]);

% Plots new slices
handles.SeedMapZ = plot_slice(handles.AvgSeedMap,get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Z',get(hObject,'Value'),handles.SeedMapZ);
handles.SubjSeedMapZ = plot_slice(handles.SeedMaps{get(handles.SubjectMenu,'Value')},get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Z',get(hObject,'Value'),handles.SubjSeedMapZ);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SeedMapSliderZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SeedMapSliderZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%% Slider controls (visualization threshold)
% We want to replot and modify the colorbar according to the visualization
% threshold that we select (and also change the text)

function TSeed_Slider_Callback(hObject, eventdata, handles)

% Clears previous plot contents
cla(handles.SeedMapX);
cla(handles.SeedMapY);
cla(handles.SeedMapZ);
cla(handles.SubjSeedMapX);
cla(handles.SubjSeedMapY);
cla(handles.SubjSeedMapZ);

% Plots new slices (average seed maps)
handles.SeedMapX = plot_slice(handles.AvgSeedMap,get(hObject,'Value'),1,...
    handles.mask{handles.ReferencePopulation},handles.brain,...
    handles.brain_info{handles.ReferencePopulation},'X',...
    get(handles.SeedMapSliderX,'Value'),handles.SeedMapX);

handles.SeedMapY = plot_slice(handles.AvgSeedMap,get(hObject,'Value'),1,...
    handles.mask{handles.ReferencePopulation},handles.brain,...
    handles.brain_info{handles.ReferencePopulation},'Y',...
    get(handles.SeedMapSliderY,'Value'),handles.SeedMapY);

handles.SeedMapZ = plot_slice(handles.AvgSeedMap,get(hObject,'Value'),1,...
    handles.mask{handles.ReferencePopulation},handles.brain,...
    handles.brain_info{handles.ReferencePopulation},'Z',...
    get(handles.SeedMapSliderZ,'Value'),handles.SeedMapZ);

% Plots subject specific slices after having selected the subject
handles.SubjSeedMapX = plot_slice(handles.SeedMaps{get(handles.SubjectMenu,'Value')},...
    get(hObject,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,...
    handles.brain_info{handles.ReferencePopulation},'X',get(handles.SeedMapSliderX,...
    'Value'),handles.SubjSeedMapX);

handles.SubjSeedMapY = plot_slice(handles.SeedMaps{get(handles.SubjectMenu,'Value')},...
    get(hObject,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,...
    handles.brain_info{handles.ReferencePopulation},'Y',get(handles.SeedMapSliderY,...
    'Value'),handles.SubjSeedMapY);

handles.SubjSeedMapZ = plot_slice(handles.SeedMaps{get(handles.SubjectMenu,'Value')},...
    get(hObject,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,...
    handles.brain_info{handles.ReferencePopulation},'Z',get(handles.SeedMapSliderZ,...
    'Value'),handles.SubjSeedMapZ);

% Modifies the text
set(handles.TSeed,'String',['Tv: ',sprintf('%.2f',get(hObject,'Value'))]);

% Clears and replots the colorbar
cla(handles.ColorbarSeed);
handles.ColorbarSeed = Create_CAP_colorbar(-1,1,0.5,get(hObject,'Value'),'',handles.ColorbarSeed,'Horizontal');

guidata(hObject,handles);

function TSeed_Slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%% Subject Menu Control
% We want that when we change the value of the subject menu, the displays
% from the graphs change too (if the graphs already contain something)

function SubjectMenu_Callback(hObject, eventdata, handles)

% If the seed maps have already been computed (i.e. if we have already
% pressed on the computation button), then we adjust the displays of the
% subject-specific maps
try
    cla(handles.SubjSeedMapX);
    cla(handles.SubjSeedMapY);
    cla(handles.SubjSeedMapZ);
    handles.SubjSeedMapX = plot_slice(handles.SeedMaps{get(handles.SubjectMenu,'Value')},get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'X',get(handles.SeedMapSliderX,'Value'),handles.SubjSeedMapX);
    handles.SubjSeedMapY = plot_slice(handles.SeedMaps{get(handles.SubjectMenu,'Value')},get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Y',get(handles.SeedMapSliderY,'Value'),handles.SubjSeedMapY);
    handles.SubjSeedMapZ = plot_slice(handles.SeedMaps{get(handles.SubjectMenu,'Value')},get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Z',get(handles.SeedMapSliderZ,'Value'),handles.SubjSeedMapZ);
catch
    errordlg('Please recompute seed maps for the presently considered population !');
end

handles.Log = CAP_AddToLog(handles.Log,'Changed displayed subjectwise seed map',{get(handles.SubjectMenu,'Value')},{'New subject index'});

guidata(hObject, handles);

function SubjectMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%PART 3: CAP ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Motion parameter entry

function TMotEdit_Callback(hObject, eventdata, handles)

% If we enter a reasonable value, it is taken as a new threshold
if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 0) && (str2double(get(hObject,'String')) <= 0.5)
    handles.Tmot = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [0.4 0.6 0.4]);
    
    handles.Log = CAP_AddToLog(handles.Log,'Valid motion threshold value entered',{handles.Tmot},{'Motion threshold value'});
    
% If we set something wrong again, we set the threshold value back to the
% default of 0.5
else
    set(hObject,'BackgroundColor', [0.93 0.84 0.84]);
    handles.Tmot = 0.5;
end

guidata(hObject, handles); 

% When clicking on the motion button
function handles = TMotEdit_ButtonDownFcn(hObject, eventdata, handles)

set(hObject,'Enable','on');
set(hObject,'String','');
set(hObject,'FontAngle','normal');
uicontrol(hObject);

guidata(hObject, handles); 

% When the object is created
function handles = TMotEdit_CreateFcn(hObject, eventdata, handles)

set(hObject,'Enable','off');
set(hObject,'String','Click to enter...');
set(hObject,'FontAngle','italic');

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','r');
end
guidata(hObject, handles); 

%% Frame selection parameter entry

function TEdit_Callback(hObject, eventdata, handles)

% If we enter a reasonable value, it is taken as the new threshold
if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 0) && (str2double(get(hObject,'String')) <= 100)
    handles.T = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [0.4 0.6 0.4]);
    
    handles.Log = CAP_AddToLog(handles.Log,'Valid (de)activation threshold entered',{handles.T},{'Threshold value'});
    
% If we set something wrong again, we set the threshold value back to the
% default of 0.5
else
    set(hObject,'BackgroundColor',[0.93 0.84 0.84]);
    handles.T = 0.5;
end

guidata(hObject, handles); 

% When the object is created
function handles = TEdit_CreateFcn(hObject, eventdata, handles)

set(hObject,'Enable','off');
set(hObject,'String','Click to enter...');
set(hObject,'FontAngle','italic');

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles); 

% When clicking on it
function TEdit_ButtonDownFcn(hObject,eventdata,handles)

set(hObject,'Enable','on');
set(hObject,'String','');
set(hObject,'FontAngle','normal');
uicontrol(hObject);

guidata(hObject,handles);

%% Frame selection mode control buttons

% If we select threshold, we update the Selmode accordingly
function TRadio_Callback(hObject, eventdata, handles)

set(handles.TText,'String','T [-]');
handles.SelMode = 'Threshold';

handles.Log = CAP_AddToLog(handles.Log,'Changed time points selection scheme',{handles.SelMode},{'Selected mode'});

guidata(hObject,handles);

% Same for percentage
function PRadio_Callback(hObject, eventdata, handles)

set(handles.TText,'String','P [%]');
handles.SelMode = 'Percentage';

handles.Log = CAP_AddToLog(handles.Log,'Changed time points selection scheme',{handles.SelMode},{'Selected mode'});

guidata(hObject,handles);

%% Time points selection control
% When clicking on the select time points button, the frames matching the
% provided thresholds for scrubbing and for frame retention are selected.
% Both activation and deactivation frames are selected. This is performed
% on all the loaded populations of subjects

% Upon clicking on the 'Select time points' button
function TPSelectionButton_Callback(hObject, eventdata, handles)

% Clears the current plot display (for the case of having already computed
% something before with other parameters)
cla(handles.TPViolin);

% Performs the analysis to extract frames of activity for all loaded
% populations (done for each dataset)
for n_ds = 1:handles.n_datasets
    % Xonp and Xonn contain the frames (deactivation frames have been
    % switched in sign, so that deactivation is positive)
    if ~get(handles.Union_Checkbox,'Value')
        [handles.Xonp{n_ds},handles.Xonn{n_ds},p,Indices] = CAP_find_activity(handles.TC{n_ds},handles.seed,handles.T,handles.FD{n_ds},handles.Tmot,handles.SelMode);
    else
        [handles.Xonp{n_ds},handles.Xonn{n_ds},p,Indices] = CAP_find_activity_twoseed(handles.TC{n_ds},handles.seed,handles.seed2,handles.T,handles.FD{n_ds},handles.Tmot);
    end
    % Percentage of retained frames for both activation and deactivation
    % cases across subjects
    handles.RetainedPercentage{n_ds} = p(4:5,:);

    % Indices of the frames that have been retained (used later for metrics
    % computations)
    handles.FrameIndices{n_ds} = Indices; 
end

% Enables to go to the next step of the analysis and cluster the extracted
% frames
set(handles.ClusterButton,'Enable','on');

tmp_toplot = ConcatMat(handles.RetainedPercentage,handles.n_datasets,2,handles.n_subjects,'FD');

% Displays the violin plot of subject scrubbing percentage for the
% reference population
[~,~,handles.TPViolin] = MakeViolin(tmp_toplot,handles.TPViolin,{'A','D'},'Frames ret. [%]',handles.PopColor,handles.n_datasets,2);
set(handles.TPViolin,'Visible','on');

handles.Log = CAP_AddToLog(handles.Log,'Time points selected',{['1 to ',num2str(handles.n_datasets)],handles.SelMode,handles.T,handles.Tmot},{'Datasets indices','Selection mode','Activation threshold','Motion threshold'});

guidata(hObject, handles);


%% Parameter control: number of clusters to use

function ClusterEdit_Callback(hObject, eventdata, handles)

if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 1) && (str2double(get(hObject,'String')) <= 12)
    handles.K = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [0.4 0.6 0.4]);

    handles.Log = CAP_AddToLog(handles.Log,'Valid number of clusters chosen',{handles.K},{'Number of clusters'});
    
else
    set(hObject,'BackgroundColor',[0.93 0.84 0.84]);
    handles.K = 5;
end

guidata(hObject, handles);

function handles = ClusterEdit_CreateFcn(hObject, eventdata, handles)

set(hObject,'Enable','off');
set(hObject,'String','Click to enter...');
set(hObject,'FontAngle','italic');

guidata(hObject, handles);

function ClusterEdit_ButtonDownFcn(hObject,eventdata,handles)

set(hObject,'Enable','on');
set(hObject,'String','');
set(hObject,'FontAngle','normal');
uicontrol(hObject);

guidata(hObject,handles);

%% Parameter control: number of replicates of the k-means algorithm

function ClusterRepEdit_Callback(hObject, eventdata, handles)

if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 0) && (str2double(get(hObject,'String')) <= 50)
    handles.n_rep = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [0.4 0.6 0.4]);

    handles.Log = CAP_AddToLog(handles.Log,'Valid number of replicates chosen',{handles.n_rep},{'Number of replicates'});
    
else
    set(hObject,'BackgroundColor',[0.93 0.84 0.84]);
    handles.n_rep = 20;
end

guidata(hObject, handles);

function handles = ClusterRepEdit_CreateFcn(hObject, eventdata, handles)

set(hObject,'Enable','off');
set(hObject,'String','Click to enter...');
set(hObject,'FontAngle','italic');

guidata(hObject, handles);

function ClusterRepEdit_ButtonDownFcn(hObject, eventdata, handles)

set(hObject,'Enable','on');
set(hObject,'String','');
set(hObject,'FontAngle','normal');
uicontrol(hObject);

guidata(hObject,handles);

%% Parameter control: percentage of positive voxels on which to cluster

function ClusterPpEdit_Callback(hObject, eventdata, handles)

if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 0) && (str2double(get(hObject,'String')) <= 100)
    handles.Pp = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [0.4 0.6 0.4]);

    handles.Log = CAP_AddToLog(handles.Log,'Valid percentage positive voxels chosen',{handles.Pp},{'Percentage positive voxels'});
    
else
    set(hObject,'BackgroundColor',[0.93 0.84 0.84]);
    handles.Pp = 20;
end

guidata(hObject, handles);

function handles = ClusterPpEdit_CreateFcn(hObject, eventdata, handles)

set(hObject,'Enable','off');
set(hObject,'String','Click to enter...');
set(hObject,'FontAngle','italic');

guidata(hObject, handles);

function ClusterPpEdit_ButtonDownFcn(hObject, eventdata, handles)

set(hObject,'Enable','on');
set(hObject,'String','');
set(hObject,'FontAngle','normal');
uicontrol(hObject);

guidata(hObject,handles);

%% Parameter control: percentage of negative voxels on which to cluster

function ClusterPnEdit_Callback(hObject, eventdata, handles)

if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 0) && (str2double(get(hObject,'String')) <= 100)
    handles.Pn = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [0.4 0.6 0.4]);
    
    handles.Log = CAP_AddToLog(handles.Log,'Valid percentage negative voxels chosen',{handles.Pn},{'Percentage negative voxels'});
    
else
    set(hObject,'BackgroundColor',[0.93 0.84 0.84]);
    handles.Pn = 20;
end

guidata(hObject, handles);

function handles = ClusterPnEdit_CreateFcn(hObject, eventdata, handles)

set(hObject,'Enable','off');
set(hObject,'String','Click to enter...');
set(hObject,'FontAngle','italic');

guidata(hObject, handles);

function ClusterPnEdit_ButtonDownFcn(hObject, eventdata, handles)

set(hObject,'Enable','on');
set(hObject,'String','');
set(hObject,'FontAngle','normal');
uicontrol(hObject);

guidata(hObject,handles);

%% Clustering control
% When pressing on the 'Cluster' button, we want to run clustering for the
% specified mode (Activation frames, Deactivation frames, or both types of
% frames together), using the previously declared parameters

% Upon clicking on 'Cluster'
function ClusterButton_Callback(hObject, eventdata, handles)

% Determines which type of clustering we want to perform, and runs
% clustering accordingly

% Clustering activation frames
if get(handles.ActRadio,'Value')
    
    handles.CAPType = 'Act';
    [handles.CAP,~,handles.STDCAP,handles.idx{handles.ReferencePopulation},...
        handles.CorrDist] = Run_Clustering(cell2mat(handles.Xonp{handles.ReferencePopulation}),...
        handles.K,handles.mask{handles.ReferencePopulation},handles.brain_info{handles.ReferencePopulation},...
        handles.Pp,handles.Pn,handles.n_rep);
    
% Clustering deactivation frames
elseif get(handles.DeactRadio,'Value')
    
    handles.CAPType = 'Deact';
    [handles.CAP,~,handles.STDCAP,handles.idx{handles.ReferencePopulation},...
        handles.CorrDist] = Run_Clustering(cell2mat(handles.Xonn{handles.ReferencePopulation}),...
        handles.K,handles.mask{handles.ReferencePopulation},handles.brain_info{handles.ReferencePopulation},...
        handles.Pp,handles.Pn,handles.n_rep);
    
% Clustering both activation and deactivation frames
elseif get(handles.BothFramesRadio,'Value')
    
    handles.CAPType = 'Both';
    
    % In this last case, we want the activation and deactivation frames to
    % be of opposite sign, so we invert the sign of Xonn
    [handles.CAP,~,handles.STDCAP,handles.idx{handles.ReferencePopulation},...
        handles.CorrDist] = Run_Clustering([cell2mat(handles.Xonp{handles.ReferencePopulation}),...
        (-1)*cell2mat(handles.Xonn{handles.ReferencePopulation})],...
        handles.K,handles.mask{handles.ReferencePopulation},handles.brain_info{handles.ReferencePopulation},...
        handles.Pp,handles.Pn,handles.n_rep);
    
% Normally we should NEVER enter this
else
    errordlg('I have no idea how you got there, but yes, you made everything crash in a way I do not get... Congrats...');
end

% Makes the sliders visible, and the related text too (CAP MNI coordinates)
set(handles.CAP_SliderX,'Visible','on');
set(handles.CAP_SliderY,'Visible','on');
set(handles.CAP_SliderZ,'Visible','on');
set(handles.CAP_XC,'Visible','on');
set(handles.CAP_YC,'Visible','on');
set(handles.CAP_ZC,'Visible','on');
set(handles.CAP_XC,'String',['X: ',sprintf('%.2f',get(handles.CAP_SliderX,'Value'))]);
set(handles.CAP_YC,'String',['Y: ',sprintf('%.2f',get(handles.CAP_SliderY,'Value'))]);
set(handles.CAP_ZC,'String',['Z: ',sprintf('%.2f',get(handles.CAP_SliderZ,'Value'))]);

% Same for the slider for the visualization threshold
set(handles.TVIS,'Visible','on');
set(handles.TVIS_Slider,'Visible','on'); 
set(handles.TVIS,'String',['Tv: ',sprintf('%.2f',get(handles.TVIS_Slider,'Value'))]);

% Makes the colorbar for the CAPs visible
handles.ColorbarCAP = Create_CAP_colorbar(-1.5,1.5,0.5,get(handles.TVIS_Slider,'Value'),'',handles.ColorbarCAP,'Horizontal');
set(handles.ColorbarCAP,'Visible','on');

% Concatenates all CAP information into metavariables for easier subsequent
% changes
tmpX = {handles.CAP1X,handles.CAP2X,handles.CAP3X,handles.CAP4X,handles.CAP5X,handles.CAP6X};
tmpY = {handles.CAP1Y,handles.CAP2Y,handles.CAP3Y,handles.CAP4Y,handles.CAP5Y,handles.CAP6Y};
tmpZ = {handles.CAP1Z,handles.CAP2Z,handles.CAP3Z,handles.CAP4Z,handles.CAP5Z,handles.CAP6Z};
tmpF = {handles.CAP1_Frames,handles.CAP2_Frames,handles.CAP3_Frames,handles.CAP4_Frames,handles.CAP5_Frames,handles.CAP6_Frames};

% For each CAP...
for i_CAP = 1:min([handles.K,6])
    
    % Clears the display for each dimension
    cla(tmpX{i_CAP});
    cla(tmpY{i_CAP});
    cla(tmpZ{i_CAP});
    
    % Plots the new slice for each dimension
    tmpX{i_CAP} = plot_slice(handles.CAP(i_CAP,:)./handles.STDCAP(:,i_CAP)',...
        get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},...
        'X',get(handles.CAP_SliderX,'Value'),tmpX{i_CAP});
    
    tmpY{i_CAP} = plot_slice(handles.CAP(i_CAP,:)./handles.STDCAP(:,i_CAP)',...
        get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},...
        'Y',get(handles.CAP_SliderY,'Value'),tmpY{i_CAP});
    
    tmpZ{i_CAP} = plot_slice(handles.CAP(i_CAP,:)./handles.STDCAP(:,i_CAP)',...
        get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},...
        'Z',get(handles.CAP_SliderZ,'Value'),tmpZ{i_CAP});

    % Sets the frame percentage text visible and at the right value (number
    % of frames from a CAP/total frame number, and then percentage that it
    % stands for)
    set(tmpF{i_CAP},'Visible','on');
    set(tmpF{i_CAP},'String',{[num2str(sum(handles.idx{handles.ReferencePopulation}==i_CAP)),'/',...
        num2str(size(handles.idx{handles.ReferencePopulation},1))],[sprintf('%.2f',...
        sum(handles.idx{handles.ReferencePopulation}==i_CAP)/size(handles.idx{handles.ReferencePopulation},1)*100),' [%]']});
end

% Fills that subject menu with the subjects from the reference population
handles = FillSubjectList(handles.SubjectMenuMetrics,handles);

% Also enables and fills the state menu
handles = FillStateList(handles.StateMenu,handles);


% Enables the Metrics button for the next part of the analysis if we
% only deal with one dataset
if handles.n_datasets == 1
    set(handles.MetricsButton,'Enable','on');
% Else, we enable the assignment before enabling the metrics computation
elseif handles.n_datasets > 1
    set(handles.AssignButton,'Enable','on');
end

handles.Log = CAP_AddToLog(handles.Log,'Clustering performed',...
    {handles.ReferencePopulation,handles.K,handles.n_rep,handles.Pp,...
    handles.Pn,handles.CAPType},{'Reference group index',...
    'Number of clusters','Number of replicates',...
    'Percentage positive voxels','Percentage negative voxels','Type of CAPs'});

% When we have done one k-means clustering run and thus have 'meaningful'
% brain states, we are allowed to go for the GMM option
set(handles.GMM_Button,'Enable','on');

guidata(hObject, handles);

%% Clustering with the more elaborate GMM (Gaussian Mixture Model) method
function GMM_Button_Callback(hObject, eventdata, handles)

if handles.K ~= size(handles.CAP,1)
    errordlg('Please choose a matching value of K with respect to the latest run k-means trial!');
end

% Initializes GMM parameters with the results from k-means
for i = 1:handles.K
    handles.Mu_GMM(:,i) = handles.CAP(i,:)';
    handles.Sigma_GMM{i} = sparse(1:handles.SubjSize.VOX,1:handles.SubjSize.VOX,handles.STDCAP(:,i));
    handles.Priors_GMM(i) = 1/handles.K;
end

% handles.Mu_GMM = 2*rand(1000,4)-ones(1000,4);
% handles.Sigma_GMM = rand(1000,1000,4);

% Clustering activation frames
if get(handles.ActRadio,'Value')
    
    handles.CAPType = 'Act';
    
    [handles.Gamma_GMM,handles.Priors_GMM,handles.Mu_GMM,handles.Sigma_GMM,handles.CorrDist,handles.idx{handles.ReferencePopulation}] =...
        CAP_GMM(cell2mat(handles.Xonp{handles.ReferencePopulation}),handles.K,...
        handles.Priors_GMM,...
        handles.Mu_GMM,...
        handles.Sigma_GMM,...
        handles.mask{handles.ReferencePopulation},handles.brain_info{handles.ReferencePopulation},...
        handles.Pp,handles.Pn);
    
% Clustering deactivation frames
elseif get(handles.DeactRadio,'Value')
    
    handles.CAPType = 'Deact';
    
    [handles.Gamma_GMM,handles.Priors_GMM,handles.Mu_GMM,handles.Sigma_GMM,handles.CorrDist,handles.idx{handles.ReferencePopulation}] =...
        CAP_GMM(cell2mat(handles.Xonn{handles.ReferencePopulation}),handles.K,...
        handles.Priors_GMM,...
        handles.Mu_GMM,...
        handles.Sigma_GMM,...
        handles.mask{handles.ReferencePopulation},handles.brain_info{handles.ReferencePopulation},...
        handles.Pp,handles.Pn);
    
% Clustering both activation and deactivation frames
elseif get(handles.BothFramesRadio,'Value')
    
    handles.CAPType = 'Both';
    
    % In this last case, we want the activation and deactivation frames to
    % be of opposite sign, so we invert the sign of Xonn
    [handles.Gamma_GMM,handles.Priors_GMM,handles.Mu_GMM,handles.Sigma_GMM,handles.CorrDist,handles.idx{handles.ReferencePopulation}] =...
        CAP_GMM([cell2mat(handles.Xonp{handles.ReferencePopulation}),...
        (-1)*cell2mat(handles.Xonn{handles.ReferencePopulation})],handles.K,...
        handles.Priors_GMM,...
        handles.Mu_GMM,...
        handles.Sigma_GMM,...
        handles.mask{handles.ReferencePopulation},handles.brain_info{handles.ReferencePopulation},...
        handles.Pp,handles.Pn);
    
% Normally we should NEVER enter this
else
    errordlg('I have no idea how you got there, but yes, you made everything crash in a way I do not get... Congrats...');
end

% Calculates 'CAP equivalents' from the GMM outcomes
handles.CAP = handles.Mu_GMM';

for i = 1:handles.K
    handles.STDCAP(:,i) = sqrt(diag(handles.Sigma_GMM{i}));
end
    
% Makes the sliders visible, and the related text too (CAP MNI coordinates)
set(handles.CAP_SliderX,'Visible','on');
set(handles.CAP_SliderY,'Visible','on');
set(handles.CAP_SliderZ,'Visible','on');
set(handles.CAP_XC,'Visible','on');
set(handles.CAP_YC,'Visible','on');
set(handles.CAP_ZC,'Visible','on');
set(handles.CAP_XC,'String',['X: ',sprintf('%.2f',get(handles.CAP_SliderX,'Value'))]);
set(handles.CAP_YC,'String',['Y: ',sprintf('%.2f',get(handles.CAP_SliderY,'Value'))]);
set(handles.CAP_ZC,'String',['Z: ',sprintf('%.2f',get(handles.CAP_SliderZ,'Value'))]);

% Same for the slider for the visualization threshold
set(handles.TVIS,'Visible','on');
set(handles.TVIS_Slider,'Visible','on'); 
set(handles.TVIS,'String',['Tv: ',sprintf('%.2f',get(handles.TVIS_Slider,'Value'))]);

% Makes the colorbar for the CAPs visible
handles.ColorbarCAP = Create_CAP_colorbar(-1.5,1.5,0.5,get(handles.TVIS_Slider,'Value'),'',handles.ColorbarCAP,'Horizontal');
set(handles.ColorbarCAP,'Visible','on');

% Concatenates all CAP information into metavariables for easier subsequent
% changes
tmpX = {handles.CAP1X,handles.CAP2X,handles.CAP3X,handles.CAP4X,handles.CAP5X,handles.CAP6X};
tmpY = {handles.CAP1Y,handles.CAP2Y,handles.CAP3Y,handles.CAP4Y,handles.CAP5Y,handles.CAP6Y};
tmpZ = {handles.CAP1Z,handles.CAP2Z,handles.CAP3Z,handles.CAP4Z,handles.CAP5Z,handles.CAP6Z};
tmpF = {handles.CAP1_Frames,handles.CAP2_Frames,handles.CAP3_Frames,handles.CAP4_Frames,handles.CAP5_Frames,handles.CAP6_Frames};

% For each CAP...
for i_CAP = 1:min([handles.K,6])
    
    % Clears the display for each dimension
    cla(tmpX{i_CAP});
    cla(tmpY{i_CAP});
    cla(tmpZ{i_CAP});
    
    % Plots the new slice for each dimension
    tmpX{i_CAP} = plot_slice(handles.CAP(i_CAP,:)./handles.STDCAP(:,i_CAP)',...
        get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},...
        'X',get(handles.CAP_SliderX,'Value'),tmpX{i_CAP});
    
    tmpY{i_CAP} = plot_slice(handles.CAP(i_CAP,:)./handles.STDCAP(:,i_CAP)',...
        get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},...
        'Y',get(handles.CAP_SliderY,'Value'),tmpY{i_CAP});
    
    tmpZ{i_CAP} = plot_slice(handles.CAP(i_CAP,:)./handles.STDCAP(:,i_CAP)',...
        get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},...
        'Z',get(handles.CAP_SliderZ,'Value'),tmpZ{i_CAP});

    % Sets the frame percentage text visible and at the right value (number
    % of frames from a CAP/total frame number, and then percentage that it
    % stands for)
    set(tmpF{i_CAP},'Visible','on');
    set(tmpF{i_CAP},'String',{[num2str(sum(handles.idx{handles.ReferencePopulation}==i_CAP)),'/',...
        num2str(size(handles.idx{handles.ReferencePopulation},1))],[sprintf('%.2f',...
        sum(handles.idx{handles.ReferencePopulation}==i_CAP)/size(handles.idx{handles.ReferencePopulation},1)*100),' [%]']});
end

% Fills that subject menu with the subjects from the reference population
handles = FillSubjectList(handles.SubjectMenuMetrics,handles);

% Also enables and fills the state menu
handles = FillStateList(handles.StateMenu,handles);

% Enables the Metrics button for the next part of the analysis if we
% only deal with one dataset
if handles.n_datasets == 1
    set(handles.MetricsButton,'Enable','on');
% Else, we enable the assignment before enabling the metrics computation
elseif handles.n_datasets > 1
    set(handles.AssignButton,'Enable','on');
end

handles.Log = CAP_AddToLog(handles.Log,'Clustering performed',...
    {handles.ReferencePopulation,handles.K,handles.n_rep,handles.Pp,...
    handles.Pn,handles.CAPType},{'Reference group index',...
    'Number of clusters','Number of replicates',...
    'Percentage positive voxels','Percentage negative voxels','Type of CAPs'});

guidata(hObject, handles);



%% Frame assignment control
% This button is only enabled after clustering has been performed on the
% reference population. It assigns frames from the other populations to the
% computed CAPs

% Happens upon clicking on the 'Assign' buttons
function AssignButton_Callback(hObject, eventdata, handles)

tmp_notref = [];
tmp_computedTPsel = [];

% For each non-reference dataset...
for n_ds = 1:handles.n_datasets
    if n_ds ~= handles.ReferencePopulation
        
        tmp_notref = [tmp_notref,n_ds];
        
        % Attempts to access the frames for a given dataset; if it fails, it
        % means we must compute activity. If it works, we do nothing because
        % activity has already been computed
        try
            justtotest = handles.Xonp{n_ds};
        catch
            if get(handles.Union_Checkbox,'Value')
                [handles.Xonp{n_ds},handles.Xonn{n_ds},p,handles.FrameIndices{n_ds}] = ...
                    CAP_find_activity_twoseed(handles.TC{n_ds},handles.seed,handles.seed2,handles.T,handles.FD{n_ds},handles.Tmot);
            else
                [handles.Xonp{n_ds},handles.Xonn{n_ds},p,handles.FrameIndices{n_ds}] = ...
                    CAP_find_activity(handles.TC{n_ds},handles.seed,handles.T,handles.FD{n_ds},handles.Tmot,handles.SelMode);
            end
            handles.RetainedPercentage{n_ds} = p(4:5,:);
            
            tmp_computedTPsel = [tmp_computedTPsel,n_ds];
        end

        try
            % Assigning activation frames
            if strcmp(handles.CAPType,'Act') && get(handles.ActRadio,'Value')
                handles.idx{n_ds} = CAP_AssignFrames(handles.CAP,cell2mat(handles.Xonp{n_ds}),handles.CorrDist,handles.percentile)';

            % Assigning deactivation frames
            elseif strcmp(handles.CAPType,'Deact') && get(handles.DeactRadio,'Value')
                handles.idx{n_ds} = CAP_AssignFrames(handles.CAP,cell2mat(handles.Xonn{n_ds}),handles.CorrDist,handles.percentile)';

            % Assigning both activation and deactivation frames
            elseif strcmp(handles.CAPType,'Both') && get(handles.BothFramesRadio,'Value')
                handles.idx{n_ds} = CAP_AssignFrames(handles.CAP,[cell2mat(handles.Xonp{n_ds}),cell2mat((-1)*handles.Xonn{n_ds})],handles.CorrDist,handles.percentile)';

            % Normally we should NEVER enter this
            else
                errordlg('I have no idea how you got there, but yes, you made everything crash in a way I do not get... Congrats...');
            end
        catch
            errordlg('You computed CAPs with a different CAP type compared to the one used now; please use the same CAP type !');
        end
    end
end

% We then enable the computation of metrics
set(handles.MetricsButton,'Enable','on');

handles.Log = CAP_AddToLog(handles.Log,'Frame assignment performed',...
    {handles.ReferencePopulation,num2str(tmp_computedTPsel),num2str(tmp_notref),...
    handles.CAPType},{'Reference group index','Group indices for which frames were computed',...
    'Group indices for which frames were assigned','Type of CAPs'});

guidata(hObject, handles);

%% Parameter control: percentile to use for frame assignment
% This asks for the percentile to use in frame assignment (i.e. the
% threshold of correlation below which frames are left unassigned)

function Percentile_Edit_Callback(hObject, eventdata, handles)

if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 0) && (str2double(get(hObject,'String')) <= 100)
    handles.percentile = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [0.4 0.6 0.4]);

    handles.Log = CAP_AddToLog(handles.Log,'Valid percentile chosen',{handles.percentile},{'Percentile'});
    
else
    set(hObject,'BackgroundColor',[0.93 0.84 0.84]);
    handles.percentile = 5;
end

guidata(hObject, handles);

function handles = Percentile_Edit_CreateFcn(hObject, eventdata, handles)

set(hObject,'Enable','off');
set(hObject,'String','Click to enter...');
set(hObject,'FontAngle','italic');

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Percentile_Edit_ButtonDownFcn(hObject, eventdata, handles)

set(hObject,'Enable','on');
set(hObject,'String','');
set(hObject,'FontAngle','normal');
uicontrol(hObject);

guidata(hObject, handles);


%% Sliders for CAP visualization (MNI coordinates)
% When changing along a slider, we want to update the graphs and the text of
% the MNI coordinate below the slider

% X dimension slider
function CAP_SliderX_Callback(hObject, eventdata, handles)

set(handles.CAP_XC,'String',['X: ',sprintf('%.2f',get(hObject,'Value'))]);
tmp_struct = {handles.CAP1X,handles.CAP2X,handles.CAP3X,handles.CAP4X,handles.CAP5X,handles.CAP6X};

for i_CAP = 1:min([handles.K,6])
    cla(tmp_struct{i_CAP});  
    tmp_struct{i_CAP} = plot_slice(handles.CAP(i_CAP,:),get(handles.TVIS_Slider,'Value'),...
        1.5,handles.mask{handles.ReferencePopulation},handles.brain,...
        handles.brain_info{handles.ReferencePopulation},'X',get(hObject,'Value'),tmp_struct{i_CAP});
end

guidata(hObject, handles); 

function CAP_SliderX_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% Y dimension slider
function CAP_SliderY_Callback(hObject, eventdata, handles)

set(handles.CAP_YC,'String',['Y: ',sprintf('%.2f',get(hObject,'Value'))]);
tmp_struct = {handles.CAP1Y,handles.CAP2Y,handles.CAP3Y,handles.CAP4Y,handles.CAP5Y,handles.CAP6Y};

for i_CAP = 1:min([handles.K,6])
    cla(tmp_struct{i_CAP});
    tmp_struct{i_CAP} = plot_slice(handles.CAP(i_CAP,:),get(handles.TVIS_Slider,'Value'),...
        1.5,handles.mask{handles.ReferencePopulation},handles.brain,...
        handles.brain_info{handles.ReferencePopulation},'Y',get(hObject,'Value'),tmp_struct{i_CAP});
end

guidata(hObject,handles);

function CAP_SliderY_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% Z dimension slider
function CAP_SliderZ_Callback(hObject, eventdata, handles)

set(handles.CAP_ZC,'String',['Z: ',sprintf('%.2f',get(hObject,'Value'))]);
tmp_struct = {handles.CAP1Z,handles.CAP2Z,handles.CAP3Z,handles.CAP4Z,handles.CAP5Z,handles.CAP6Z};

for i_CAP = 1:min([handles.K,6])
   
    cla(tmp_struct{i_CAP});
    tmp_struct{i_CAP} = plot_slice(handles.CAP(i_CAP,:),get(handles.TVIS_Slider,'Value'),...
        1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Z',get(hObject,'Value'),tmp_struct{i_CAP});
end

guidata(hObject,handles);

function CAP_SliderZ_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%% Sliders for threshold visualization (CAP analysis)
% Again, we want to update the slices and the text if we change those
% sliders

function TVIS_Slider_Callback(hObject, eventdata, handles)

% The text is changed
set(handles.TVIS,'String',['Tv: ',sprintf('%.2f',get(hObject,'Value'))]);

% The colorbar graph is modified to suit the new threshold value
cla(handles.ColorbarCAP);
handles.ColorbarCAP = Create_CAP_colorbar(-1.5,1.5,0.5,get(hObject,'Value'),'',handles.ColorbarCAP,'Horizontal');

% The brain slices are replotted
tmpX = {handles.CAP1X,handles.CAP2X,handles.CAP3X,handles.CAP4X,handles.CAP5X,handles.CAP6X};
tmpY = {handles.CAP1Y,handles.CAP2Y,handles.CAP3Y,handles.CAP4Y,handles.CAP5Y,handles.CAP6Y};
tmpZ = {handles.CAP1Z,handles.CAP2Z,handles.CAP3Z,handles.CAP4Z,handles.CAP5Z,handles.CAP6Z};

for i_CAP = 1:min([handles.K,6])
    
    cla(tmpX{i_CAP});
    cla(tmpY{i_CAP});
    cla(tmpZ{i_CAP});
    
    tmpX{i_CAP} = plot_slice(handles.CAP(i_CAP,:),get(hObject,'Value'),1.5,...
        handles.mask{handles.ReferencePopulation},handles.brain,...
        handles.brain_info{handles.ReferencePopulation},'X',get(handles.CAP_SliderX,'Value'),tmpX{i_CAP});
    
    tmpY{i_CAP} = plot_slice(handles.CAP(i_CAP,:),get(hObject,'Value'),1.5,...
        handles.mask{handles.ReferencePopulation},handles.brain,...
        handles.brain_info{handles.ReferencePopulation},'Y',get(handles.CAP_SliderY,'Value'),tmpY{i_CAP});
    
    tmpZ{i_CAP} = plot_slice(handles.CAP(i_CAP,:),get(hObject,'Value'),1.5,...
        handles.mask{handles.ReferencePopulation},handles.brain,...
        handles.brain_info{handles.ReferencePopulation},'Z',get(handles.CAP_SliderZ,'Value'),tmpZ{i_CAP});
end

guidata(hObject,handles);

function TVIS_Slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PART 4: METRICS COMPUTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Metrics computation control
% When pressing on the 'Compute metrics' button, the different metrics for
% CAP analysis are computed, including:
% - Similarity between CAPs
% - Transition probabilities from a state to the other (average + subject)
% - Sequence of states for each subject
% - Cumulative state sequence for all subjects
% - Counts (number of frames in a state)
% - Number of times entering a state, and duration spent in a state

function MetricsButton_Callback(hObject, eventdata, handles)

% All the metrics are computed for all the datasets
for n_ds = 1:handles.n_datasets
    if n_ds == handles.ReferencePopulation
        tmp_nclust = handles.K;
    else
        tmp_nclust = handles.K+1;
    end
    
    try
        [handles.TPM{n_ds},handles.Counts{n_ds},...
            handles.Number{n_ds},handles.Avg_Duration{n_ds},...
            handles.Duration{n_ds},handles.TM{n_ds}] =...
            Compute_Metrics(handles.idx{n_ds},handles.FrameIndices{n_ds}.kept.active,...
            handles.FrameIndices{n_ds}.kept.deactive,handles.FrameIndices{n_ds}.scrubbed,...
            tmp_nclust,handles.TR,handles.CAPType);   
    catch
        
        errordlg('You tried computing metrics using parameter values different from the ones that were employed to generate CAPs; please check !');
    end
end

handles.Log = CAP_AddToLog(handles.Log,'Metrics computed',...
    {handles.n_datasets,handles.K,handles.TR,handles.CAPType},...
    {'Number of datasets','Number of clusters','TR','CAP type'});

% 1. Display of the similarity matrix between CAPs

% Computation of the similarity
SimMat = corr(handles.CAP',handles.CAP');

% Graph set visible, and plotting
handles = ResetGraphDisplay(handles.CAP_Mat,handles);
set(handles.CAP_Mat,'Visible','on');
imagesc(SimMat,'Parent',handles.CAP_Mat);

% Correlation ranges from -1 to 1, so this is what we make the graph
% colorbar vary within. We also make the graph square and remove the axes
caxis(handles.CAP_Mat,[-1 1]);
axis(handles.CAP_Mat,'square','on');
axis(handles.CAP_Mat,'off');

% Addition of the colorbar just below
set(handles.ColorbarSimMat,'Visible','on');
handles.ColorbarSimMat = Create_CAP_colorbar(-1,1,0.5,0,'',...
    handles.ColorbarSimMat,'Horizontal');


% 2. Transition matrix for all subjects together

tmp_toplot = squeeze(mean(handles.TM{handles.ReferencePopulation},3));
tmp_toplot = tmp_toplot(4:end,4:end);

% Make graph visible and plotting
handles = ResetGraphDisplay(handles.TMGraph,handles);
set(handles.TMGraph,'Visible','on');
imagesc(tmp_toplot,'Parent',handles.TMGraph);

clear tmp_toplot

% Arbitrary setting of probability scale from 0 to 0.03
caxis(handles.TMGraph,[0 0.03]);
axis(handles.TMGraph,'square','on');
axis(handles.TMGraph,'off');

% 3. Transition matrix for one subject

tmp_toplot = squeeze(handles.TM{handles.ReferencePopulation}(:,:,get(handles.SubjectMenuMetrics,'Value')));
tmp_toplot = tmp_toplot(4:end,4:end);

% makes graph visible and plots the information given by the Subject popup
handles = ResetGraphDisplay(handles.TM_Subject,handles);
set(handles.TM_Subject,'Visible','on');
imagesc(tmp_toplot,...
    'Parent',handles.TM_Subject);

clear tmp_toplot

% Same setting for the axes as for the average graph
caxis(handles.TM_Subject,[0 0.03]);
axis(handles.TM_Subject,'square','on');
axis(handles.TM_Subject,'off');

% We then create the colorbar for both cases
set(handles.ColorbarTransMat,'Visible','on');
handles.ColorbarTransMat = Create_CAP_colorbar(0,0.03,0.01,0,'',...
    handles.ColorbarTransMat,'Horizontal');

% Makes the subject menu visible
set(handles.SubjectMenuMetrics,'Visible','on');

% 4. Dynamic state plotting

%%%%%%%% TO SOLVE SOMEHOW: for the moment, the states are displayed with
%%%%%%%% the classical colorbar. The problem is that using another colorbar
%%%%%%%% instead will make all the figures use that colorbar. For now, I
%%%%%%%% could not find a solution to have only that specific figure be
%%%%%%%% plotted with different colors

% Makes the graph visible
handles = ResetGraphDisplay(handles.DynStates,handles);
set(handles.DynStates,'Visible','on');

% Concatenates information from the different datasets
tmp_toplot = [];

for i = 1:handles.n_datasets
    tmp_toplot = [tmp_toplot; handles.TPM{i}; -2*ones(5,handles.SubjSize.TP)];
end
tmp_toplot = tmp_toplot(1:end-5,:);

% If the TR has been properly entered, the x-axis is time; else, it depicts
% time index. In any case, we plot the states
if handles.isTROK
    imagesc(tmp_toplot,'Parent',handles.DynStates);
    xlabel(handles.DynStates,'Time [s]');
else
    imagesc(tmp_toplot,'Parent',handles.DynStates);
    xlabel(handles.DynStates,'Time index [-]');
end

ylabel(handles.DynStates,'Subjects [-]');
axis(handles.DynStates,'off');
caxis(handles.DynStates,[-2.5,handles.K+0.5]);

clear tmp_toplot

% 5. Cumulative state distributions

% Makes the graph visible
handles = ResetGraphDisplay(handles.CumStates,handles);
set(handles.CumStates,'Visible','on');

for i = 1:handles.n_datasets
    % Cumulative distribution for the state that we want to be displayed (i.e.
    % the state from the popup menu)
    handles.TPMCum{i} = cumsum(handles.TPM{i} == get(handles.StateMenu,'Value'),2);

    % Average of the considered state across subjects
    tmp_TPMCum{i} = mean(handles.TPMCum{i},1);
end
    
% Similarly as above, we plot time if we have a valid TR; else, we plot
% 'time index'
if handles.isTROK == false

    for i = 1:handles.n_datasets
        % We first plot each subject curve
        for j = 1:size(handles.TPMCum{i},1)
            plot(handles.TPMCum{i}(j,:),'Color',handles.PopColor{1}(i,:),...
                'Parent',handles.CumStates);
            hold(handles.CumStates,'on');
        end
    end

    for i = 1:handles.n_datasets
        % Then, we plot a bold average across subjects
        plot(tmp_TPMCum{i},'Color',handles.PopColor{2}(i,:),...
            'LineWidth',2,'Parent',handles.CumStates); 
        xlabel(handles.CumStates,'Time index [-]','FontSize',8);
    end
else
    for i = 1:handles.n_datasets
        for j = 1:size(handles.TPMCum{i},1)
            plot(((1:size(handles.TPM{i},2))-1)*handles.TR,...
                handles.TPMCum{i}(j,:),...
                'Color',handles.PopColor{1}(i,:),'Parent',handles.CumStates);
            hold(handles.CumStates,'on');
        end
    end

    for i = 1:handles.n_datasets
        plot(((1:size(handles.TPM{i},2))-1)*handles.TR,...
            tmp_TPMCum{i},...
            'LineWidth',2,'Color',handles.PopColor{2}(i,:),'Parent',handles.CumStates);
        xlabel(handles.CumStates,'Time [s]','FontSize',8);
    end
end


ylabel(handles.CumStates,'Cumul. sum [-]','FontSize',8);
set(handles.CumStates,'Box','off');

% Makes the state menu visible
set(handles.StateMenu,'Visible','on');

% 6. Violin plots
% Below, we plot violin plots depicting:
% - Raw counts of state excursions
% - Fractional counts of state excursions
% - Number of times entering a state
% - Duration of state excursions

% We build the legend used to plot the violins
% leg_viol = cell(handles.K);
for i = 1:handles.K
    leg_viol{i} = num2str(i);
end

% Makes graphs ready
handles = ResetGraphDisplay(handles.ViolinCounts,handles);
set(handles.ViolinCounts,'Visible','on');

handles = ResetGraphDisplay(handles.ViolinCountsFrac,handles);
set(handles.ViolinCountsFrac,'Visible','on');

handles = ResetGraphDisplay(handles.ViolinNumber,handles);
set(handles.ViolinNumber,'Visible','on');

handles = ResetGraphDisplay(handles.ViolinDuration,handles);
set(handles.ViolinDuration,'Visible','on');

% Concatenates the values from the different populations
tmp_toplot = ConcatMat(handles.Counts,handles.n_datasets,handles.K,handles.n_subjects,'Raw counts');

% Plots the raw count values
[~,~,handles.ViolinCounts] = MakeViolin(tmp_toplot,...
    handles.ViolinCounts,leg_viol,'Raw counts [-]',handles.PopColor,handles.n_datasets,handles.K);

clear tmp_toplot

tmp_toplot = ConcatMat(handles.Counts,handles.n_datasets,handles.K,handles.n_subjects,'Normalized counts');

% Plots the normalized count values
[~,~,handles.ViolinCountsFrac] = MakeViolin(tmp_toplot,...
    handles.ViolinCountsFrac,leg_viol,'Norm counts [-]',handles.PopColor,handles.n_datasets,handles.K);

clear tmp_toplot

tmp_toplot = ConcatMat(handles.Number,handles.n_datasets,handles.K,handles.n_subjects,'Number');

% Plots the number of times a state is entered
[~,~,handles.ViolinNumber] = MakeViolin(tmp_toplot,...
    handles.ViolinNumber,leg_viol,'Number [-]',handles.PopColor,handles.n_datasets,handles.K);

clear tmp_toplot


% Plots the duration graph
if handles.isTROK
    
    tmp_toplot = ConcatMat(handles.Avg_Duration,handles.n_datasets,handles.K,handles.n_subjects,'Duration');

    [~,~,handles.ViolinDuration] = MakeViolin(tmp_toplot,...
        handles.ViolinDuration,leg_viol,'Dur. [s]',handles.PopColor,handles.n_datasets,handles.K);
    
    clear tmp_toplot
   
else
    errordlg('Not showing duration violin plot as TR was not entered !');
    set(handles.ViolinDuration,'Visible','off');
end
 
guidata(hObject,handles);

%% Subject popup menu control (metrics computation)
% When a new subject is chosen, the display of the transition matrix graph
% is changed

function SubjectMenuMetrics_Callback(hObject, eventdata, handles)

% In the case when we have something to plot...
try
    % ... we reset the graph display, make the graph visible, and plot
    % again
    handles = ResetGraphDisplay(handles.TM_Subject,handles);
    set(handles.TM_Subject,'Visible','on');
    
    tmp_toplot = squeeze(handles.TM{handles.ReferencePopulation}(:,:,get(handles.SubjectMenuMetrics,'Value')));
    tmp_toplot = tmp_toplot(4:end,4:end);
    
    imagesc(tmp_toplot,'Parent',handles.TM_Subject);
    caxis(handles.TM_Subject,[0 0.06]);
    axis(handles.TM_Subject,'square','on');
    axis(handles.TM_Subject,'off');
    
    clear tmp_toplot
    
    handles.Log = CAP_AddToLog(handles.Log,'Subject index changed (metrics)',...
    {get(hObject,'Value')},...
    {'Subject index'});

catch
    errordlg('Please recompute metrics for the presently considered population !');
end

guidata(hObject,handles);

function SubjectMenuMetrics_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% State menu popup control
% When we change our state of interest, we will change the display of the
% cumulative state being displayed

function StateMenu_Callback(hObject, eventdata, handles)

handles = ResetGraphDisplay(handles.CumStates,handles);
set(handles.CumStates,'Visible','on');

% In the case of a non-null matrix...
if ~isempty(handles.TPM)
    
    for i = 1:handles.n_datasets
        % Cumulative distribution for the state that we want to be displayed (i.e.
        % the state from the popup menu)
        handles.TPMCum{i} = cumsum(handles.TPM{i} == get(handles.StateMenu,'Value'),2);

        % Average of the considered state across subjects
        tmp_TPMCum{i} = mean(handles.TPMCum{i},1);
    end

    % Similarly as above, we plot time if we have a valid TR; else, we plot
    % 'time index'
    if handles.isTROK == false

        for i = 1:handles.n_datasets
            % We first plot each subject curve
            for j = 1:size(handles.TPMCum{i},1)
                plot(handles.TPMCum{i}(j,:),'Color',handles.PopColor{1}(i,:),...
                    'Parent',handles.CumStates);
                hold(handles.CumStates,'on');
            end
        end

        for i = 1:handles.n_datasets
            % Then, we plot a bold average across subjects
            plot(tmp_TPMCum{i},'Color',handles.PopColor{2}(i,:),...
                'LineWidth',2,'Parent',handles.CumStates); 
            xlabel(handles.CumStates,'Time index [-]','FontSize',8);
        end
    else
        for i = 1:handles.n_datasets
            for j = 1:size(handles.TPMCum{i},1)
                plot(((1:size(handles.TPM{i},2))-1)*handles.TR,...
                    handles.TPMCum{i}(j,:),...
                    'Color',handles.PopColor{1}(i,:),'Parent',handles.CumStates);
                hold(handles.CumStates,'on');
            end
        end

        for i = 1:handles.n_datasets
            plot(((1:size(handles.TPM{i},2))-1)*handles.TR,...
                tmp_TPMCum{i},...
                'LineWidth',2,'Color',handles.PopColor{2}(i,:),'Parent',handles.CumStates);
            xlabel(handles.CumStates,'Time [s]','FontSize',8);
        end
    end


ylabel(handles.CumStates,'Cumul. sum [-]','FontSize',8);
set(handles.CumStates,'Box','off');
    
end

guidata(hObject,handles);

function StateMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% GENERAL UTILITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Save Folder Control
% Upon clicking on the 'select save folder' button, the user can choose the
% directory where data will be saved. 

function SaveFolderButton_Callback(hObject, eventdata, handles)

% Selection of a directory
[dirname]=uigetdir('*.*','Please select a save directory');
handles.savedir = dirname;

% If the user has indeed chosen a directory, we set it as the new save
% folder
if ~isequal(dirname,0)
    set(handles.SaveFolderText,'String',handles.savedir);
    set(hObject,'BackgroundColor', [0.4 0.6 0.4]);
    
    handles.Log = CAP_AddToLog(handles.Log,'Save folder changed',...
    {handles.savedir},...
    {'Save folder'});
else
    errordlg('Please select a directory !');
end


%% Save Button Control
% Upon clicking on the 'SAVE' button, the data will be saved entirely under 
% a file name partly chosen by the user and partly depending on the present 
% date and time

function SaveButton_Callback(hObject, eventdata, handles)

% Upon pressing the save button, we want to save all the important
% information into a big matlab structure
SAVED = [];

% General information on the project
SAVED.ProjectInfo.title = handles.project_title;
SAVED.ProjectInfo.date = date;

% Name of the files that were loaded
SAVED.SubjData.SubjFileNames = handles.SubjNames;
SAVED.SubjData.MotFileNames = handles.MotName;

% Dimension over time and voxels of the files analyzed
SAVED.SubjData.Dimensions.TP = handles.SubjSize.TP;
SAVED.SubjData.Dimensions.VOX = handles.SubjSize.VOX;

% Number of subjects considered
SAVED.SubjData.n_subjects = handles.n_subjects;

% TR of the experiment
SAVED.SubjData.TR = handles.TR;

% Information about the NIFTI files used (dimensions, mapping between real
% world and index)
SAVED.BrainData.brain_info = handles.brain_info;

% Mask that was used on the considered data
SAVED.BrainData.mask = handles.mask;

% Seed used for the analysis
SAVED.BrainData.seed = handles.seed;
SAVED.BrainData.seed2 = handles.seed2;

% Motion threshold and activation threshold used in time points selection
SAVED.TPSelData.Tmot = handles.Tmot;
SAVED.TPSelData.T = handles.T;

% Type of frame selection used
SAVED.TPSelData.SelMode = handles.SelMode;

% Frames that were considered in the clustering process
SAVED.TPSelData.Act = handles.Xonp;
SAVED.TPSelData.Deact = handles.Xonn;

% Percentage frames retained for the clustering
SAVED.TPSelData.PercRetained = handles.RetainedPercentage;

% Computed seed maps (average and subject-wise)
SAVED.SeedMap.AvgMap = handles.AvgSeedMap;
SAVED.SeedMap.SubjMaps = handles.SeedMaps;

% Parameters used for clustering
SAVED.ClusterData.N = handles.n_rep;
SAVED.ClusterData.K = handles.K;
SAVED.ClusterData.Pp = handles.Pp;
SAVED.ClusterData.Pn = handles.Pn;
SAVED.ClusterData.CAPType = handles.CAPType;

% CAP data
SAVED.ClusterData.CAPs = handles.CAP;
SAVED.ClusterData.StdCAPs = handles.STDCAP;
SAVED.ClusterData.idx = handles.idx;

% GMM data
SAVED.GMMData.Mu = handles.Mu_GMM;
SAVED.GMMData.Sigma = handles.Sigma_GMM;
SAVED.GMMData.Priors = handles.Priors_GMM;
SAVED.GMMData.Gamma = handles.Gamma_GMM;

% Computed metrics
SAVED.Metrics.TPM = handles.TPM;
SAVED.Metrics.Counts = handles.Counts;
SAVED.Metrics.Number = handles.Number;
SAVED.Metrics.Avg_Duration = handles.Avg_Duration;
SAVED.Metrics.Duration = handles.Duration;
SAVED.Metrics.TM = handles.TM;

[tmp_date,tmp_date2] = strtok(date,'-');
[tmp_date2,tmp_date3] = strtok(tmp_date2,'-');
tmp_date3 = strtok(tmp_date3,'-');

% Name that will be given to the saved files
fancy_name = [handles.project_title,'_',tmp_date,'_',tmp_date2,'_',tmp_date3,'_',...
    num2str(hour(now)),'_',num2str(minute(now)),'_',...
    num2str(round(second(now)))];

% Saves NIFTI files storing the CAPs (already 'normalized'), in MNI space
CAPToNIFTI(handles.CAP./handles.STDCAP',...
    handles.mask{handles.ReferencePopulation},handles.brain_info{handles.ReferencePopulation},...
    handles.savedir,['CAP_NIFTI_',fancy_name]);

% Saves the different variables from the program
save(fullfile(handles.savedir,fancy_name),'SAVED','-v7.3');

% Adds the save process to the log
handles.Log = CAP_AddToLog(handles.Log,'Data saved');

% Writes a log .txt file with what has been done so far
file_ID = fopen(fullfile(handles.savedir,[fancy_name,'.txt']),'wt');

for i = 1:length(handles.Log)
    for j = 1:length(handles.Log{i})
        fprintf(file_ID,[handles.Log{i}{j},'\n']);
    end
    fprintf(file_ID,'\n');
end

fclose(file_ID);

% Clears the structure now that it has been saved
clear SAVED


%% Clear functions
% Executes when pressing on the 'CLEAR' button for data loading; supposed
% to set everything back to normal (when the window opened)

function handles = ClearDataButton_Callback(hObject, eventdata, handles)
     
%%%%%%%%%% Putting the loading part (top) back to normal %%%%%%%%%%%

% Makes 'A. Load data' red again
set(handles.DataButton,'BackgroundColor',[0.93 0.84 0.84]);

% Resets the time point and voxel parameters
handles.SubjSize.TP = -inf;
handles.SubjSize.VOX = -inf;

% Resets the TR
handles.TR = -inf;
handles.isTROK = false;

% Resets the reference population
handles.ReferencePopulation = 1;
set(handles.RefPop_Text,'String','Reference group');

% Also resets the number of subjects variable and associated text
set(handles.Dimensionality_Text, 'String','_ frames x _ voxels (_)');
handles.n_subjects = {};

% Resets the number of datasets entered to 0
handles.n_datasets = 0;

% Makes the reference population list invisible again
set(handles.RefPop,'Visible','off');

% Empties the data, motion, brain information and mask variables
handles.TC = {};
handles.FD = {};
handles.mask = {};
handles.brain_info = {};

%%%%%%%%%% Putting the loading part (bottom) back to normal %%%%%%%%%%%

% We also want to set the TR textbox back to its initial state
handles = TR_Entry_CreateFcn(handles.TR_Entry, eventdata, handles);

% Puts back the seed buttons information to original state
handles.seed = [];
handles.seed2 = [];
set(handles.SeedButton,'BackgroundColor',[0.93 0.84 0.84]);
set(handles.SeedButton,'Enable','off');
set(handles.PlotSeedButton,'Enable','off');

% Removes graph display for the seed
cla(handles.SeedGraphX);
cla(handles.SeedGraphY);
cla(handles.SeedGraphZ);
set(handles.SeedGraphX,'Visible','off');
set(handles.SeedGraphY,'Visible','off');
set(handles.SeedGraphZ,'Visible','off');
set(handles.SliderX,'Visible','off');
set(handles.SliderY,'Visible','off');
set(handles.SliderZ,'Visible','off');
set(handles.XCoordText,'Visible','off');
set(handles.YCoordText,'Visible','off');
set(handles.ZCoordText,'Visible','off');

%%%%%%%%%%%% Putting the seed map part back to normal %%%%%%%%%%%%%%%%%%%

% Resets the variable containing the seed maps of the subjects
handles.SeedMaps = {};
handles.AvgSeedMap = [];

% Not clickable anymore
set(handles.SeedMapPushButton,'Enable','off');

% No more subject list visible
set(handles.SubjectMenu,'Visible','off');

% Resets colorbar display
handles = ResetGraphDisplay(handles.ColorbarSeed,handles);

% Makes the slider and the text linked to slider of the seed map threshold
% back to invisible
set(handles.TSeed_Slider,'Visible','off');
set(handles.TSeed,'Visible','off');

% Resets graphs with seed map plots
handles = ResetGraphDisplay(handles.SeedMapX,handles);
handles = ResetGraphDisplay(handles.SeedMapY,handles);
handles = ResetGraphDisplay(handles.SeedMapZ,handles);
handles = ResetGraphDisplay(handles.SubjSeedMapX,handles);
handles = ResetGraphDisplay(handles.SubjSeedMapY,handles);
handles = ResetGraphDisplay(handles.SubjSeedMapZ,handles);

% Resets associated sliders
set(handles.SeedMapSliderX,'Visible','off');
set(handles.SeedMapSliderY,'Visible','off');
set(handles.SeedMapSliderZ,'Visible','off');

% Resets associated slider texts
set(handles.SeedMap_SliderX,'Visible','off');
set(handles.SeedMap_SliderY,'Visible','off');
set(handles.SeedMap_SliderZ,'Visible','off');

%%%%%%%%%% Putting the CAP analysis part back to normal %%%%%%%%%%%%%%%%%%

% Buttons not clickable anymore (frame selection and clustering)
set(handles.TPSelectionButton,'Enable','off');
set(handles.ClusterButton,'Enable','off');
set(handles.GMM_Button,'Enable','off');
set(handles.AssignButton,'Enable','off');

% Reinitializes motion and the motion box
handles.Tmot = 0.5;
handles = TMotEdit_CreateFcn(handles.TMotEdit,eventdata,handles);

% Reinitializes frame selection threshold and the linked box
handles.T = 0.5;
handles = TEdit_CreateFcn(handles.TEdit,eventdata,handles);

% Resets the frame and percentage retention variables
handles.Xonp = {};
handles.Xonn = {};
handles.RetainedPercentage = {};
handles.FrameIndices = {};

% Resets the violin plot with percentage retained frames
handles = ResetGraphDisplay(handles.TPViolin,handles);

% Resets the parameter input boxes
handles = ClusterEdit_CreateFcn(handles.ClusterEdit,eventdata,handles);
handles = ClusterRepEdit_CreateFcn(handles.ClusterRepEdit,eventdata,handles);
handles = ClusterPpEdit_CreateFcn(handles.ClusterPpEdit,eventdata,handles);
handles = ClusterPnEdit_CreateFcn(handles.ClusterPnEdit,eventdata,handles);
handles = Percentile_Edit_CreateFcn(handles.Percentile_Edit,eventdata,handles);

% Resets the parameters themselves
handles.K = 5;
handles.n_rep = 20;
handles.Pp = 100;
handles.Pn = 100;
handles.percentile = 5;

% Resets the type of performed clustering
handles.CAPType = '';

% Resets the CAP parameters (CAPs, standard deviation within CAPs and
% indices of the CAPs to which all retained frames were assigned)
handles.CAP = [];
handles.STDCAP = [];
handles.idx = {};

% Resets the graph display of the CAP colorbar
handles = ResetGraphDisplay(handles.ColorbarCAP,handles);

% Reset all graph displays for the CAPs
tmpX = {handles.CAP1X,handles.CAP2X,handles.CAP3X,handles.CAP4X,handles.CAP5X,handles.CAP6X};
tmpY = {handles.CAP1Y,handles.CAP2Y,handles.CAP3Y,handles.CAP4Y,handles.CAP5Y,handles.CAP6Y};
tmpZ = {handles.CAP1Z,handles.CAP2Z,handles.CAP3Z,handles.CAP4Z,handles.CAP5Z,handles.CAP6Z};
tmpF = {handles.CAP1_Frames,handles.CAP2_Frames,handles.CAP3_Frames,handles.CAP4_Frames,handles.CAP5_Frames,handles.CAP6_Frames};


for i_CAP = 1:6
    set(tmpF{i_CAP},'Visible','off');
    handles = ResetGraphDisplay(tmpX{i_CAP},handles);
    handles = ResetGraphDisplay(tmpY{i_CAP},handles);
    handles = ResetGraphDisplay(tmpZ{i_CAP},handles);
end

% Resets the sliders and the textboxes for the CAPs
set(handles.CAP_SliderX,'Visible','off');
set(handles.CAP_SliderY,'Visible','off');
set(handles.CAP_SliderZ,'Visible','off');
set(handles.CAP_XC,'Visible','off');
set(handles.CAP_YC,'Visible','off');
set(handles.CAP_ZC,'Visible','off');

% Resets the slider and textbox for the CAPs visualization threshold
set(handles.TVIS_Slider,'Visible','off');
set(handles.TVIS,'Visible','off');

%%%%%%%%%%% Putting the metrics part back to normal %%%%%%%%%%%%%%%%%%%%%%

set(handles.MetricsButton,'Enable','off');

% Resets the metrics variables
handles.TPM = {};
handles.Counts = {};
handles.Number = {};
handles.Avg_Duration = {};
handles.Duration = {};
handles.TM = {};
handles.TPMCum = {};

% Set the sliding lists of subjects invisible again
set(handles.SubjectMenuMetrics,'Visible','off');
set(handles.StateMenu,'Visible','off');

% Resets the colorbars from the metrics part
handles = ResetGraphDisplay(handles.ColorbarSimMat,handles);
handles = ResetGraphDisplay(handles.ColorbarTransMat,handles);

% Resets all the graphs from the metrics part
handles = ResetGraphDisplay(handles.CAP_Mat,handles);
handles = ResetGraphDisplay(handles.TMGraph,handles);
handles = ResetGraphDisplay(handles.TM_Subject,handles);
handles = ResetGraphDisplay(handles.DynStates,handles);
handles = ResetGraphDisplay(handles.CumStates,handles);
handles = ResetGraphDisplay(handles.ViolinCounts,handles);
handles = ResetGraphDisplay(handles.ViolinCountsFrac,handles);
handles = ResetGraphDisplay(handles.ViolinNumber,handles);
handles = ResetGraphDisplay(handles.ViolinDuration,handles);

handles.Log = CAP_AddToLog(handles.Log,'Data cleared');

guidata(hObject, handles); 

% Resets the display of a graph object (called within the clear function)
function handles = ResetGraphDisplay(Graph,handles)
cla(Graph);
set(Graph,'Visible','off');

%% Popup Filling Utilities

% Fills the entries of a pop-up menu with 'Subject _' entries from the
% reference population
function handles = FillSubjectList(ToFill,handles)

tmp_string = {};

for ns = 1:handles.n_subjects{handles.ReferencePopulation}
    tmp_string{ns} = ['Subject ',num2str(ns)];
end

set(ToFill,'String',tmp_string);

clear tmp_string

% Fills the entries of a pop-up menu with the different population entries
function handles = FillPopulationList(ToFill,handles)

tmp_string = {};

for ns = 1:handles.n_datasets
    tmp_string{ns} = [handles.SubjNames{ns}];
end

set(ToFill,'String',tmp_string);

clear tmp_string

% Fills the entries of a pop-up menu with the different population entries
function handles = FillStateList(ToFill,handles)

tmp_string = {};

for ns = 1:handles.K
    tmp_string{ns} = ['State ',num2str(ns)];
end

set(ToFill,'String',tmp_string);

clear tmp_string

%% Project Title Controls
% Selects the title to give to the project

function ProjectTitleText_ButtonDownFcn(hObject, eventdata, handles)

set(hObject,'Enable','on');
set(hObject,'String','');
set(hObject,'FontAngle','normal');
uicontrol(hObject);

function ProjectTitleText_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'Enable','off');
set(hObject,'String','Click to enter...');
set(hObject,'FontAngle','italic');

guidata(hObject, handles);

function ProjectTitleText_Callback(hObject, eventdata, handles)

% If we have entered a valid string, then we name the project as such
if ~isempty((get(hObject,'String')))
    handles.project_title = get(hObject,'String');
    set(hObject,'BackgroundColor', [0.4 0.6 0.4]);
    
    handles.Log = CAP_AddToLog(handles.Log,'Valid project title entered',{handles.project_title},{'New project title'});
    
% If we haven't entered anything, the project is just named 'untitled'
else
    handles.project_title = 'Untitled';
    set(hObject,'BackgroundColor',[0.93 0.84 0.84]);
end

guidata(hObject, handles);

%% Matrix modulation utilities

% Removes NaN-containing lines from a matrix (used for the plotting of
% duration violin plots)
function M2 = DiscardNaN(M)

% We initialize the output matrix as a void one
M2 = [];

% For each row, we count the amount of NaN entries; if not equal to zero,
% then we discard the line
for i = 1:size(M,1)
    if sum(isnan(M(i,:))) > 0
    else
        M2 = [M2;M(i,:)];
    end
end

% Concatenates populations appropriately for Violin plotting
function M2 = ConcatMat(M,n_pop,n_states,n_subjects,type)

% Creates the data matrix (nan values are used to have the same amount of
% data for each group)
M2 = nan(n_pop*n_states,max(cell2mat(n_subjects)));

for i = 1:n_pop
    
    switch type
        case 'Raw counts'
            
            tmp = M{i}.raw.state(:,1:n_states)';
            
            for j = 1:n_states
                M2(i+(j-1)*n_pop,1:size(tmp,2)) = tmp(j,:);
            end
            
            clear tmp
            
        case 'Normalized counts'

            tmp = M{i}.frac.state(:,1:n_states)';

            for j = 1:n_states
                M2(i+(j-1)*n_pop,1:size(tmp,2)) = tmp(j,:);
            end
            
            clear tmp
            
        case 'Number'

            tmp = M{i}(:,4:4+n_states-1)';

            for j = 1:n_states
                M2(i+(j-1)*n_pop,1:size(tmp,2)) = tmp(j,:);
            end
            
            clear tmp
            
        case 'Duration'
            tmp = DiscardNaN(M{i}(:,4:4+n_states-1))';
            
            for j = 1:n_states
                M2(i+(j-1)*n_pop,1:size(tmp,2)) = tmp(j,:);
            end
            
            clear tmp
        case 'FD'
            tmp = M{i};
            
            for j = 1:n_states
                M2(i+(j-1)*n_pop,1:size(tmp,2)) = tmp(j,:);
            end
    end
end
