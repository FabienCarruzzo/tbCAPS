%% This is the main script containing the routines necessary for the use
% of the co-activation pattern analysis (CAP) toolbox
%
% Implemented and written by Thomas Bolton, Medical Image Processing 
% Laboratory (MIP:Lab)
%
% Version 1.0: November 9th 2018: fixing the last remaining issues
%
function varargout = CAP_TB(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CAP_TB_OpeningFcn, ...
                   'gui_OutputFcn',  @CAP_TB_OutputFcn, ...
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





%% Executes when the window opens
function CAP_TB_OpeningFcn(hObject, eventdata, handles, varargin)

%%%%%%%%%%%%%%%%%%%%
% Path and other miscellaneous settings

% Adds the paths to the subfolders of the toolbox that will be important
% for the plotting and the analysis
addpath(genpath('./Plotting'));
addpath(genpath('./Analysis'));
addpath(genpath('./DefaultData'));

% Sets warnings off
warning('off');

% Choose default command line output for CAP_TB
handles.output = hObject;


%%%%%%%%%%%%%%%%%%%%
% Data loading

% TC will contain the time courses of the subjects from the different
% populations (cell array, each cell with size n_TP x n_masked_voxels)
handles.TC = {};

% FD contains the traces of framewise displacement for the subjects (n_TP x
% n_subj per cell of the array, one cell per dataset)
handles.FD = {};

% Information on the NIFTI files from which the data originate
handles.brain_info = {};

% Mask used prior to CAP analysis
handles.mask = {};

% Number of datasets added to the interface. A dataset is defined as a
% population of subjects from the same experimental group (e.g., an 
% ensemble of subjects suffering from the same disorder)
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
% Note: for the moment, this parameter is fixed to 1 (no functionality for
% modifying it yet)
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

% Handles for the TR and whether it is a reasonable value
handles.TR = -inf;
handles.isTROK = false;


%%%%%%%%%%%%%%%%%%%%
% Seed selection and seed maps

% Seed used for the analysis
handles.seed = [];

% Because there can be more than one seed, we create a vector that will
% encompass the multiple information in several colors
handles.seed_display = [];

% Handle to verify whether the amount of seeds has been entered well, how
% many seeds there are, and the type
handles.isSeedOK = false;

% Option to have a subject-specific seed
handles.isSeedSubjectSpecific = 0;

% Number of different seeds
handles.n_seed = 1;

handles.SeedType = 'Average';

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

% Denotes the type of frames (activation, deactivation or both) to use for
% selecting time points
handles.SignMatrix = [1 0; 1 0];

% Activation and deactivation frames kept for all datasets
handles.Xonp = {};
handles.Xonn = {};

% Percentage of frames retained for CAP analysis (discarding both the
% baseline time points and the scrubbed time points)
handles.RetainedPercentage = {};

% Indices of the frames that have been retained (i.e. when do they occur in
% the full time course), of baseline frames, and of scrubbed frames
handles.FrameIndices = {};

handles.idx_sep_seeds = {};

%%%%%%%%%%%%%%%%%%%%
% CAP analysis

% Max number of clusters to verify with consensus clustering
handles.Kmax = 12;

% Percentage of items to use for the consensus clustering folds
handles.PCC = 80;

% Number of times that clustering is run
handles.n_rep = 20;

% Percentage voxels to keep for clustering (positive - Pp - and negative -
% Pn - ones)
handles.Pp = 100;
handles.Pn = 100;

% Number of clusters to use in the analysis
handles.K = 5;

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

% Seed fractions
handles.sfrac = [];

%%%%%%%%%%%%%%%%%%%%
% General utilities

% Log containing the different events summoned from the toolbox
handles.Log = {};

% Colors used in plotting of all populations
handles.PopColor{1} = [255,255,180; 219,224,252; 188,252,188; 230,230,230]/255;
handles.PopColor{2} = [130,48,48; 51,75,163; 59,113,86; 0, 0, 0]/255;

% Project title, by default 'Untitled'
handles.project_title = 'Untitled';

% Directory to which data is to be saved (initially loaded as ./SavedData)
handles.savedir = fullfile(pwd,'SavedData');
set(handles.SaveFolderText,'String',handles.savedir);

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = CAP_TB_OutputFcn(hObject, eventdata, handles) 

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

                % We can now enable the seed selection
                set(handles.SeedButton,'Enable','on');

                % Also, we can now color the button in green
                set(hObject,'BackgroundColor', [101,140,196]/255);

                % If we are loading the first dataset, we convert the underlay
                % to the resolution of the functional data for plotting
                if handles.n_datasets == 1

                    % If we are loading the first dataset, we can get rid
                    % of anything else that we have so far
                    handles = ClearSection3(eventdata,handles);
                    handles = ClearSection4(eventdata,handles);
                    
                    % The brain variable now contains a good resolution
                    % underlay that can directly be overlapped with the
                    % functional data
                    handles.brain = CAP_V2V(handles.brain,handles.Underlay_info.dim,...
                        handles.Underlay_info.mat,handles.brain_info{1}.dim,handles.brain_info{1}.mat);

                elseif handles.n_datasets > 1 && handles.n_datasets < 5 && ~isempty(handles.CAP)
                    
                    % If we are loading one more dataset on top of the
                    % first one, then we rather want to solely clear the
                    % Metrics content; but we may want to keep the rest
                    handles = ClearSection4(eventdata,handles);
                    
                    set(handles.AssignButton,'Enable','on');
    
                    set(handles.CAP_TP,'Visible','on');
                    set(handles.Percentile_Edit,'Visible','on');
                    
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

    % If the TR takes a reasonable value, then we validate it; we enable values
    % between 0.5 and 5s
    if (~isempty(str2double(get(hObject,'String')))) && ...
            (str2double(get(hObject,'String')) > 0.5) && ...
            (str2double(get(hObject,'String')) <= 5)  

        handles.TR = str2double(get(hObject,'String'));    
        set(hObject,'BackgroundColor', [101,140,196]/255);
        handles.isTROK = true;

        handles.Log = CAP_AddToLog(handles.Log,'Correct value of TR entered',{handles.TR},{'TR'});

    % Else, the TR value is not accepted
    else 
        set(hObject,'BackgroundColor', [204,146,146]/255);
        handles.isTROK = false;
    end

guidata(hObject, handles); 



% Executes during creation of the TR textbox
function handles = TR_Entry_CreateFcn(hObject, eventdata, handles)

    set(hObject,'Enable','off');
    set(hObject,'String','Click to enter...');
    set(hObject,'FontAngle','italic');

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','r');
    end

guidata(hObject, handles); 



% Executes when clicking on the TR text space
function TR_Entry_ButtonDownFcn(hObject, eventdata, handles)

    set(hObject,'Enable','on');
    set(hObject,'String','');
    set(hObject,'FontAngle','normal');
    uicontrol(hObject);

guidata(hObject, handles); 





%% Seed Button Controls
% We want to define what happens when loading data (SeedButton) or when
% attempting to plot them (PlotSeedButton)

% --- Executes on button press in CheckS1POS.
function CheckS1POS_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    set(handles.CheckS1NEG,'Value',0);
else
    set(handles.CheckS1POS,'Value',1);
end

handles.SignMatrix(2,:) = [1 0];

rectangle('Position',[0 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);

% Modifies the axis parameters to have a big enough seed display
axis(handles.FancyCircles,'square');

switch handles.n_seed
    case 1
        set(handles.FancyCircles,'xlim',[-10 10]);
        set(handles.FancyCircles,'ylim',[-10 10]);
    case 2
        set(handles.FancyCircles,'xlim',[-10 30]);
        set(handles.FancyCircles,'ylim',[-10 30]);
    case 3
        set(handles.FancyCircles,'xlim',[-10 50]);
        set(handles.FancyCircles,'ylim',[-10 50]);
end

guidata(hObject,handles);



% --- Executes on button press in CheckS1NEG.
function CheckS1NEG_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    handles.SignMatrix(1,:) = [0 1];
    set(handles.CheckS1POS,'Value',0);
    rectangle('Position',[0 0 10 10],'Curvature',[1 1],'FaceColor',[51,75,163]/255,'EdgeColor','none','Parent',handles.FancyCircles);
else
    handles.SignMatrix(1,:) = [1 0];
    set(handles.CheckS1POS,'Value',1);
    rectangle('Position',[0 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);
end

% Modifies the axis parameters to have a big enough seed display
axis(handles.FancyCircles,'square');

switch handles.n_seed
    case 1
        set(handles.FancyCircles,'xlim',[-10 10]);
        set(handles.FancyCircles,'ylim',[-10 10]);
    case 2
        set(handles.FancyCircles,'xlim',[-10 30]);
        set(handles.FancyCircles,'ylim',[-10 30]);
    case 3
        set(handles.FancyCircles,'xlim',[-10 50]);
        set(handles.FancyCircles,'ylim',[-10 50]);
end

guidata(hObject,handles);


% --- Executes on button press in CheckS2POS.
function CheckS2POS_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    set(handles.CheckS2NEG,'Value',0);
else
    set(handles.CheckS2POS,'Value',1);
end

handles.SignMatrix(2,:) = [1 0];
rectangle('Position',[20 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);

% Modifies the axis parameters to have a big enough seed display
axis(handles.FancyCircles,'square');

switch handles.n_seed
    case 1
        set(handles.FancyCircles,'xlim',[-10 10]);
        set(handles.FancyCircles,'ylim',[-10 10]);
    case 2
        set(handles.FancyCircles,'xlim',[-10 30]);
        set(handles.FancyCircles,'ylim',[-10 30]);
    case 3
        set(handles.FancyCircles,'xlim',[-10 50]);
        set(handles.FancyCircles,'ylim',[-10 50]);
end

guidata(hObject,handles);


% --- Executes on button press in CheckS2NEG.
function CheckS2NEG_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    handles.SignMatrix(2,:) = [0 1];
    set(handles.CheckS2POS,'Value',0);
    rectangle('Position',[20 0 10 10],'Curvature',[1 1],'FaceColor',[51,75,163]/255,'EdgeColor','none','Parent',handles.FancyCircles);
else
    handles.SignMatrix(2,:) = [1 0];
    set(handles.CheckS2POS,'Value',1);
    rectangle('Position',[20 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);
end

% Modifies the axis parameters to have a big enough seed display
axis(handles.FancyCircles,'square');

switch handles.n_seed
    case 1
        set(handles.FancyCircles,'xlim',[-10 10]);
        set(handles.FancyCircles,'ylim',[-10 10]);
    case 2
        set(handles.FancyCircles,'xlim',[-10 30]);
        set(handles.FancyCircles,'ylim',[-10 30]);
    case 3
        set(handles.FancyCircles,'xlim',[-10 50]);
        set(handles.FancyCircles,'ylim',[-10 50]);
end

guidata(hObject,handles);



% --- Executes on button press in CheckS3POS.
function CheckS3POS_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    set(handles.CheckS3NEG,'Value',0);
else
    set(handles.CheckS3POS,'Value',1);
end

handles.SignMatrix(3,:) = [1 0];

rectangle('Position',[40 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);

% Modifies the axis parameters to have a big enough seed display
axis(handles.FancyCircles,'square');

switch handles.n_seed
    case 1
        set(handles.FancyCircles,'xlim',[-10 10]);
        set(handles.FancyCircles,'ylim',[-10 10]);
    case 2
        set(handles.FancyCircles,'xlim',[-10 30]);
        set(handles.FancyCircles,'ylim',[-10 30]);
    case 3
        set(handles.FancyCircles,'xlim',[-10 50]);
        set(handles.FancyCircles,'ylim',[-10 50]);
end

guidata(hObject,handles);

% --- Executes on button press in CheckS3NEG.
function CheckS3NEG_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    handles.SignMatrix(3,:) = [0 1];
    set(handles.CheckS3POS,'Value',0);
    rectangle('Position',[40 0 10 10],'Curvature',[1 1],'FaceColor',[51,75,163]/255,'EdgeColor','none','Parent',handles.FancyCircles);
else
    handles.SignMatrix(3,:) = [1 0];
    set(handles.CheckS3POS,'Value',1);
    rectangle('Position',[40 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);
end

% Modifies the axis parameters to have a big enough seed display
axis(handles.FancyCircles,'square');

switch handles.n_seed
    case 1
        set(handles.FancyCircles,'xlim',[-10 10]);
        set(handles.FancyCircles,'ylim',[-10 10]);
    case 2
        set(handles.FancyCircles,'xlim',[-10 30]);
        set(handles.FancyCircles,'ylim',[-10 30]);
    case 3
        set(handles.FancyCircles,'xlim',[-10 50]);
        set(handles.FancyCircles,'ylim',[-10 50]);
end

guidata(hObject,handles);



% Executes when the user changes the type of seed relationship desired
% (Average, Intersection, Union of seed signals)
function SeedPopup_Callback(hObject, eventdata, handles)

% Clears and makes the circles display reflecting seed types unvisible
handles = ResetGraphDisplay(handles.FancyCircles,handles);

% We consider the value selected by the user for the first seed between
% '+' or '-', and update the FancyCircles plot accordingly
if get(handles.CheckS1POS,'Value')
    rectangle('Position',[0 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);
elseif get(handles.CheckS1NEG,'Value')
    rectangle('Position',[0 0 10 10],'Curvature',[1 1],'FaceColor',[51,75,163]/255,'EdgeColor','none','Parent',handles.FancyCircles);
end

% The same is done for the other seeds, if more than one has been loaded
% (or if more than two have been loaded)
if handles.n_seed > 1

    if get(handles.CheckS2POS,'Value')
        rectangle('Position',[20 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);
    elseif get(handles.CheckS2NEG,'Value')
        rectangle('Position',[20 0 10 10],'Curvature',[1 1],'FaceColor',[51,75,163]/255,'EdgeColor','none','Parent',handles.FancyCircles);
    end
    
    if handles.n_seed > 2
        
        if get(handles.CheckS3POS,'Value')
            rectangle('Position',[40 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);
        elseif get(handles.CheckS3NEG,'Value')
            rectangle('Position',[40 0 10 10],'Curvature',[1 1],'FaceColor',[51,75,163]/255,'EdgeColor','none','Parent',handles.FancyCircles);
        end
    end
end

% The second step in the plotting is to add the "arrows" that link circles
% differently depending on whether we consider a Union or an Intersection
switch get(hObject,'Value')
    
    % If we have the "average" pick
    case 1
        % If the user selects "average" but the number of seeds entered
        % exceeds 1, then by default we revert back to the "union" case,
        % since the user should enter one averaged seed file to have an
        % averaged seed output
        if handles.n_seed > 1
            set(hObject,'Value',2);
            handles.SeedType = 'Union';
            rectangle('Position',[12 -10 6 8],'Curvature',[0.8 0.8],'Parent',handles.FancyCircles);
            rectangle('Position',[11 -5 8 4],'Curvature',[0.8 0.8],'EdgeColor','none','FaceColor','w','Parent',handles.FancyCircles);
            if handles.n_seed > 2
                rectangle('Position',[32 -10 6 8],'Curvature',[0.8 0.8],'Parent',handles.FancyCircles);
                rectangle('Position',[31 -5 8 4],'Curvature',[0.8 0.8],'EdgeColor','none','FaceColor','w','Parent',handles.FancyCircles);
            end
        else
            handles.SeedType = 'Average';
        end
        
    % If we have the union pick
    case 2
        handles.SeedType = 'Union';
        rectangle('Position',[12 -10 6 8],'Curvature',[0.8 0.8],'Parent',handles.FancyCircles);
        rectangle('Position',[11 -5 8 4],'Curvature',[0.8 0.8],'EdgeColor','none','FaceColor','w','Parent',handles.FancyCircles);
        if handles.n_seed > 2
            rectangle('Position',[32 -10 6 8],'Curvature',[0.8 0.8],'Parent',handles.FancyCircles);
            rectangle('Position',[31 -5 8 4],'Curvature',[0.8 0.8],'EdgeColor','none','FaceColor','w','Parent',handles.FancyCircles);
        end
        
    % If we have the intersection pick
    case 3
        handles.SeedType = 'Intersection';
        rectangle('Position',[12 10 6 8],'Curvature',[0.8 0.8],'Parent',handles.FancyCircles);
        rectangle('Position',[11 9 8 4],'Curvature',[0.8 0.8],'EdgeColor','none','FaceColor','w','Parent',handles.FancyCircles);
        if handles.n_seed > 2
            rectangle('Position',[32 10 6 8],'Curvature',[0.8 0.8],'Parent',handles.FancyCircles);
            rectangle('Position',[31 9 8 4],'Curvature',[0.8 0.8],'EdgeColor','none','FaceColor','w','Parent',handles.FancyCircles);
        end
end

% Modifies the axis parameters to have a big enough seed display
axis(handles.FancyCircles,'square');

switch handles.n_seed
    case 1
        set(handles.FancyCircles,'xlim',[-10 10]);
        set(handles.FancyCircles,'ylim',[-10 10]);
    case 2
        set(handles.FancyCircles,'xlim',[-10 30]);
        set(handles.FancyCircles,'ylim',[-10 30]);
    case 3
        set(handles.FancyCircles,'xlim',[-10 50]);
        set(handles.FancyCircles,'ylim',[-10 50]);
end

handles.Log = CAP_AddToLog(handles.Log,'Seed union status changed',{handles.SeedType},{'Status'});

guidata(hObject, handles);



% Executes during object creation, after setting all properties.
function SeedPopup_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

guidata(hObject, handles);



% Executes when clicking on 'B. Select a seed'
function SeedButton_Callback(hObject, eventdata, handles)

% We want to clear everything from and after this stage, since we will
% compute results with a new seed
handles = ClearSection3(eventdata,handles);
handles = ClearSection4(eventdata,handles);

% Multiselection is on, so that as many as three files can be picked
[filename_seed,pathname_seed]=uigetfile({'*.*','All Files'},...
  'Select Seed File...','MultiSelect','on');

% Stores the type of frames to retain for each seed
handles.SignMatrix = [1 0; 1 0; 1 0];

% Resets the graph displays as well as the seed-related parameters
handles = ResetGraphDisplay(handles.FancyCircles,handles);

set(handles.Seed1Text,'Visible','off');
set(handles.Seed2Text,'Visible','off');
set(handles.Seed3Text,'Visible','off');

set(handles.SeedPlusText,'Visible','off');
set(handles.SeedMinusText,'Visible','off');

set(handles.CheckS1POS,'Visible','off');
set(handles.CheckS1NEG,'Visible','off');

set(handles.CheckS1POS,'Value',1);
set(handles.CheckS1NEG,'Value',0);

set(handles.CheckS2POS,'Visible','off');
set(handles.CheckS2NEG,'Visible','off');

set(handles.CheckS2POS,'Value',1);
set(handles.CheckS2NEG,'Value',0);

set(handles.CheckS3POS,'Visible','off');
set(handles.CheckS3NEG,'Visible','off');

set(handles.CheckS3POS,'Value',1);
set(handles.CheckS3NEG,'Value',0);

handles = ResetGraphDisplay(handles.SeedGraphX,handles);
handles = ResetGraphDisplay(handles.SeedGraphZ,handles);
handles = ResetGraphDisplay(handles.SeedMapX,handles);
handles = ResetGraphDisplay(handles.SeedMapZ,handles);

% If the user has indeed entered files
if ~isequal(filename_seed,0) || ~isequal(pathname_seed,0)
    % There should be three or less selected files. In this switch, we test
    % for the amount of entered files
   
    % In the case in which only one seed file is entered ('char' type),
    % we convert into an array
    if strcmp(class(filename_seed),'char')
        filename_seed = {filename_seed};
    end
    
    % If we enter that statement, it means that we have only one seed type
    % across subjects (that is, no subject-specific data)
    if length(filename_seed) <= 3

        % Number of entered seeds
        handles.n_seed = length(filename_seed);

        for myindex = 1:length(filename_seed)
            File_seed = fullfile(pathname_seed, filename_seed{myindex});
            tmp = load(File_seed);
            assignin('base','tmp', tmp);
            tmp = struct2array(tmp);

            % If the file is of suitable dimensions
            if islogical(tmp) && size(tmp,2) == 1 && size(tmp,1) == sum(handles.mask{1})

                % Then we put it in the handles, enable the plotting button, and
                % make the seed selection button green
                handles.seed(:,myindex) = tmp;
                
                handles.Log = CAP_AddToLog(handles.Log,'Seed chosen',{File_seed},{'Seed file'});
            else
                errordlg('The file you entered appears to be of wrong dimensions...');
                handles = ClearDataButton_Callback(handles.ClearDataButton, eventdata, handles);
            end

        end
                
        % If we survive all the above, we can see that the seed files are
        % good
        handles.isSeedOK = true;

        % We also enable to plot the seed results
        set(handles.PlotSeedButton,'Enable','on');
        set(handles.PlotSeedButton,'Visible','on');
        
        % We make the text legends visible
        set(handles.S_SEED1,'Visible','on');
        
        set(handles.SeedButton,'BackgroundColor', [101,140,196]/255);

        % We can now go through the next parts of the analysis, so we
        % enable the related buttons
        set(handles.TPSelectionButton,'Enable','on');
        set(handles.SeedMapPushButton,'Enable','on');
        
        set(handles.TPSelectionButton,'Visible','on');
        set(handles.SeedMapPushButton,'Visible','on');
        
        % Makes other TP selection utilities visible
        set(handles.PRadio,'Visible','on');
        set(handles.TRadio,'Visible','on');
        set(handles.uibuttongroup7,'Visible','on');
        set(handles.TText,'Visible','on');
        set(handles.TMotText,'Visible','on');
        set(handles.TEdit,'Visible','on');
        set(handles.TMotEdit,'Visible','on');
        
        % We also see the displays for entering seed specifications
        set(handles.SeedPopup,'Visible','on');
        
        handles.Log = CAP_AddToLog(handles.Log,'Correct amount of seeds entered',{handles.n_seed},{'Seed amount'});
    
        handles.seed_display = zeros(length(handles.seed(:,1)),1);
        
        set(handles.CheckS1POS,'Visible','on');
        set(handles.CheckS1NEG,'Visible','on');
        CheckS1POS_Callback(handles.CheckS1POS,eventdata,handles);
        
        set(handles.Seed1Text,'Visible','on');
        
        set(handles.SeedPlusText,'Visible','on');
        set(handles.SeedMinusText,'Visible','on');
        
        % If there are more than one seed, then we allow the popup button to be
        % changed for a more complex seed use choice
        if handles.n_seed > 1
            set(handles.SeedPopup,'Enable','on');
            set(handles.SeedPopup,'Value',2);
            handles.SeedType = 'Union';
            
            % We also update the circles in the seed illustration
            rectangle('Position',[12 -10 6 8],'Curvature',[0.8 0.8],'Parent',handles.FancyCircles);
            rectangle('Position',[11 -5 8 4],'Curvature',[0.8 0.8],'EdgeColor','none','FaceColor','w','Parent',handles.FancyCircles);
            
            % We fill seed_display with one scalar value across seed voxels
            % per seed (to have different colors plotted in the seed choice
            % graph
            useless_vector = [0.25,0.75,1];
            
            for se = 1:handles.n_seed
                handles.seed_display = handles.seed_display + useless_vector(se)*handles.seed(:,se);
            end
            
            set(handles.S_SEED2,'Visible','on');
            
            set(handles.CheckS2POS,'Visible','on');
            set(handles.CheckS2NEG,'Visible','on');
            CheckS2POS_Callback(handles.CheckS2POS,eventdata,handles);
            
            set(handles.Seed2Text,'Visible','on');
        
            % Same for 3 seeds
            if handles.n_seed > 2
                set(handles.CheckS3POS,'Visible','on');
                set(handles.CheckS3NEG,'Visible','on');
                CheckS3POS_Callback(handles.CheckS3POS,eventdata,handles);
                
                set(handles.Seed3Text,'Visible','on');
                
                rectangle('Position',[32 -10 6 8],'Curvature',[0.8 0.8],'Parent',handles.FancyCircles);
                rectangle('Position',[31 -5 8 4],'Curvature',[0.8 0.8],'EdgeColor','none','FaceColor','w','Parent',handles.FancyCircles);
                
                set(handles.S_SEED3,'Visible','on');
            end
         
        % Entered if we have just one seed
        else
            set(handles.SeedPopup,'Enable','off');
            set(handles.SeedPopup,'Value',1);
            handles.SeedType = 'Average';
            
            handles.seed_display = handles.seed;
            
            set(handles.CheckS1POS,'Visible','on');
            set(handles.CheckS1NEG,'Visible','on');
        end
        
        rectangle('Position',[0 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);
        
        if handles.n_seed > 1
            rectangle('Position',[20 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);
            
            if handles.n_seed > 2
                rectangle('Position',[40 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);
            end
        end

    % If the seed file is of length n_subjects, then we want a
    % subject-specific seed scheme to run
    elseif length(filename_seed) == handles.n_subjects{1}
        
        handles.isSeedSubjectSpecific = 1;
        
        % Then, we only allow one seed to be considered
        handles.n_seed = 1;

        % For each entry (i.e., each subject data), we load the seed
        for myindex = 1:length(filename_seed)
            File_seed = fullfile(pathname_seed, filename_seed{myindex});
            tmp = load(File_seed);
            assignin('base','tmp', tmp);
            tmp = struct2array(tmp);

            % If the file is of suitable dimensions
            if islogical(tmp) && size(tmp,2) == 1 && size(tmp,1) == sum(handles.mask{1})

                % Then we put it in the handles, enable the plotting button, and
                % make the seed selection button green
                handles.seed(:,myindex) = tmp;
                
                handles.Log = CAP_AddToLog(handles.Log,'Seed chosen',{File_seed},{'Seed file'});
            else
                errordlg('The file you entered appears to be of wrong dimensions...');
                handles = ClearDataButton_Callback(handles.ClearDataButton, eventdata, handles);
            end

        end
        
        handles.seed_display = zeros(length(handles.seed(:,1)),1);
                
        % If we survive all the above, we can see that the seed files are
        % good
        handles.isSeedOK = true;

        % We also enable to plot the seed results
        set(handles.PlotSeedButton,'Enable','on');
        set(handles.PlotSeedButton,'Visible','on');
        
        set(handles.S_SEED1,'Visible','on');
        
        set(handles.SeedButton,'BackgroundColor', [101,140,196]/255);

        set(handles.SeedPopup,'Enable','off');
        set(handles.SeedPopup,'Value',1);
        handles.SeedType = 'Average';
        set(handles.CheckS1POS,'Visible','on');
        set(handles.CheckS1NEG,'Visible','on');
        
        % Creates the seed information to plot in the subject-specific case
        for idx_seed = 1:size(handles.seed,2)
            
            handles.seed_display = handles.seed_display + handles.seed(:,idx_seed);
        end
        
        handles.seed_display = handles.seed_display/handles.n_subjects{1};
        
        % We can now go through the next parts of the analysis, so we
        % enable the related buttons
        set(handles.TPSelectionButton,'Enable','on');
        set(handles.SeedMapPushButton,'Enable','on');
        
        set(handles.TPSelectionButton,'Visible','on');
        set(handles.SeedMapPushButton,'Visible','on');
        
        set(handles.PRadio,'Visible','on');
        set(handles.TRadio,'Visible','on');
        set(handles.uibuttongroup7,'Visible','on');
        set(handles.TText,'Visible','on');
        set(handles.TMotText,'Visible','on');
        set(handles.TEdit,'Visible','on');
        set(handles.TMotEdit,'Visible','on');
        
        % We also see the displays for entering seed specifications
        set(handles.SeedPopup,'Visible','on');
        
        handles.Log = CAP_AddToLog(handles.Log,'Correct subject-specific seed data entered',{handles.n_seed},{'Seed amount'});
    
        handles.seed_display = zeros(length(handles.seed(:,1)),1);
        
        set(handles.CheckS1POS,'Visible','on');
        set(handles.CheckS1NEG,'Visible','on');
        CheckS1POS_Callback(handles.CheckS1POS,eventdata,handles);
        
        set(handles.Seed1Text,'Visible','on');
        
        set(handles.SeedPlusText,'Visible','on');
        set(handles.SeedMinusText,'Visible','on');
        
        rectangle('Position',[0 0 10 10],'Curvature',[1 1],'FaceColor',[150,48,48]/255,'EdgeColor','none','Parent',handles.FancyCircles);
 
        set(handles.FancyCircles,'xlim',[-10 10]);
        set(handles.FancyCircles,'ylim',[-10 10]);

    else
        errordlg('Problem with the amount of seed files entered !');
        handles = ClearDataButton_Callback(handles.ClearDataButton, eventdata, handles);
    end
    
    % Updates the limits of the plot
    switch handles.n_seed
        case 1
            set(handles.FancyCircles,'xlim',[-10 10]);
            set(handles.FancyCircles,'ylim',[-10 10]);
        case 2
            set(handles.FancyCircles,'xlim',[-10 30]);
            set(handles.FancyCircles,'ylim',[-10 30]);
        case 3
            set(handles.FancyCircles,'xlim',[-10 50]);
            set(handles.FancyCircles,'ylim',[-10 50]);
    end
    
else
    errordlg('You did not enter a seed file !');
    handles = ClearDataButton_Callback(handles.ClearDataButton, eventdata, handles);
end
        
guidata(hObject, handles);



% Executes when clicking on 'Plot Seed'
function PlotSeedButton_Callback(hObject, eventdata, handles)

% Clears the present graph content
cla(handles.SeedGraphX);
cla(handles.SeedGraphZ);

% Plots the slices within the graph windows
handles.SeedGraphX = plot_slice(handles.seed_display,get(handles.TVIS_Slider,...
    'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,...
    handles.brain_info{handles.ReferencePopulation},'X',get(handles.SliderX,...
    'Value'),handles.SeedGraphX);

handles.SeedGraphZ = plot_slice(handles.seed_display,get(handles.TVIS_Slider,...
    'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,...
    handles.brain_info{handles.ReferencePopulation},'Z',get(handles.SliderZ,...
    'Value'),handles.SeedGraphZ);

% Sets the sliders to visible
set(handles.SliderX,'Visible','on');
set(handles.SliderZ,'Visible','on');

% Sets the text values at the ones of the sliders
set(handles.XCoordText,'String',['X: ',sprintf('%.2f',get(handles.SliderX,'Value'))]);
set(handles.ZCoordText,'String',['Z: ',sprintf('%.2f',get(handles.SliderZ,'Value'))]);

% Sets the visibility of the slider texts to on
set(handles.XCoordText,'Visible','on');
set(handles.ZCoordText,'Visible','on');

handles.Log = CAP_AddToLog(handles.Log,'Seed plots activated');

guidata(hObject, handles);





%% Seed sliders interactions
% For the below functions, the goal is to change the value of the slider
% textboxes when the sliders are moved, and to update the graph display
% accordingly. For this purpose, cla is used to clear graph content prior
% to a new display

% Executes on slider movement.
function SliderX_Callback(hObject, eventdata, handles)

% Clears the content of the graph
cla(handles.SeedGraphX);

% Gets the MNI slice coordinate value associated to the new display
set(handles.XCoordText,'String',['X: ',sprintf('%.2f',get(hObject,'Value'))]);

% Slice plotting itself; 1.5 is the color past which the display saturates
handles.SeedGraphX = plot_slice(handles.seed_display,get(handles.TVIS_Slider,...
    'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,...
    handles.brain_info{handles.ReferencePopulation},'X',get(hObject,'Value'),...
    handles.SeedGraphX);

guidata(hObject, handles);



% Executes during object creation, after setting all properties.
function SliderX_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

guidata(hObject, handles); 



% Executes on slider movement.
function SliderZ_Callback(hObject, eventdata, handles)

cla(handles.SeedGraphZ);
set(handles.ZCoordText,'String',['Z: ',sprintf('%.2f',get(hObject,'Value'))]);
handles.SeedGraphZ = plot_slice(handles.seed_display,get(handles.TVIS_Slider,...
    'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,...
    handles.brain_info{handles.ReferencePopulation},'Z',get(hObject,'Value'),...
    handles.SeedGraphZ);

guidata(hObject, handles); 



% Executes during object creation, after setting all properties.
function SliderZ_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

guidata(hObject, handles); 





%% Seed Map Computation Button
% When pressing on this button, classical seed maps are computed for the
% population of subjects chosen as the reference one in the loading part.
% The last entered seed is used

function SeedMapPushButton_Callback(hObject, eventdata, handles)

% Computes seed maps for each subject and for the population, using the
% data from the chosen reference population
[~,handles.AvgSeedMap] = CAP_Compute_SeedMap(handles.TC{handles.ReferencePopulation},handles.seed,handles.isSeedSubjectSpecific);

% Graphical displays

% Making the plots, texts and sliders visible
set(handles.SeedMapX,'Visible','on');
set(handles.SeedMapZ,'Visible','on');

set(handles.SeedMap_SliderX,'Visible','on');
set(handles.SeedMap_SliderZ,'Visible','on');
set(handles.SeedMapSliderX,'Visible','on');
set(handles.SeedMapSliderZ,'Visible','on');

set(handles.TSeed_Slider,'Visible','on');
set(handles.TSeed,'Visible','on');
set(handles.ColorbarSeed,'Visible','on');

% Writing down the text with current MNI coordinates
set(handles.SeedMap_SliderX,'String',['X: ',sprintf('%.2f',get(handles.SeedMapSliderX,'Value'))]);
set(handles.SeedMap_SliderZ,'String',['Z: ',sprintf('%.2f',get(handles.SeedMapSliderZ,'Value'))]);

% Clears previous plot contents (in case we want to re-plot after changing
% the seed)
cla(handles.SeedMapX);
cla(handles.SeedMapZ);

% Plots new slices
handles.SeedMapX = plot_slice(handles.AvgSeedMap,0.25,1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'X',get(handles.SeedMapSliderX,'Value'),handles.SeedMapX);
handles.SeedMapZ = plot_slice(handles.AvgSeedMap,0.25,1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Z',get(handles.SeedMapSliderZ,'Value'),handles.SeedMapZ);

% Adds the colorbar for the seed maps (between -1 and 1)
handles.ColorbarSeed = Create_CAP_colorbar(-1,1,0.5,get(handles.TSeed_Slider,'Value'),'',handles.ColorbarSeed,'Horizontal','div','RdBu',1000);

handles.Log = CAP_AddToLog(handles.Log,'Seed maps displayed');

guidata(hObject,handles);





%% Slider Controls (MNI coordinates)
% We want to reload the seed images with the new parameters when changing a
% slider, so we clear the previous display, change the text summarizing the
% MNI coordinate where we stand, and plot the new image

function SeedMapSliderX_Callback(hObject, eventdata, handles)

% Clears graphs
cla(handles.SeedMapX);

% Changes slider texts
set(handles.SeedMap_SliderX,'String',['X: ',sprintf('%.2f',get(hObject,'Value'))]);

% Plots new slices
handles.SeedMapX = plot_slice(handles.AvgSeedMap,get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'X',get(hObject,'Value'),handles.SeedMapX);

guidata(hObject, handles);



function SeedMapSliderX_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function SeedMapSliderZ_Callback(hObject, eventdata, handles)

cla(handles.SeedMapZ);

% Changes slider texts
set(handles.SeedMap_SliderZ,'String',['Z: ',sprintf('%.2f',get(hObject,'Value'))]);

% Plots new slices
handles.SeedMapZ = plot_slice(handles.AvgSeedMap,get(handles.TSeed_Slider,'Value'),1,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Z',get(hObject,'Value'),handles.SeedMapZ);

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
cla(handles.SeedMapZ);

% Plots new slices (average seed maps)
handles.SeedMapX = plot_slice(handles.AvgSeedMap,get(hObject,'Value'),1,...
    handles.mask{handles.ReferencePopulation},handles.brain,...
    handles.brain_info{handles.ReferencePopulation},'X',...
    get(handles.SeedMapSliderX,'Value'),handles.SeedMapX);

handles.SeedMapZ = plot_slice(handles.AvgSeedMap,get(hObject,'Value'),1,...
    handles.mask{handles.ReferencePopulation},handles.brain,...
    handles.brain_info{handles.ReferencePopulation},'Z',...
    get(handles.SeedMapSliderZ,'Value'),handles.SeedMapZ);

% Modifies the text
set(handles.TSeed,'String',['Tv: ',sprintf('%.2f',get(hObject,'Value'))]);

% Clears and replots the colorbar
cla(handles.ColorbarSeed);
handles.ColorbarSeed = Create_CAP_colorbar(-1,1,0.5,get(hObject,'Value'),'',handles.ColorbarSeed,'Horizontal','div','RdBu',1000);

guidata(hObject,handles);



function TSeed_Slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end





%% Motion parameter entry
% In the following, we have the functions controling the motion threshold
% (Power's FD) value

function TMotEdit_Callback(hObject, eventdata, handles)

% If we enter a reasonable value, it is taken as a new threshold
if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 0) && (str2double(get(hObject,'String')) <= 0.5)
    handles.Tmot = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [101,140,196]/255);
    
    handles.Log = CAP_AddToLog(handles.Log,'Valid motion threshold value entered',{handles.Tmot},{'Motion threshold value'});
    
% If we set something wrong again, we set the threshold value back to the
% default of 0.5
else
    set(hObject,'BackgroundColor', [203,146,146]/255);
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
% What comes below describes the control of the frame selection threshold
% value entered by the user

% Threshold for selection of the frames to keep
function TEdit_Callback(hObject, eventdata, handles)

% If we enter a reasonable value, it is taken as the new threshold
if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 0) && (str2double(get(hObject,'String')) <= 100)
    handles.T = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [101,140,196]/255);
    
    handles.Log = CAP_AddToLog(handles.Log,'Valid (de)activation threshold entered',{handles.T},{'Threshold value'});
    
% If we set something wrong again, we set the threshold value back to the
% default of 0.5
else
    set(hObject,'BackgroundColor',[203,146,146]/255);
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
% The below code is run when the user wishes to change between frame
% selection modes

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
    [handles.Xonp{n_ds},p,Indices,...
        handles.idx_sep_seeds{n_ds}] = CAP_find_activity(handles.TC{n_ds},...
        handles.seed,handles.T,handles.FD{n_ds},handles.Tmot,...
        handles.SelMode,handles.SeedType,handles.SignMatrix,handles.isSeedSubjectSpecific);
    
    % Percentage of retained frames across subjects
    handles.RetainedPercentage{n_ds} = p(3,:);

    % Indices of the frames that have been retained (used later for metrics
    % computations)
    handles.FrameIndices{n_ds} = Indices; 
end

% Enables to go to the next step of the analysis and cluster the extracted
% frames
set(handles.ClusterButton,'Enable','on');
set(handles.CCButton,'Enable','on');

set(handles.ClusterButton,'Visible','on');
set(handles.CCButton,'Visible','on');

% Sets related displays to visible
set(handles.CAP_Kmax,'Visible','on');
set(handles.CAP_PPC,'Visible','on');
set(handles.CAP_N,'Visible','on');
set(handles.CAP_K,'Visible','on');
set(handles.CAP_PP,'Visible','on');
set(handles.CAP_PN,'Visible','on');

set(handles.KRange_Edit,'Visible','on');
set(handles.PCC_Edit,'Visible','on');
set(handles.ClusterEdit,'Visible','on');
set(handles.ClusterRepEdit,'Visible','on');
set(handles.ClusterPpEdit,'Visible','on');
set(handles.ClusterPnEdit,'Visible','on');

tmp_toplot = ConcatMat(handles.RetainedPercentage,handles.n_datasets,1,handles.n_subjects,'FD');

% Displays the violin plot of subject scrubbing percentage for the
% reference population
[~,~,handles.TPViolin] = MakeViolin(tmp_toplot,handles.TPViolin,{' '},'Frames ret. [%]',handles.PopColor,handles.n_datasets,1);
set(handles.TPViolin,'Visible','on');

clear tmp_toplot

handles.Log = CAP_AddToLog(handles.Log,'Time points selected',{['1 to ',num2str(handles.n_datasets)],handles.SelMode,handles.T,handles.Tmot},{'Datasets indices','Selection mode','Activation threshold','Motion threshold'});

guidata(hObject, handles);





%% Data saving and loading
% The functions below are summoned when the user wishes to save his/her
% data into a MATLAB structure, or to load previously processed data (until
% the end of spatio-temporal selection) and attempt clustering computations
% and CAP generation

% Save folder change function
function SaveFolderButton_Callback(hObject, eventdata, handles)

% Selection of a directory
[dirname]=uigetdir('*.*','Please select a save directory');
handles.savedir = dirname;

% If the user has indeed chosen a directory, we set it as the new save
% folder
if ~isequal(dirname,0)
    set(handles.SaveFolderText,'String',handles.savedir);
    set(hObject,'BackgroundColor', [101,140,196]/255);
    
    handles.Log = CAP_AddToLog(handles.Log,'Save folder changed',...
    {handles.savedir},...
    {'Save folder'});
else
    errordlg('Please select a directory !');
end

guidata(hObject,handles);

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

% Motion threshold and activation threshold used in time points selection
SAVED.TPSelData.Tmot = handles.Tmot;
SAVED.TPSelData.T = handles.T;

% Type of frame selection used
SAVED.TPSelData.SelMode = handles.SelMode;

% Frames that were considered in the clustering process
% SAVED.TPSelData.Act = handles.Xonp;

% Type of seed computation chosen ('Avg','Union','Intersection')
SAVED.TPSelData.SeedType = handles.SeedType;

% Indices, for each seed used, of the retained frames across subjects
SAVED.TPSelData.idx_sep_seeds = handles.idx_sep_seeds;

% Percentage frames retained for the clustering
SAVED.TPSelData.PercRetained = handles.RetainedPercentage;

% Computed seed maps (average and subject-wise)
SAVED.SeedMap.AvgMap = handles.AvgSeedMap;

% Parameters used for clustering
SAVED.ClusterData.N = handles.n_rep;
SAVED.ClusterData.K = handles.K;
SAVED.ClusterData.Pp = handles.Pp;
SAVED.ClusterData.Pn = handles.Pn;

% CAP data
SAVED.ClusterData.CAPs = handles.CAP;
SAVED.ClusterData.StdCAPs = handles.STDCAP;
SAVED.ClusterData.idx = handles.idx;

% Computed metrics
SAVED.Metrics.TPM = handles.TPM;
SAVED.Metrics.Counts = handles.Counts;
SAVED.Metrics.Number = handles.Number;
SAVED.Metrics.Avg_Duration = handles.Avg_Duration;
SAVED.Metrics.Duration = handles.Duration;
SAVED.Metrics.TM = handles.TM;
SAVED.Metrics.SeedFrac = handles.sfrac;

%[tmp_date,tmp_date2] = strtok(date,'-');
%[tmp_date2,tmp_date3] = strtok(tmp_date2,'-');
%tmp_date3 = strtok(tmp_date3,'-');

% Name that will be given to the saved files
% fancy_name = [handles.project_title,'_',tmp_date,'_',tmp_date2,'_',tmp_date3,'_',...
%     num2str(hour(now)),'_',num2str(minute(now)),'_',...
%     num2str(round(second(now)))];

fancy_name = [handles.project_title];

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

guidata(hObject,handles);



% This function is summoned when the user wishes to load previously
% computed data and pursue an analysis
function LoadButton_Callback(hObject, eventdata, handles)

% We want to verify 
[filename1,pathname1]=uigetfile({'*.*','All Files'},...
      'Select struct file to load...','MultiSelect','on');

    % If the user has indeed entered a file
    if ~isequal(filename1,0) || ~isequal(pathname1,0)
        
        Loaded_File = fullfile(pathname1, filename1);
        tmp = load(Loaded_File);
        assignin('base','tmp', tmp);
        tmp = struct2array(tmp);
    end
    
    % Checks that the right fields do exist
    try 
        if tmp.SubjData.n_subjects{1} > 3 &&...
        ~isempty(tmp.BrainData.brain_info) &&...
        sum(tmp.BrainData.mask{1}) == size(tmp.BrainData.seed,1) && ...
        length(tmp.TPSelData.Act{1}) == tmp.SubjData.n_subjects{1} &&...
        size(tmp.TPSelData.Act{1}{1},1) == sum(tmp.BrainData.mask{1})

            handles.brain_info = tmp.BrainData.brain_info;
            handles.Xonp = tmp.TPSelData.Act{1};
            handles.mask = tmp.BrainData.mask;
            handles.SeedType = tmp.TPSelData.SeedType;
            handles.idx_sep_seeds = tmp.TPSelData.idx_sep_seeds;
            
            % Makes the LOAD button green to indicate that loading worked
            set(handles.LoadButton,'BackgroundColor',[59 113 86]/255);
            
            % Enables to go to the next step of the analysis and cluster the extracted
            % frames
            set(handles.ClusterButton,'Enable','on');
            set(handles.CCButton,'Enable','on');

            set(handles.ClusterButton,'Visible','on');
            set(handles.CCButton,'Visible','on');

            % Sets related displays to visible
            set(handles.CAP_Kmax,'Visible','on');
            set(handles.CAP_PPC,'Visible','on');
            set(handles.CAP_N,'Visible','on');
            set(handles.CAP_K,'Visible','on');
            set(handles.CAP_PP,'Visible','on');
            set(handles.CAP_PN,'Visible','on');

            set(handles.KRange_Edit,'Visible','on');
            set(handles.PCC_Edit,'Visible','on');
            set(handles.ClusterEdit,'Visible','on');
            set(handles.ClusterRepEdit,'Visible','on');
            set(handles.ClusterPpEdit,'Visible','on');
            set(handles.ClusterPnEdit,'Visible','on');
        end
    catch
        errordlg('Problem in loading the data...');
    end    
                        
guidata(hObject,handles);



% Executes when pressing on the 'CLEAR' button for data loading; supposed
% to set everything back to normal (when the window opened)
function handles = ClearDataButton_Callback(hObject, eventdata, handles)
     
    handles = ClearSection1(eventdata,handles);
    handles = ClearSection2(eventdata,handles);
    handles = ClearSection3(eventdata,handles);
    handles = ClearSection4(eventdata,handles);
    
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

    handles.Log = CAP_AddToLog(handles.Log,'Data cleared');

guidata(hObject, handles); 



% Clears the content of section 1 only
function handles = ClearSection1(eventdata, handles)

% Makes 'A. Load data' red again
set(handles.DataButton,'BackgroundColor',[204,146,146]/255);

% Same for Save folder button
set(handles.SaveFolderButton,'BackgroundColor',[204,146,146]/255);

% Resets the time point and voxel parameters
handles.SubjSize.TP = -inf;
handles.SubjSize.VOX = -inf;

% Resets the TR
handles.TR = -inf;
handles.isTROK = false;

% Resets the reference population
handles.ReferencePopulation = 1;

handles = ProjectTitleText_CreateFcn(handles.ProjectTitleText,eventdata,handles);

% Also resets the number of subjects variable and associated text
set(handles.Dimensionality_Text, 'String','_ frames x _ voxels (_)');
handles.n_subjects = {};

% Resets the number of datasets entered to 0
handles.n_datasets = 0;

% Empties the data, motion, brain information and mask variables
handles.TC = {};
handles.FD = {};
handles.mask = {};
handles.brain_info = {};

% Resets the text related to motion and data files
handles.SubjNames = {};
handles.MotName = {};

% Resets the title and save folder information
handles.Log = {};

% Project title, by default 'Untitled'
handles.project_title = 'Untitled';

set(handles.LoadButton,'BackgroundColor',[51 75 163]/255);

% Directory to which data is to be saved (initially loaded as ./SavedData)
handles.savedir = fullfile(pwd,'SavedData');
set(handles.SaveFolderText,'String',handles.savedir);

%%%%%%%%%% Putting the loading part (bottom) back to normal %%%%%%%%%%%

% We also want to set the TR textbox back to its initial state
handles = TR_Entry_CreateFcn(handles.TR_Entry, eventdata, handles);



% Clears the content of section 1 only
function handles = ClearSection2(eventdata, handles)

% Puts back the seed buttons information to original state
handles.seed = [];
set(handles.SeedButton,'BackgroundColor',[204,146,146]/255);
set(handles.SeedButton,'Enable','off');
set(handles.PlotSeedButton,'Enable','off');
set(handles.PlotSeedButton,'Visible','off');

% Seed label entries set invisible
set(handles.S_SEED1,'Visible','off');
set(handles.S_SEED2,'Visible','off');
set(handles.S_SEED3,'Visible','off');

% Puts back the logical defining the type of seed information entered
handles.isSeedSubjectSpecific = 0;

% Removes graph display for the seed
cla(handles.SeedGraphX);
cla(handles.SeedGraphZ);
set(handles.SeedGraphX,'Visible','off');
set(handles.SeedGraphZ,'Visible','off');
set(handles.SliderX,'Visible','off');
set(handles.SliderZ,'Visible','off');
set(handles.XCoordText,'Visible','off');
set(handles.ZCoordText,'Visible','off');

%%%%%%%%%%%% Putting the seed map part back to normal %%%%%%%%%%%%%%%%%%%

% Resets the variable containing the seed maps of the subjects
handles.AvgSeedMap = [];

% Not clickable anymore
set(handles.SeedMapPushButton,'Enable','off');
set(handles.SeedMapPushButton,'Visible','off');

% Resets colorbar display
handles = ResetGraphDisplay(handles.ColorbarSeed,handles);

% Makes the slider and the text linked to slider of the seed map threshold
% back to invisible
set(handles.TSeed_Slider,'Visible','off');
set(handles.TSeed,'Visible','off');

% Resets graphs with seed map plots
handles = ResetGraphDisplay(handles.SeedMapX,handles);
handles = ResetGraphDisplay(handles.SeedMapZ,handles);

% Resets associated sliders
set(handles.SeedMapSliderX,'Visible','off');
set(handles.SeedMapSliderZ,'Visible','off');

% Resets associated slider texts
set(handles.SeedMap_SliderX,'Visible','off');
set(handles.SeedMap_SliderZ,'Visible','off');

% Resets the circles plot
handles = ResetGraphDisplay(handles.FancyCircles,handles);

% Sets the associated text back to invisible
set(handles.SeedPlusText,'Visible','off');
set(handles.SeedMinusText,'Visible','off');
set(handles.Seed1Text,'Visible','off');
set(handles.Seed2Text,'Visible','off');
set(handles.Seed3Text,'Visible','off');

% Puts the seed boxes back to not visible
set(handles.CheckS1POS,'Visible','off');
set(handles.CheckS2POS,'Visible','off');
set(handles.CheckS3POS,'Visible','off');
set(handles.CheckS1NEG,'Visible','off');
set(handles.CheckS2NEG,'Visible','off');
set(handles.CheckS3NEG,'Visible','off');

set(handles.TPSelectionButton,'Enable','off');
set(handles.TPSelectionButton,'Visible','off');

set(handles.PRadio,'Visible','off');
set(handles.TRadio,'Visible','off');
set(handles.uibuttongroup7,'Visible','off');
set(handles.TText,'Visible','off');
set(handles.TMotText,'Visible','off');
set(handles.TEdit,'Visible','off');
set(handles.TMotEdit,'Visible','off');

% Resets the frame selection mode
handles.SelMode = 'Threshold';

% Invisible Seed type list
set(handles.SeedPopup,'Visible','off');

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

% Resets the variables indexing seed selection time points retained
handles.idx_sep_seeds = {};
handles.sfrac = [];

% Resets the violin plot with percentage retained frames
handles = ResetGraphDisplay(handles.TPViolin,handles);



% Clears the content of section 3 only
function handles = ClearSection3(eventdata, handles)

set(handles.ClusterButton,'Enable','off');
set(handles.CCButton,'Enable','off');
set(handles.AssignButton,'Enable','off');
set(handles.AssignButton,'Visible','off');
    
set(handles.CCPlot,'Visible','off');
cla(handles.CCPlot);

set(handles.CAP_TP,'Visible','off');
set(handles.Percentile_Edit,'Visible','off');

set(handles.CAP_Kmax,'Visible','off');
set(handles.CAP_PPC,'Visible','off');
set(handles.CAP_N,'Visible','off');
set(handles.CAP_K,'Visible','off');
set(handles.CAP_PP,'Visible','off');
set(handles.CAP_PN,'Visible','off');

set(handles.KRange_Edit,'Visible','off');
set(handles.PCC_Edit,'Visible','off');
set(handles.ClusterEdit,'Visible','off');
set(handles.ClusterRepEdit,'Visible','off');
set(handles.ClusterPpEdit,'Visible','off');
set(handles.ClusterPnEdit,'Visible','off');

% Resets the consensus clustering parameter input boxes
handles = KRange_Edit_CreateFcn(handles.KRange_Edit,eventdata,handles);
handles = PCC_Edit_CreateFcn(handles.PCC_Edit,eventdata,handles);

% Resets the parameter input boxes
handles = ClusterEdit_CreateFcn(handles.ClusterEdit,eventdata,handles);
handles = ClusterRepEdit_CreateFcn(handles.ClusterRepEdit,eventdata,handles);
handles = ClusterPpEdit_CreateFcn(handles.ClusterPpEdit,eventdata,handles);
handles = ClusterPnEdit_CreateFcn(handles.ClusterPnEdit,eventdata,handles);
handles = Percentile_Edit_CreateFcn(handles.Percentile_Edit,eventdata,handles);

% Resets the consensus clustering parameters themselves
handles.Kmax = 12;
handles.PCC = 80;

set(handles.CCButton,'Visible','off');
set(handles.ClusterButton,'Visible','off');

set(handles.PIE_S1,'Visible','off');
set(handles.PIE_S2,'Visible','off');
set(handles.PIE_S3,'Visible','off');
set(handles.PIE_S1S2,'Visible','off');
set(handles.PIE_S2S3,'Visible','off');
set(handles.PIE_S1S3,'Visible','off');
set(handles.PIE_S1S2S3,'Visible','off');

% Resets the parameters themselves
handles.K = 5;
handles.n_rep = 20;
handles.Pp = 100;
handles.Pn = 100;
handles.percentile = 5;

% Resets the CAP parameters (CAPs, standard deviation within CAPs and
% indices of the CAPs to which all retained frames were assigned)
handles.CAP = [];
handles.STDCAP = [];
handles.idx = {};

% Resets the graph display of the CAP colorbar
handles = ResetGraphDisplay(handles.ColorbarCAP,handles);

% Reset all graph displays for the CAPs
tmpX = {handles.CAP1X,handles.CAP2X,handles.CAP3X,handles.CAP4X,handles.CAP5X};
tmpY = {handles.CAP1Y,handles.CAP2Y,handles.CAP3Y,handles.CAP4Y,handles.CAP5Y};
tmpZ = {handles.CAP1Z,handles.CAP2Z,handles.CAP3Z,handles.CAP4Z,handles.CAP5Z};
tmpF = {handles.CAP1_Frames,handles.CAP2_Frames,handles.CAP3_Frames,handles.CAP4_Frames,handles.CAP5_Frames};


for i_CAP = 1:5
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

% Resets the pie charts
handles = ResetGraphDisplay(handles.pie1,handles);
handles = ResetGraphDisplay(handles.pie2,handles);
handles = ResetGraphDisplay(handles.pie3,handles);
handles = ResetGraphDisplay(handles.pie4,handles);
handles = ResetGraphDisplay(handles.pie5,handles);

handles = ResetGraphDisplay(handles.ColorbarSimMat,handles);
handles = ResetGraphDisplay(handles.CAP_Mat,handles);



% Clears the content of section 4 only
function handles = ClearSection4(eventdata, handles)

set(handles.MetricsButton,'Enable','off');
set(handles.MetricsButton,'Visible','off');

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
handles = ResetGraphDisplay(handles.ColorbarTransMat,handles);

% Resets all the graphs from the metrics part
handles = ResetGraphDisplay(handles.TMGraph,handles);
handles = ResetGraphDisplay(handles.TM_Subject,handles);
handles = ResetGraphDisplay(handles.DynStates,handles);
handles = ResetGraphDisplay(handles.CumStates,handles);
handles = ResetGraphDisplay(handles.ViolinCounts,handles);
handles = ResetGraphDisplay(handles.ViolinCountsFrac,handles);
handles = ResetGraphDisplay(handles.ViolinNumber,handles);
handles = ResetGraphDisplay(handles.ViolinDuration,handles);

% Removes all the labels linked to Metrics displays
set(handles.DS_Scrubbed,'Visible','off');
set(handles.DS_NotSelected,'Visible','off');
set(handles.DS_Unassigned,'Visible','off');

tmp = {handles.DS_CAP1,handles.DS_CAP2,handles.DS_CAP3,handles.DS_CAP4,...
    handles.DS_CAP5,handles.DS_CAP6,handles.DS_CAP7,handles.DS_CAP8,...
    handles.DS_CAP9,handles.DS_CAP10,handles.DS_CAP11,handles.DS_CAP12};

for i = 1:length(tmp)
    set(tmp{i},'Visible','off');
end

clear tmp

tmp = {handles.V_POP1,handles.V_POP2,handles.V_POP3,handles.V_POP4};

for i = 1:length(tmp)
    set(tmp{i},'Visible','off');
end

clear tmp


% The following functions enable to modify the text of the project title
function ProjectTitleText_ButtonDownFcn(hObject, eventdata, handles)

set(hObject,'Enable','on');
set(hObject,'String','');
set(hObject,'FontAngle','normal');
uicontrol(hObject);

guidata(hObject,handles);


function handles = ProjectTitleText_CreateFcn(hObject, eventdata, handles)

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
    set(hObject,'BackgroundColor', [101,140,196]/255);
    
    handles.Log = CAP_AddToLog(handles.Log,'Valid project title entered',{handles.project_title},{'New project title'});
    
% If we haven't entered anything, the project is just named 'untitled'
else
    handles.project_title = 'Untitled';
    set(hObject,'BackgroundColor',[204,146,146]/255);
end

guidata(hObject, handles);





%% Consensus clustering functions to compute optimal K
% The functions below enable to determine the optimal values of cluster
% number into which to disentangle the CAPs data

% Executes when pressing on the button to run consensus clustering
function CCButton_Callback(hObject, eventdata, handles)

    % Computes the consensus results
    [Consensus] = CAP_ConsensusClustering(handles.Xonp{handles.ReferencePopulation},2:handles.Kmax,'items',handles.PCC/100,handles.n_rep,'correlation');

    % Calculates the quality metrics
    [~,Lorena] = ComputeClusteringQuality(Consensus,2:handles.Kmax);
    
    set(handles.CCPlot,'Visible','on');
    tmp_plot = bar(2:handles.Kmax,1-Lorena,'Parent',handles.CCPlot);
    xlabel(get(tmp_plot(1),'Parent'),'Cluster number K');
    ylabel(get(tmp_plot(1),'Parent'),'Stability');
    xlim(get(tmp_plot(1),'Parent'),[2-0.6,handles.Kmax+0.6]);
    ylim(get(tmp_plot(1),'Parent'),[0,1]);
    custom_cm = cbrewer('seq','Reds',25);
    colormap(handles.CCPlot,custom_cm(6:25,:));
    
guidata(hObject,handles);



function KRange_Edit_Callback(hObject, eventdata, handles)

if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 1) && (str2double(get(hObject,'String')) <= 12)
    handles.Kmax = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [101,140,196]/255);

    handles.Log = CAP_AddToLog(handles.Log,'Valid number of Kmax chosen',{handles.K},{'Max cluster number'});
    
else
    set(hObject,'BackgroundColor',[204,146,146]/255);
    handles.Kmax = 12;
end

guidata(hObject, handles);



function handles = KRange_Edit_CreateFcn(hObject, eventdata, handles)

set(hObject,'Enable','off');
set(hObject,'String','Click to enter...');
set(hObject,'FontAngle','italic');

guidata(hObject, handles);



function KRange_Edit_ButtonDownFcn(hObject, eventdata, handles)

set(hObject,'Enable','on');
set(hObject,'String','');
set(hObject,'FontAngle','normal');
uicontrol(hObject);

guidata(hObject,handles);



function PCC_Edit_Callback(hObject, eventdata, handles)

if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 50) && (str2double(get(hObject,'String')) <= 100)
    handles.PCC = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [101,140,196]/255);

    handles.Log = CAP_AddToLog(handles.Log,'Valid percentage of items chosen',{handles.K},{'Percentage items to cluster'});
    
else
    set(hObject,'BackgroundColor',[204,146,146]/255);
    handles.PCC = 12;
end

guidata(hObject, handles);



function handles = PCC_Edit_CreateFcn(hObject, eventdata, handles)

set(hObject,'Enable','off');
set(hObject,'String','Click to enter...');
set(hObject,'FontAngle','italic');

guidata(hObject, handles);



function PCC_Edit_ButtonDownFcn(hObject, eventdata, handles)

set(hObject,'Enable','on');
set(hObject,'String','');
set(hObject,'FontAngle','normal');
uicontrol(hObject);

guidata(hObject,handles);





%% Editing of all the CAP generation parameters

% Number of clusters
function ClusterEdit_Callback(hObject, eventdata, handles)

if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 1) && (str2double(get(hObject,'String')) <= 12)
    handles.K = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [101,140,196]/255);

    handles.Log = CAP_AddToLog(handles.Log,'Valid number of clusters chosen',{handles.K},{'Number of clusters'});
    
else
    set(hObject,'BackgroundColor',[204,146,146]/255);
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



% Number of k-means repetitions
function ClusterRepEdit_Callback(hObject, eventdata, handles)

if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 0) && (str2double(get(hObject,'String')) <= 50)
    handles.n_rep = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [101,140,196]/255);

    handles.Log = CAP_AddToLog(handles.Log,'Valid number of replicates chosen',{handles.n_rep},{'Number of replicates'});
    
else
    set(hObject,'BackgroundColor',[204,146,146]/255);
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



% Percentage of positive-valued voxels to keep
function ClusterPpEdit_Callback(hObject, eventdata, handles)

if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 0) && (str2double(get(hObject,'String')) <= 100)
    handles.Pp = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [101,140,196]/255);

    handles.Log = CAP_AddToLog(handles.Log,'Valid percentage positive voxels chosen',{handles.Pp},{'Percentage positive voxels'});
    
else
    set(hObject,'BackgroundColor',[204,146,146]/255);
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



% Percentage of negative-valued voxels to keep
function ClusterPnEdit_Callback(hObject, eventdata, handles)

if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 0) && (str2double(get(hObject,'String')) <= 100)
    handles.Pn = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [101,140,196]/255);
    
    handles.Log = CAP_AddToLog(handles.Log,'Valid percentage negative voxels chosen',{handles.Pn},{'Percentage negative voxels'});
    
else
    set(hObject,'BackgroundColor',[204,146,146]/255);
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

% We perform clustering
[handles.CAP,~,handles.STDCAP,handles.idx{handles.ReferencePopulation},...
    handles.CorrDist,handles.sfrac] = Run_Clustering(cell2mat(handles.Xonp{handles.ReferencePopulation}),...
    handles.K,handles.mask{handles.ReferencePopulation},handles.brain_info{handles.ReferencePopulation},...
    handles.Pp,handles.Pn,handles.n_rep,handles.idx_sep_seeds{handles.ReferencePopulation},handles.SeedType);
    
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

% Computation of the similarity
SimMat = corr(handles.CAP',handles.CAP');
SimMat(isnan(SimMat))=0;

% Graph set visible, and plotting
handles = ResetGraphDisplay(handles.CAP_Mat,handles);
set(handles.CAP_Mat,'Visible','on');
imagesc(SimMat,'Parent',handles.CAP_Mat);

tmp_cb2 = cbrewer('div','RdBu',1000);

colormap(handles.CAP_Mat,flipud(tmp_cb2));

% Correlation ranges from -1 to 1, so this is what we make the graph
% colorbar vary within. We also make the graph square and remove the axes
caxis(handles.CAP_Mat,[-1 1]);
axis(handles.CAP_Mat,'square','on');
axis(handles.CAP_Mat,'off');

% Addition of the colorbar just below
set(handles.ColorbarSimMat,'Visible','on');
handles.ColorbarSimMat = Create_CAP_colorbar(-1,1,0.5,0,'',...
    handles.ColorbarSimMat,'Vertical','div','RdBu',1000);

% If using the 'Intersection' option...
if strcmp(handles.SeedType,'Intersection')
    
    % Custom colormap
    custom_cm = 1/255*[211,36,36;11,170,65;51,75,163;255,255,180;58,221,221;186,59,204;242,242,242];
    
    % Graph displays are stored in a common tmp_sfrac cell array
    tmp_sfrac = {handles.pie1,handles.pie2,handles.pie3,handles.pie4,...
        handles.pie5};
    
    % The pie charts for each cluster are created
    for cc = 1:min([handles.K,5])
        
        % Pie charts
        set(tmp_sfrac{cc},'Visible','on');
        for tt = 1:size(handles.sfrac,3)
            lab{tt} = '';
        end
        
        pie(tmp_sfrac{cc},realmin*ones(size(handles.sfrac,3),1)+squeeze(mean(handles.sfrac(:,cc,:),1)),lab);
        
        switch handles.n_seed
            case 1
                errordlg('You managed the impossible, congratulations!');
            case 2
                colormap(tmp_sfrac{cc},flipud(custom_cm));
                set(handles.PIE_S1,'Visible','on');
                set(handles.PIE_S2,'Visible','on');
                set(handles.PIE_S3,'Visible','on');
            case 3
                colormap(tmp_sfrac{cc},flipud(custom_cm));
                set(handles.PIE_S1,'Visible','on');
                set(handles.PIE_S2,'Visible','on');
                set(handles.PIE_S3,'Visible','on');
                set(handles.PIE_S2S3,'Visible','on');
                set(handles.PIE_S1S2,'Visible','on');
                set(handles.PIE_S1S3,'Visible','on');
                set(handles.PIE_S1S2S3,'Visible','on');
        end
    end
end

% Same for the slider for the visualization threshold
set(handles.TVIS,'Visible','on');
set(handles.TVIS_Slider,'Visible','on'); 
set(handles.TVIS,'String',['Tv: ',sprintf('%.2f',get(handles.TVIS_Slider,'Value'))]);

% Makes the colorbar for the CAPs visible
handles.ColorbarCAP = Create_CAP_colorbar(-1.5,1.5,0.5,get(handles.TVIS_Slider,'Value'),'',handles.ColorbarCAP,'Horizontal','div','RdBu',1000);
set(handles.ColorbarCAP,'Visible','on');

% Concatenates all CAP information into metavariables for easier subsequent
% changes
tmpX = {handles.CAP1X,handles.CAP2X,handles.CAP3X,handles.CAP4X,handles.CAP5X};
tmpY = {handles.CAP1Y,handles.CAP2Y,handles.CAP3Y,handles.CAP4Y,handles.CAP5Y};
tmpZ = {handles.CAP1Z,handles.CAP2Z,handles.CAP3Z,handles.CAP4Z,handles.CAP5Z};
tmpF = {handles.CAP1_Frames,handles.CAP2_Frames,handles.CAP3_Frames,handles.CAP4_Frames,handles.CAP5_Frames};

% For each CAP...
for i_CAP = 1:min([handles.K,5])
    
    % Clears the display for each dimension
    cla(tmpX{i_CAP});
    cla(tmpY{i_CAP});
    cla(tmpZ{i_CAP});
    
    % Plots the new slice for each dimension
    tmpX{i_CAP} = plot_slice(handles.CAP(i_CAP,:),...
        get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},...
        'X',get(handles.CAP_SliderX,'Value'),tmpX{i_CAP});

    tmpY{i_CAP} = plot_slice(handles.CAP(i_CAP,:),...
        get(handles.TVIS_Slider,'Value'),1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},...
        'Y',get(handles.CAP_SliderY,'Value'),tmpY{i_CAP});
    
    tmpZ{i_CAP} = plot_slice(handles.CAP(i_CAP,:),...
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
    set(handles.MetricsButton,'Visible','on');

% Else, we enable the assignment before enabling the metrics computation
elseif handles.n_datasets > 1
    set(handles.AssignButton,'Enable','on');
    set(handles.AssignButton,'Visible','on');
    
    set(handles.CAP_TP,'Visible','on');
    set(handles.Percentile_Edit,'Visible','on');
end

handles.Log = CAP_AddToLog(handles.Log,'Clustering performed',...
    {handles.ReferencePopulation,handles.K,handles.n_rep,handles.Pp,...
    handles.Pn},{'Reference group index',...
    'Number of clusters','Number of replicates',...
    'Percentage positive voxels','Percentage negative voxels'});

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

            [handles.Xonp{n_ds},handles.Xonn{n_ds},p,handles.FrameIndices{n_ds},handles.idx_sep_seeds{n_ds}] = ...
                CAP_find_activity(handles.TC{n_ds},handles.seed,handles.T,handles.FD{n_ds},handles.Tmot,handles.SelMode,handles.SeedType);

            handles.RetainedPercentage{n_ds} = p(4:5,:);
            
            tmp_computedTPsel = [tmp_computedTPsel,n_ds];
        end

        try
            handles.idx{n_ds} = CAP_AssignFrames(handles.CAP,cell2mat(handles.Xonp{n_ds}),handles.CorrDist,handles.percentile)';
        catch
            errordlg('You computed CAPs with a different CAP type compared to the one used now; please use the same CAP type !');
        end
    end
end

% We then enable the computation of metrics
set(handles.MetricsButton,'Enable','on');
set(handles.MetricsButton,'Visible','on');

handles.Log = CAP_AddToLog(handles.Log,'Frame assignment performed',...
    {handles.ReferencePopulation,num2str(tmp_computedTPsel),...
    num2str(tmp_notref)},{'Reference group index','Group indices for which frames were computed',...
    'Group indices for which frames were assigned'});

guidata(hObject, handles);





%% Parameter control: percentile to use for frame assignment
% This asks for the percentile to use in frame assignment (i.e. the
% threshold of correlation below which frames are left unassigned)

function Percentile_Edit_Callback(hObject, eventdata, handles)

if ~isempty(str2double(get(hObject,'String'))) && (str2double(get(hObject,'String')) > 0) && (str2double(get(hObject,'String')) <= 100)
    handles.percentile = str2double(get(hObject,'String'));
    set(hObject,'BackgroundColor', [101,140,196]/255);

    handles.Log = CAP_AddToLog(handles.Log,'Valid percentile chosen',{handles.percentile},{'Percentile'});
    
else
    set(hObject,'BackgroundColor',[204,146,146]/255);
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
tmp_struct = {handles.CAP1X,handles.CAP2X,handles.CAP3X,handles.CAP4X,handles.CAP5X};

for i_CAP = 1:min([handles.K,5])
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

guidata(hObject,handles);



% Y dimension slider
function CAP_SliderY_Callback(hObject, eventdata, handles)

set(handles.CAP_YC,'String',['Y: ',sprintf('%.2f',get(hObject,'Value'))]);
tmp_struct = {handles.CAP1Y,handles.CAP2Y,handles.CAP3Y,handles.CAP4Y,handles.CAP5Y};

for i_CAP = 1:min([handles.K,5])
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

guidata(hObject,handles);



% Z dimension slider
function CAP_SliderZ_Callback(hObject, eventdata, handles)

set(handles.CAP_ZC,'String',['Z: ',sprintf('%.2f',get(hObject,'Value'))]);
tmp_struct = {handles.CAP1Z,handles.CAP2Z,handles.CAP3Z,handles.CAP4Z,handles.CAP5Z};

for i_CAP = 1:min([handles.K,5])
   
    cla(tmp_struct{i_CAP});
    tmp_struct{i_CAP} = plot_slice(handles.CAP(i_CAP,:),get(handles.TVIS_Slider,'Value'),...
        1.5,handles.mask{handles.ReferencePopulation},handles.brain,handles.brain_info{handles.ReferencePopulation},'Z',get(hObject,'Value'),tmp_struct{i_CAP});
end

guidata(hObject,handles);



function CAP_SliderZ_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

guidata(hObject,handles);





%% Sliders for threshold visualization (CAP analysis)
% Again, we want to update the slices and the text if we change those
% sliders

function TVIS_Slider_Callback(hObject, eventdata, handles)

% The text is changed
set(handles.TVIS,'String',['Tv: ',sprintf('%.2f',get(hObject,'Value'))]);

% The colorbar graph is modified to suit the new threshold value
cla(handles.ColorbarCAP);
handles.ColorbarCAP = Create_CAP_colorbar(-1.5,1.5,0.5,get(hObject,'Value'),'',handles.ColorbarCAP,'Horizontal','div','RdBu',1000);

% The brain slices are replotted
tmpX = {handles.CAP1X,handles.CAP2X,handles.CAP3X,handles.CAP4X,handles.CAP5X};
tmpY = {handles.CAP1Y,handles.CAP2Y,handles.CAP3Y,handles.CAP4Y,handles.CAP5Y};
tmpZ = {handles.CAP1Z,handles.CAP2Z,handles.CAP3Z,handles.CAP4Z,handles.CAP5Z};

for i_CAP = 1:min([handles.K,5])
    
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
            handles.FrameIndices{n_ds}.scrubbed,...
            tmp_nclust,handles.TR);   
    catch
        
        errordlg('You tried computing metrics using parameter values different from the ones that were employed to generate CAPs; please check !');
    end
end

handles.Log = CAP_AddToLog(handles.Log,'Metrics computed',...
    {handles.n_datasets,handles.K,handles.TR},...
    {'Number of datasets','Number of clusters','TR'});

tmp_cb = cbrewer('seq','Greys',1000);

% 2. Transition matrix for all subjects together

tmp_toplot = squeeze(mean(handles.TM{handles.ReferencePopulation},3));
tmp_toplot = tmp_toplot(3:end,3:end);

% Make graph visible and plotting
handles = ResetGraphDisplay(handles.TMGraph,handles);
set(handles.TMGraph,'Visible','on');
imagesc(tmp_toplot,'Parent',handles.TMGraph);

colormap(handles.TMGraph,flipud(tmp_cb));

clear tmp_toplot

% Arbitrary setting of probability scale from 0 to 0.03
caxis(handles.TMGraph,[0 0.03]);
axis(handles.TMGraph,'square','on');
axis(handles.TMGraph,'off');

% 3. Transition matrix for one subject

tmp_toplot = squeeze(handles.TM{handles.ReferencePopulation}(:,:,get(handles.SubjectMenuMetrics,'Value')));
tmp_toplot = tmp_toplot(3:end,3:end);

% makes graph visible and plots the information given by the Subject popup
handles = ResetGraphDisplay(handles.TM_Subject,handles);
set(handles.TM_Subject,'Visible','on');
imagesc(tmp_toplot,...
    'Parent',handles.TM_Subject);

colormap(handles.TM_Subject,flipud(tmp_cb));

clear tmp_toplot

% Same setting for the axes as for the average graph
caxis(handles.TM_Subject,[0 0.03]);
axis(handles.TM_Subject,'square','on');
axis(handles.TM_Subject,'off');

% We then create the colorbar for both cases
set(handles.ColorbarTransMat,'Visible','on');
handles.ColorbarTransMat = Create_CAP_colorbar(0,0.03,0.01,0,'',...
    handles.ColorbarTransMat,'Vertical','seq','Greys',1000);

% Makes the subject menu visible
set(handles.SubjectMenuMetrics,'Visible','on');

% 4. Dynamic state plotting

% Makes the graph visible
handles = ResetGraphDisplay(handles.DynStates,handles);
set(handles.DynStates,'Visible','on');

% Concatenates information from the different datasets
tmp_toplot = [];

for i = 1:handles.n_datasets
    tmp_toplot = [tmp_toplot; handles.TPM{i}; 0*ones(5,handles.SubjSize.TP)];
end
tmp_toplot = tmp_toplot(1:end-5,:);

custom_cm = cbrewer('qual','Set1',handles.K+1);
custom_cm = [0.05,0.05,0.05;1,1,1;custom_cm];

% If the TR has been properly entered, the x-axis is time; else, it depicts
% time index. In any case, we plot the states
if handles.isTROK
    imagesc(tmp_toplot,'Parent',handles.DynStates);
    colormap(handles.DynStates,(custom_cm));
    xlabel(handles.DynStates,'Time [s]');
else
    imagesc(tmp_toplot,'Parent',handles.DynStates);
    colormap(handles.DynStates,(custom_cm));
    xlabel(handles.DynStates,'Time index [-]');
end

ylabel(handles.DynStates,'Subjects [-]');
axis(handles.DynStates,'off');
caxis(handles.DynStates,[-1,handles.K+1]);

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
            plot(1:size(handles.TPM{i},2),handles.TPMCum{i}(j,:),'Color',handles.PopColor{1}(i,:),...
                'Parent',handles.CumStates);
            hold(handles.CumStates,'on');
        end
    end

    for i = 1:handles.n_datasets
        % Then, we plot a bold average across subjects
        plot(1:size(handles.TPM{i},2),tmp_TPMCum{i},'Color',handles.PopColor{2}(i,:),...
            'LineWidth',2,'Parent',handles.CumStates); 
        xlabel(handles.CumStates,'Time index [-]','FontSize',10);
        xlim(handles.CumStates,[1,size(handles.TPM{i},2)]);
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
        xlabel(handles.CumStates,'Time [s]','FontSize',10);
        xlim(handles.CumStates,[0,(size(handles.TPM{i},2)-1)*handles.TR]);
    end
end


ylabel(handles.CumStates,'Cumul. sum [-]','FontSize',10);
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

% Makes the displays visible
set(handles.DS_Scrubbed,'Visible','on');
set(handles.DS_Scrubbed,'ForegroundColor',custom_cm(1,:));

set(handles.DS_NotSelected,'Visible','on');
set(handles.DS_NotSelected,'ForegroundColor',[0.9,0.9,0.9]);

set(handles.DS_Unassigned,'Visible','on');
set(handles.DS_Unassigned,'ForegroundColor',custom_cm(handles.K+3,:));

tmp = {handles.DS_CAP1,handles.DS_CAP2,handles.DS_CAP3,handles.DS_CAP4,...
    handles.DS_CAP5,handles.DS_CAP6,handles.DS_CAP7,handles.DS_CAP8,...
    handles.DS_CAP9,handles.DS_CAP10,handles.DS_CAP11,handles.DS_CAP12};

for i = 1:handles.K
    set(tmp{i},'Visible','on');
    set(tmp{i},'ForegroundColor',custom_cm(2+i,:));
end

clear tmp

tmp = {handles.V_POP1,handles.V_POP2,handles.V_POP3,handles.V_POP4};

for i = 1:handles.n_datasets
    set(tmp{i},'Visible','on');
end

clear tmp
 
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
    tmp_toplot = tmp_toplot(3:end,3:end);
    
    imagesc(tmp_toplot,'Parent',handles.TM_Subject);
    caxis(handles.TM_Subject,[0 0.03]);
    axis(handles.TM_Subject,'square','on');
    axis(handles.TM_Subject,'off');
    
    tmp_cb = cbrewer('seq','Greys',1000);
    colormap(handles.TM_Subject,flipud(tmp_cb));
    
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





%% General utilities

% Resets the display of a graph object
function handles = ResetGraphDisplay(Graph,handles)

cla(Graph);
set(Graph,'Visible','off');



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

            tmp = M{i}(:,3:3+n_states-1)';

            for j = 1:n_states
                M2(i+(j-1)*n_pop,1:size(tmp,2)) = tmp(j,:);
            end
            
            clear tmp
            
        case 'Duration'
            tmp = DiscardNaN(M{i}(:,3:3+n_states-1))';
            
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
