%% This script generates artificial data to test the outputs of the toolbox

% Number of subjects
S = 10;

% Number of voxels
V = 1000;

% Number of time points
T = 200;

% Number of CAPs
K = 4;

N_seed_voxels = 50;

% CAP patterns
CAP = [[ones(250,1);zeros(750,1)], [zeros(250,1);ones(250,1);zeros(500,1)],...
    [zeros(500,1);ones(250,1);zeros(250,1)], [zeros(750,1);ones(250,1)]];

% The first 50 voxels are the seed voxels
CAP = [ones(N_seed_voxels,4);CAP];
CAP = [zeros(V+N_seed_voxels,1),CAP];

% Transition probabilities betwen CAPs and baseline (row/column 1 =
% baseline, the other ones = CAPs
TM = [0.8,0.2,0,0,0; ...
    0,0,1,0,0;...
    0,0,0,0.5,0.5;...
    0.4,0,0,0.6,0;...
    0,0,0,0.2,0.8];

FD = zeros(T,S);
mask = logical(ones(V+N_seed_voxels,1));
brain_info = [];
Seed = logical([ones(N_seed_voxels,1);zeros(V,1)]);


% Generates the state in which we are
for s = 1:S
    s
    State(s,1) = 1;
    
    for t = 2:T
        State(s,t) = CAP_PickState(State(s,t-1),TM);
    end
end

% Generates the associated time courses
for s = 1:S
    for t = 1:T
        TC{s}(t,:) = CAP(:,State(s,t))';
    end
end

