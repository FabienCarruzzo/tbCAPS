%% This script intends to validate the CAP toolbox proper functioning with
% simulated data of different types

%% 1. Determination of ground truth data parameters

% Number of states that we define
K = 4;

% Number of voxels (10x10x10)
V = [10 10 10];

% For now, the mask keeps all the voxels
mask = logical(ones(prod(V),1));

% Number of time points
T = 100;

% Number of subjects
S = 4;

% FD is said to be null for now
FD = zeros(T,S);

% We try a simple brain_info
brain_info.dim = [10 10 10];
brain_info.mat = [20 0 0 -90; 0 20 0 -90; 0 0 20 -90; 0 0 0 1];


%% 2. Determination of the ground truth co-activation patterns
% For the moment, we select four partially spatially overlapping states: a
% cylinder along Z, a cube along Z, and two stripe patterns


% We must first define the four states
Circle = zeros(10,10);
for i = -4.5:4.5
    for j = -4.5:4.5
        if sqrt((i*i+j*j))<3
            Circle(i+5.5,j+5.5)=1;
        else
            Circle(i+5.5,j+5.5)=0;
        end
    end
end

for z = 1:10
    State1(:,:,z) = Circle;
end

Square = zeros(10,10);
Square(1:5,1:5) = ones(5,5);

for z = 1:10
    State2(:,:,z) = Square;
end

Vertical = zeros(10,10);
for j = 1:10

    if mod(j,2)==0
        Vertical(:,j)=ones(10,1);
    end
end

for z = 1:10
    State3(:,:,z) = Vertical;
end

Horizontal = zeros(10,10);
for i = 1:10

    if mod(i,2)==0
        Horizontal(i,:)=ones(1,10);
    end
end

for z = 1:10
    State4(:,:,z) = Horizontal;
end

State1 = State1(:);
State2 = State2(:);
State3 = State3(:);
State4 = State4(:);

%% 3. Determination of seed parameters

% The seeda are created spatially
seed1 = zeros(10,10,10);
seed1(9,9,5) = 1;
seed1 = seed1(:);
seed1 = seed1(mask);
seed1 = logical(seed1);
seed2 = zeros(10,10,10);
seed2(1,5,8) = 1;
seed2 = seed2(:);
seed2 = seed2(mask);
seed2 = logical(seed2);
seed3 = zeros(10,10,10);
seed3(5,5,5) = 1;
seed3 = seed3(:);
seed3 = seed3(mask);
seed3 = logical(seed3);

State1(seed1) = 1;
State1(seed2) = 0;
State1(seed3) = 1;
State2(seed1) = 0;
State2(seed2) = 1;
State2(seed3) = 0;
State3(seed1) = 1;
State3(seed2) = 1;
State3(seed3) = 1;
State4(seed1) = 0;
State4(seed2) = 1;
State4(seed3) = 1;

States = {State1,State2,State3,State4};

% Now that we defined the four states, we must define time courses for the
% different subjects. We assume that at each time point, only one state can
% be entered
idx = [ones(1,10),3*ones(1,5),4*ones(1,5),ones(1,10),2*ones(1,10),ones(1,10),3*ones(1,5),4*ones(1,5),3*ones(1,5),4*ones(1,5),ones(1,10),2*ones(1,10),ones(1,10);...
    ones(1,10),3*ones(1,5),4*ones(1,5),ones(1,10),2*ones(1,10),ones(1,10),3*ones(1,5),4*ones(1,5),3*ones(1,5),4*ones(1,5),ones(1,10),2*ones(1,10),ones(1,10);...
    ones(1,10),3*ones(1,5),4*ones(1,5),ones(1,10),2*ones(1,10),ones(1,10),3*ones(1,5),4*ones(1,5),3*ones(1,5),4*ones(1,5),ones(1,10),2*ones(1,10),ones(1,10);...
    ones(1,10),3*ones(1,5),4*ones(1,5),ones(1,10),2*ones(1,10),ones(1,10),3*ones(1,5),4*ones(1,5),3*ones(1,5),4*ones(1,5),ones(1,10),2*ones(1,10),ones(1,10)];


% TC contains the above time courses
TC = cell(1,4);

for s = 1:S
    for t = 1:T
        TC{s}(t,:) = States{idx(s,t)}(:)';% + 0.1*rand(1,prod(V));
    end
end

save('TC','TC');
save('FD','FD');
save('mask','mask');
save('seed1','seed1');
save('seed2','seed2');
save('seed3','seed3');
save('brain_info','brain_info');