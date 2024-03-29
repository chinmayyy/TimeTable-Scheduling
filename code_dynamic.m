%% COURSE NAMES
course_names = ["BIO310", "BIO311", "BIO313", "BIO314", "BIO315", "BIO316", "BIO320", "BIO321", "BIO410", "BIO411", "BIO415", "CHM311", "CHM312", "CHM320", "ECS310", "ECS312", "ECS313", "ECS335", "ECS420", "HSS334", "MTH310", "MTH311", "MTH312", "MTH314", "MTH315", "MTH318", "MTH435", "MTH436", "PHY310", "PHY311", "PHY312", "PHY313", "BIO301", "BIO355", "BIO431", "BIO454", "CHM301", "CHM331", "CHM332", "CHM340", "ECS301", "ECS317", "ECS330", "ECS331", "HSS333", "HSS342", "HSS351", "MTH301", "PHY301", "PHY330", "PHY334", "PHY340", "PHY462"]'
%% SELECTIONS OF COURSES

%Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /media/backslash/2TB/Matlab/Chinmay/Timetable/Final Project/SemData_MATRIX.xlsx
%    Worksheet: SemReg-b12-5 (3)
%
% Auto-generated by MATLAB on 25-Jun-2019 17:37:14

% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 53);

% Specify sheet and range
opts.Sheet = "SemReg-b12-5 (3)";
opts.DataRange = "A3:BA208";

% Specify column names and types
opts.VariableNames = ["BIO310", "BIO311", "BIO313", "BIO314", "BIO315", "BIO316", "BIO320", "BIO321", "BIO410", "BIO411", "BIO415", "CHM311", "CHM312", "CHM320", "ECS310", "ECS312", "ECS313", "ECS335", "ECS420", "HSS334", "MTH310", "MTH311", "MTH312", "MTH314", "MTH315", "MTH318", "MTH435", "MTH436", "PHY310", "PHY311", "PHY312", "PHY313", "BIO301", "BIO355", "BIO431", "BIO454", "CHM301", "CHM331", "CHM332", "CHM340", "ECS301", "ECS317", "ECS330", "ECS331", "HSS333", "HSS342", "HSS351", "MTH301", "PHY301", "PHY330", "PHY334", "PHY340", "PHY462"];
opts.SelectedVariableNames = ["BIO310", "BIO311", "BIO313", "BIO314", "BIO315", "BIO316", "BIO320", "BIO321", "BIO410", "BIO411", "BIO415", "CHM311", "CHM312", "CHM320", "ECS310", "ECS312", "ECS313", "ECS335", "ECS420", "HSS334", "MTH310", "MTH311", "MTH312", "MTH314", "MTH315", "MTH318", "MTH435", "MTH436", "PHY310", "PHY311", "PHY312", "PHY313", "BIO301", "BIO355", "BIO431", "BIO454", "CHM301", "CHM331", "CHM332", "CHM340", "ECS301", "ECS317", "ECS330", "ECS331", "HSS333", "HSS342", "HSS351", "MTH301", "PHY301", "PHY330", "PHY334", "PHY340", "PHY462"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Setup rules for import
opts = setvaropts(opts, [1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53], "FillValue", 0);

% Import the data
course_selections = readtable("/media/backslash/2TB/Matlab/Chinmay/Timetable/Final Project/SemData_MATRIX.xlsx", opts, "UseExcel", false);

% Convert to output type
course_selections = table2array(course_selections);

%Clear temporary variables
clear opts


%% CREDITS AND COURSES

%Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /media/backslash/2TB/Matlab/Chinmay/Timetable/Final Project/SemData_MATRIX.xlsx
%    Worksheet: SemReg-b12-5 (3)
%
% Auto-generated by MATLAB on 25-Jun-2019 17:49:33

%Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 53);

% Specify sheet
opts.Sheet = "SemReg-b12-5 (3)";

% Specify column names and types
opts.VariableNames = ["BIO310", "BIO311", "BIO313", "BIO314", "BIO315", "BIO316", "BIO320", "BIO321", "BIO410", "BIO411", "BIO415", "CHM311", "CHM312", "CHM320", "ECS310", "ECS312", "ECS313", "ECS335", "ECS420", "HSS334", "MTH310", "MTH311", "MTH312", "MTH314", "MTH315", "MTH318", "MTH435", "MTH436", "PHY310", "PHY311", "PHY312", "PHY313", "BIO301", "BIO355", "BIO431", "BIO454", "CHM301", "CHM331", "CHM332", "CHM340", "ECS301", "ECS317", "ECS330", "ECS331", "HSS333", "HSS342", "HSS351", "MTH301", "PHY301", "PHY330", "PHY334", "PHY340", "PHY462"];
opts.SelectedVariableNames = ["BIO310", "BIO311", "BIO313", "BIO314", "BIO315", "BIO316", "BIO320", "BIO321", "BIO410", "BIO411", "BIO415", "CHM311", "CHM312", "CHM320", "ECS310", "ECS312", "ECS313", "ECS335", "ECS420", "HSS334", "MTH310", "MTH311", "MTH312", "MTH314", "MTH315", "MTH318", "MTH435", "MTH436", "PHY310", "PHY311", "PHY312", "PHY313", "BIO301", "BIO355", "BIO431", "BIO454", "CHM301", "CHM331", "CHM332", "CHM340", "ECS301", "ECS317", "ECS330", "ECS331", "HSS333", "HSS342", "HSS351", "MTH301", "PHY301", "PHY330", "PHY334", "PHY340", "PHY462"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
credits = table;
ranges = ["A1:BA1", "A3:BA3"];
for idx = 1:length(ranges)
    opts.DataRange = ranges(idx);
    tb = readtable("/media/backslash/2TB/Matlab/Chinmay/Timetable/Final Project/SemData_MATRIX.xlsx", opts, "UseExcel", false);
    credits = [credits; tb]; %#ok<AGROW>
end

% Convert to output type
credits = table2array(credits);

% Clear temporary variables
clear idx opts ranges tb
credits = credits'

%% GENERATE EVENTS ARRAY

n = sum(credits(:,2))
e = [1:n]
cred_select = []
creditnames = []
for i = 1:length(course_selections(1,:))
    credit = course_selections(1,i);
    for j= 1:credit
        cred_select = [cred_select course_selections(2:end,i)];
        cred_select;
        creditnames = [creditnames course_names(i)];
    end
end
creditnames = creditnames'
%% GENERATE TIMESLOTS ARRAY
% Assigned number of time slots per day = 6
% Hence, total number of timeslots in a week (k)= 30
% These are labelled from 1 to k, 1 being the first time slot of Monday,
% and k being the last time slot of Friday.
k = 30
t = [1:k]

%% GENERATE STUDENTS ARRAY
% Total number of students in the batch is given by |s|.
s_number = length(course_selections(2:end,1))
s = [1:s_number]

%% GENERATE CLASS-ROOMS ARRAY along with individual capacities
% (Estimated Guess of the number of rooms) 
% 4 tutorial rooms on the ground floor of LHC, 4 Class rooms on the first
% floor, 2 class rooms on the 2nd floor. Each of these have an approximate
% capacity of 80. LHC101 has a capacity of 150, and LHC 103 has a capacity
% of 250.
% Capacities of each room are given by array 'cr'. 
% Temporarily, assume all the capacities to be equal. 

numberof_rooms = 12
r = [1:numberof_rooms]

capacity = 80 .* ones(1,numberof_rooms)

%% P MATRICES
% Check out the user guide for details

% Attends Matrix (s x n)
P1 = cred_select

% Room features
%P2 = zeros(r,1) %% NOT NEEDED

% Event Features 
% P3     || NOT NEEDED

% Event Availability (n x k)
P4 = ones(n,k)

% Precedence
P5 = zeros(n,n)

%% MATRICES TO CHECK HARD CONSTRAINT VIOLATIONS

% Room Suitability Matrix (n x r)
R = zeros(n,numberof_rooms);
ss = sum(P1)
for i = 1:length(ss)
    for j = 1:numberof_rooms
        if ss(i) <= capacity(j)
            R(i,j) = 1;
        end
    end
end
R
% Conflicts Matrix (n , n)

C = zeros(n,n)
for l=1:s_number
    for i=1:n
        if P1(l,i)==1 
            for j=1:n
                if P1(l,j) == 1
                    C(i,j) = 1;
                    C(j,i) = 1;
                end
            end
        end
    end
end

for l = 1:numberof_rooms
    for i = 1:n
        for j = 1:n
            if R(i,l) == 1
                if R(j,l) == 1
                    if sum(R(i,:)')== 1
                        if sum(R(j,:)') == 1
                            C(i,j)=1;
                            C(j,i) = 1;
                        end
                    end
                end
            end
        end
    end
end

for i =1:n
    for j = 1:n
        if P5(i,j) ~= 0
            C(i,j) =1;
            C(j,i) = 1;
        end
    end
end

C = C .* (1 - eye(n,n));

spy(C)
title('Conflict Matrix')

G = graph(C);
deg = degree(G);

%C = C./repmat(sum(C),191,1)
%% DYNAMIC PHASE GROUPING MODEL

%rng(1)

Pinhi = C

l=length(Pinhi);

Pex=(Pinhi==0);
Pex = Pex .* (1-eye(l,l))

%Specify the initial phases of the vertices
z=rand(l,1)

% Perturbation Pre-requisites
normal_pert_size=.001;
instant_pert=.0;
store_coup=zeros(l,it);

%Specify epsilon 
e_ij_ex = ones(l,1).* (0.016);
e_ij_inhi = ones(l,1).* (-(0.02));

%Create an empty matrix to use in a graph later
v=zeros(l,1);

it=90000;           % number of iterations

spikes_out_aftercoupling = zeros(l,1);

for i = 1:it

    EVENT = 'THRESHOLD';

    threshold_time = max ( 1- max(z));       % the phase will be increased by this amount to ensure that atleast one of the phases fire.
    z = z + threshold_time;
    spikes_out = (z>=1)+ spikes_out_aftercoupling;

    z = z.* (spikes_out==0);                % neurons receiving spikes

    v(:,i) = z;

    EVENT = 'RECEPTION';

    spikes_in = spikes_out;         % access incoming spikes at this event
    Pinhi(1:4,1:4)=5*(1-eye(4,4));
    coupling_ex = e_ij_ex .* ((Pex*spikes_in));
    coupling_inhi = e_ij_inhi .* ((Pinhi*spikes_in));
    if i>90000-10
        spikes_in
        Pex*spikes_in
        coupling_ex
        coupling_inhi
    end
    coupling_final = coupling_ex + coupling_inhi;

    temp1 = (z.^(0.5) + coupling_final)>=0;
    temp2 = (z.^(0.5) + coupling_final)<0;
    temp3 = (z.^(0.5) + coupling_final)<1;
    % Check for spikes in this event
    temp4 = (z.^(0.5) + coupling_final)>=1;
    temp5 = (z.^(0.5) + coupling_final) .* ((temp1.*temp3)~=0);
    spikes_out_aftercoupling = temp4;

    z = (temp5).^2;
    
    EVENT = 'PERTURBATION';
    
    if mod(i,1)==0    % If i is divisible by __
        rn=normrnd(0,normal_pert_size,[l,1]);
        z=(z+rn).*((z+ rn)<1).*((z+ rn)>0)+z.*((z+rn)>1)+z.*((z+rn)<0);
    end
     
    %% Adding Perturbation (Code from Sandeep)
    %rng(1)
    %{
    x=rand(l,1);
    %x=x_dd;
    x_initial=x;

    normal_pert_size=.000;
    instant_pert=.0;
    store_coup=zeros(l,it);

    for j=1:it
         if mod(j,1)==0
            rn=normrnd(0,normal_pert_size,[l,1]);
              x=(x+rn).*((x+ rn)<1).*((x+ rn)>0)+x.*((x+rn)>1)+x.*((x+rn)<0);
         end
    end
    %}
end

figure
plot(sort(z))
xlabel('Vertex #')
ylabel('Phase')

z = [creditnames z]
clusters = size(unique(z(:,2)))


%{
%%
Dist_credits = zeros(length(z),length(z));
for i = 1:length(z)
    for j = 1:length(z)
        Dist_credits(i,j) = min(abs(z(i) - z(j)), 1 - abs(z(i) - z(j)));
    end
end

Z=linkage(Dist_credits);
%heatmap(A_sudoku_coarse)
%plot(0)
dendrogram(Z);
%axis([0 10 2 8])
T = cluster(Z,'cutoff',.1);
T_=[(1:length(z))' T];
T__=sortrows(T_,2)
T___=T__(:,1);
%}