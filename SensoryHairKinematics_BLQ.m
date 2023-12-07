%% Summary

% Load and smooth xyz coordinate data
% Make some plots to make sure things look okay
% Make a body centric coordinate system (CS) and a global CS
% Calculate pitch, velocity, wingbeat frequency, camber

% Code requires Statistics and Machine Learning Toolbox, Curve Fitting
% Toolbox, etc

%% Load and smooth xyz and other input data
clear
close all

% Specify which computer you're using 
comp = 1; %For Brooke's Computer use comp=1

dir1 = 'G:/My Drive/Swartz lab materials/SensHairPilotData/xyzpts/data/'; %Brooke's laptop path %Computer 1 directory with input .csv files;  %MAKE SURE THERE IS A "/" at end of path

if comp==1
    dinfo    = dir([dir1 '*.csv']); %define path to xyz data (folder of .csv files - one for each trial)
    trls     = {dinfo.name};    % Extract the filenames for each .csv file
elseif comp==2
    dinfo    = dir([dir2 '*.csv']);
    trls     = {dinfo.name};
end

numpts   = 16; % number of tracked points
numtrls  = length(trls); % The total number of trials to analyze
addpath(genpath('G:/Shared drives/Bat Flight Lab/Sensory Hairs 2021/raw/'))
master = readtable('SensHairNovTrials.xlsx');
d_unfilt = NaN(1000, numpts*3, numtrls); % Create an array of NaNs for the unfiltered xyz data - will be filled when data is loaded
d_sm     = NaN(size(d_unfilt)); % Create an array of NaNs for the smoothed xyz data - will be filled when data is smoothed
addpath(genpath('G:/My Drive/Swartz lab materials/SensHairPilotData/xyzpts/'))
y_axis_pts_11_11 = readmatrix('yaxisupstreamxyzpts.csv'); % the digitized points are 2x3 matrix with first row upstream point, second row downstream point, xyz coords
y_axis_pts_11_12 = readmatrix('yaxisupstream_11_12xyzpts.csv'); % choosing random video set from each day for this, SensHair_07BFB_frez_Y20211112H131211.149379000_F005.cine
%windspd = [3, 3, 3, 5, 5, 5, 5, 5, 3, 3, 5, 5, 5, 5, 3]; % manually assigned wind speeds associated with each trial in current order, must change when trial input changes
% need to add and incorporate wind speed per trial for velocity and other calcs
%bat = 'Bat6B535'; % or 'Bat07BFB' or whatever bat ID is


for i=1:size(trls,2) % Load and smooth data for one trial at a time, for all the trials that are in the folder
    % Load raw xyz data
if comp==1    
    d_temp          = readmatrix([dir1 trls{i} ]); % Store the xyz data for trial(i) in a temporary array
elseif comp==2
    d_temp          = readmatrix([dir2 trls{i} ]); 
end
    if size(d_temp,2) > 48
    disp('there are more columns/points here than there should be! look into following trials:')
    disp(i)
    end
    d_unfilt(1:length(d_temp),1:size(d_temp,2),i) = d_temp; % Store trial(i) in the d_unfilt (:,:,i)
    %d_unfilt(d_unfilt==0) = NaN; % Replace any zeros with NaNs

    clear d_temp f1 f2 p1 p2 i
    
end

% Set up condition categories for reference and use throughout code
HINT3 = []; % HINT = hair intact, no turbulence, 3 m/s
HINT5 = []; % HINT = hair intact, no turbulence, 5 m/s
HINT = []; % both 3 and 5 together
HIT3 = []; % HINT = hair intact, no turbulence, 3 m/s
HIT5 = []; % HIT = hair intact, turbulence, 5 m/s
HIT = []; % both 3 and 5 together
HRT3 = []; % HRT = hair removed, turbulence, 3 m/s
HRT5 = []; % HRT = hair removed, turbulence, 5 m/s
HRT = []; % both 3 and 5 together
HRNT3 = []; % HRNT = hair removed, no turbulence, 3 m/s
HRNT5 = []; % HRNT = hair removed, no turbulence, 5 m/s
date_11_11 = []; % November 11 2022
date_11_12 = []; % November 12 2022

for i = 1:numtrls
    tmptrl = char(trls(i));
    tmptrla = tmptrl(15:35);
    for j = 1:4:size(master,1)
        strmast = char(master{j,1});
        mycheck = strcmp(string(tmptrla),string(strmast(20:40)));  
        if mycheck == 1 && string(char(master{j,4}))=='intact' && string(char(master{j,7}))=='noTurb' && master{j,6}==3;
            HINT3 = [HINT3 i];
            HINT = [HINT i];
        elseif mycheck == 1 && string(char(master{j,4}))=='intact' && string(char(master{j,7}))=='noTurb' && master{j,6}==5;
            HINT5 = [HINT5 i];
            HINT = [HINT i];
        elseif mycheck == 1 && string(char(master{j,4}))=='intact' && string(char(master{j,7}))=='turb' && master{j,6}==3;
            HIT3 = [HIT3 i];
            HIT = [HIT i];
        elseif mycheck == 1 && string(char(master{j,4}))=='intact' && string(char(master{j,7}))=='turb' && master{j,6}==5;
            HIT5 = [HIT5 i];
            HIT = [HIT i];
        elseif mycheck == 1 && string(char(master{j,4}))=='removed' && string(char(master{j,7}))=='turb' && master{j,6}==3;
            HRT3 = [HRT3 i];
            HRT = [HRT i];
        elseif mycheck == 1 && string(char(master{j,4}))=='removed' && string(char(master{j,7}))=='turb' && master{j,6}==5;
            HRT5 = [HRT5 i];
            HRT = [HRT i];
        elseif mycheck == 1 && string(char(master{j,4}))=='removed' && string(char(master{j,7}))=='noTurb' && master{j,6}==3;
            HRNT3 = [HRNT3 i];
            HRNT = [HRNT i];
        elseif mycheck == 1 && string(char(master{j,4}))=='removed' && string(char(master{j,7}))=='noTurb' && master{j,6}==5;
            HRNT5 = [HRNT5 i];
            HRNT = [HRNT i];
        else
            % do nothing
        end
    end
end
col = [repmat('m',1,length(HINT)), repmat('c',1,length(HIT)), repmat('b',1,length(HRT))]; % specifying colors for HINT then HIT then HRT (change w/ num trls)


% set which trials correspond with which dates for later use
for i = 1:numtrls
    tmptrl = char(trls(i));
    tmptrla = tmptrl(15:35);
    for j = 1:4:size(master,1)
        strmast = char(master{j,1});
        mycheck = strcmp(string(tmptrla),string(strmast(20:40)));  
        if mycheck == 1 && string(char(master{j,11}))=='11';
            date_11_11 = [date_11_11 i];
        elseif mycheck == 1 && string(char(master{j,11}))=='12';
            date_11_12 = [date_11_12 i];
        else
            % do nothing
        end
    end
end

% Smooth raw xyz data 
for i = 1:size(trls,2)
    for k=1:numpts*3
        startvals(k) = find(~isnan(d_unfilt(:,k,i))==1,1,'first'); % determine first non-NaN row in the trial (NaNs mess up the smoothing) 
        endvals(k)   = find(~isnan(d_unfilt(:,k,i))==1,1,'last');  % determine last non-NaN row trial 
        d_sm(startvals(k):endvals(k),k,i) = smooth(d_unfilt(startvals(k):endvals(k),k,i)); % smooth and fill gaps of unfiltered data, save as d_sm
    end
    f1(i) = min(startvals);
    f2(i) = max(endvals);
    clear startvals endvals
end

% set frame numbers for start of the second downstroke, since we now have
% start of first downstroke (f1) and start of third downstroke (f2)
ssd = [];
for i = 1:numtrls
    tmptrl = char(trls(i));
    tmptrla = tmptrl(15:35);
    for j = 1:4:size(master,1)
        strmast = char(master{j,1});
        mycheck = strcmp(string(tmptrla),string(strmast(20:40)));  
        if mycheck == 1;
            ssd(i) = master{j,17};        
        else
            % do nothing
        end
    end
end

%% Test plots 
% Plot unfiltered raw data (in black) and smoothed data (in red) to see the results
pl = 'r'; % 'all' plots all trials; 'r' plots a range of trials 

if pl=='y'
    for i=1:3%size(trls,2)
        figure,plot(d_unfilt(:,:,i))%, hold on, plot(d_sm(:,:,i),'r')
        title(trls{i})
    end
elseif pl=='r'
    for i=[10:11];
        spec = 1:3;
        q = d_unfilt(:,spec,i); %can substitute in q and qsin plot function as first input to plot subset of columns, probably most useful to plot in sets of three like 1:3 to get x y and z for one point
        qs = d_sm(:,spec,i); 
        figure,plot(q,'-*'), hold on, plot(qs,'s')
        legend
        title(trls{i})
    end 
elseif pl=='3D'
    for j=1:size(trls,2)
        figure, grid on, axis equal, hold on
        title({'Kinematic Data'; trls{j}})
        d=d_unfilt;
            for i=1:numpts
                hold on
                plot3(d(:,i*3-2,j),d(:,i*3-1,j),d(:,i*3,j),'-*')
            end
        %scatter3(0,0,0,'k')
        xlabel('X-Axis'), ylabel('Y-Axis'), zlabel('Z-Axis')
    end
end

% check out left wingtip over time for trial i
for i=11
        figure
        %plot3(d(:,1,i),d(:,2,i),d(:,3,i),'-*')
        plot(d_sm(:,1:3,i))
        title(trls{i})
end
% yellow = (up and down), red = (forward and back),
% blue =  (left and right)

%% Point Matrix Legend
%   1. 1:3     Left wingtip
%   2. 4:6     Left wrist
%   3. 7:9     Left 5th digit 
%   4. 10:12   Left elbow
%   5. 13:15   Left shoulder
%   6. 16:18   Left ankle
%   7. 19:21   Noseleaf tip
%   8. 22:24   Sternum (midline, ~1 cm inferior to shoulders)
%   9. 25:27   Tail base (midline, where fur meets uropatagium)
%   10. 28:30  Right shoulder
%   11. 31:33  Right ankle
%   12. 34:36  Right elbow
%   12. 37:39  Right wrist
%   14. 40:42  Right 5th digit
%   15. 43:45  Right wingtip
%   16. 46:48  Midline lower (midline, approx halfway between 8 [sternum]
%   and 9 [tail base])
%   17. 49:51  Midline shoulder (midline, halfway between tracked shoulder
%   points, calculated in matlab, not digitized, not in orig csv file)


%% Tranlate data to body-centered coordinate system (BCS)
% INPUT: 
%   - A: input matrix of order [n, 3]
%   - i/j/k: vectors of order [n, 3] that establish the old coordinate system
%            these must be orthogonal to each other and obey the right hand rule.
%   - u/v/w: vectors of order [n, 3] that establish the new coordinate system
%            these must be orthogonal to each other and obey the right hand rule.
%   - tijk: translation vector of order [n,3] based on the ijk coordinate system

% using shoulders and sternum to define axes currently
% toward right wing is + x
% perpendicular to vector through shoulders toward bat's head is + y 
% up from the bat's back, orthogonal to x and y (basically against gravity) is + z

addpath(genpath('G:\My Drive\Swartz lab materials\SensHairPilotData'));

A = NaN(1000,3,numtrls); %create empty vector to fill in for loop
tijk = NaN(1000,3,numtrls); %create empty vector to fill in for loop

bodycs = zeros(size(d_sm));
for t=1:size(trls,2)
    for jj=1:3:size(d_sm,2)
%         f1 = find(~isnan(d_unfilt(:,1,t))==1,1,'first'); % determine first non-NaN row in the trial (NaNs mess up the smoothing) 
%         f2 = find(~isnan(d_unfilt(:,1,t))==1,1,'last');  % determine last non-NaN row trial 

        A = d_sm(f1(t):f2(t),jj:jj+2,t); % n,3 % matrix of order [n,3] with n = f2-f1 (frames without NaNs)
        P1 = d_sm(f1(t):f2(t),13:15,t); %left shoulder pt
        P2 = d_sm(f1(t):f2(t),28:30,t); %right shoulder pt
        P3 = d_sm(f1(t):f2(t),22:24,t); %sternum (superior) pt

        u = (P2 - P1)./vecnorm(P2-P1,2,2);
        P4 = (P1+P2)/2;
        vs = (P4-P3)./vecnorm(P4-P3,2,2);
        w = cross(u,vs,2);
        w = w./vecnorm(w,2,2);
        v = cross(w,u,2);
        v = v./vecnorm(v,2,2);

        tijk = P3;
        
        bodycs(f1(t):f2(t),jj:jj+2,t) = changeOfCoordinateSys(A, u, v, w, tijk);
    end
end


% % make subsets of hair conditions and turbulence conditions for BCS
% d_HINT_BCS = NaN(1000,numpts*3,numtrls);
% d_HIT_BCS = NaN(1000,numpts*3,numtrls);
% d_HRT_BCS = NaN(1000,numpts*3,numtrls);
% %d_HRNT_BCS = NaN(1000,numpts*3,numtrls); % nothing in pilot test in this category
% 
% for t = HINT
%     d_HINT_BCS(:,:,t) = d_sm(:,:,t);
% end
% for t = HIT
%     d_HIT_BCS(:,:,t) = d_sm(:,:,t);
% end
% for t = HRT
%     d_HRT_BCS(:,:,t) = d_sm(:,:,t);
% end

%% Translate data to global coordinate system (GCS)
% +z = against gravity
% +y = upstream
% +x = toward right window of wind tunnel

A_HINT = NaN(1000,3,length(HINT)); %create empty vectors to fill in for loop
A_HIT = NaN(1000,3,length(HIT));
A_HRT = NaN(1000,3,length(HRT));
tijk_HINT = zeros(1000,3,length(HINT)); 
tijk_HIT = zeros(1000,3,length(HIT)); 
tijk_HRT = zeros(1000,3,length(HRT)); 
globalcs_HINT = zeros(1000,numpts*3,length(HINT));
globalcs_HIT = zeros(1000,numpts*3,length(HIT));
globalcs_HRT = zeros(1000,numpts*3,length(HRT));

% Create separate datasets for global coordinate system data by condition,
% and with correct y axis pts for the date of recordings of the conditions
for t=HINT
    for jj=1:3:size(d_sm,2)
        A_HINT = d_sm(f1(t):f2(t),jj:jj+2,t); % n,3 % matrix of order [n,3] with n = f2-f1 (frames without NaNs) and 3 as xyz
        w = repmat([0,0,1],1000,1);
        v = repmat((y_axis_pts_11_11(1,2:4)-y_axis_pts_11_11(2,2:4))./norm(y_axis_pts_11_11(1,2:4)-y_axis_pts_11_11(2,2:4)),1000,1);
        u = cross(v,w);
        tijk_HINT = repmat([0,0,0],length(f1(t):f2(t)),1); %sternum (superior) pt
        globalcs_HINT(f1(t):f2(t),jj:jj+2,t) = changeOfCoordinateSys(A_HINT, u, v, w, tijk_HINT);
        tmptrl = char(trls(t));
        tmptrla = tmptrl(15:35);
        for j = 1:4:size(master,1)
            strmast = char(master{j,1});
            mycheck = strcmp(string(tmptrla),string(strmast(20:40)));  
            if mycheck == 1 && string(char(master{j,11}))~='11'
            disp('dates are not as expected for these conditions! fix!')
            end
        end
    end
end
globalcs_HINT = globalcs_HINT(:,:,HINT);
for t=HIT
    for jj=1:3:size(d_sm,2)
        A_HIT = d_sm(f1(t):f2(t),jj:jj+2,t); % n,3 % matrix of order [n,3] with n = f2-f1 (frames without NaNs)
        w = repmat([0,0,1],1000,1);
        v = repmat((y_axis_pts_11_12(1,2:4)-y_axis_pts_11_12(2,2:4))./norm(y_axis_pts_11_12(1,2:4)-y_axis_pts_11_12(2,2:4)),1000,1);
        u = cross(v,w);
        tijk_HIT = repmat([0,0,0],length(f1(t):f2(t)),1); %sternum (superior) pt  
        globalcs_HIT(f1(t):f2(t),jj:jj+2,t) = changeOfCoordinateSys(A_HIT, u, v, w, tijk_HIT);
        tmptrl = char(trls(t));
        tmptrla = tmptrl(15:35);
        for j = 1:4:size(master,1)
            strmast = char(master{j,1});
            mycheck = strcmp(string(tmptrla),string(strmast(20:40)));  
            if mycheck == 1 && string(char(master{j,11}))~='12'
            disp('dates are not as expected for these conditions! fix!')
            end
        end
    end
end
globalcs_HIT = globalcs_HIT(:,:,HIT);
for t=HRT
    for jj=1:3:size(d_sm,2)
        A_HRT = d_sm(f1(t):f2(t),jj:jj+2,t); % n,3 % matrix of order [n,3] with n = f2-f1 (frames without NaNs)
        w = repmat([0,0,1],1000,1);
        v = repmat((y_axis_pts_11_12(1,2:4)-y_axis_pts_11_12(2,2:4))./norm(y_axis_pts_11_12(1,2:4)-y_axis_pts_11_12(2,2:4)),1000,1);
        u = cross(v,w);
        tijk_HRT = repmat([0,0,0],length(f1(t):f2(t)),1); %sternum (superior) pt     
        globalcs_HRT(f1(t):f2(t),jj:jj+2,t) = changeOfCoordinateSys(A_HRT, u, v, w, tijk_HRT);
        tmptrl = char(trls(t));
        tmptrla = tmptrl(15:35);
        for j = 1:4:size(master,1)
            strmast = char(master{j,1});
            mycheck = strcmp(string(tmptrla),string(strmast(20:40)));  
            if mycheck == 1 && string(char(master{j,11}))~='12'
            disp('dates are not as expected for these conditions! fix!')
            end
        end
    end
end
globalcs_HRT = globalcs_HRT(:,:,HRT);

%Add xyz columns with calculated midpoint between tracked shoulder markers on midline
for i = 1:length(HINT)  % global CS with midline shoulder point added
    for k=1:3
    globalcs_HINT(:,48+k,i) = (globalcs_HINT(:,12+k,i) + globalcs_HINT(:,27+k,i))/2; 
    end
end
for i = 1:length(HIT)  % global CS with midline shoulder point added
    for k=1:3
    globalcs_HIT(:,48+k,i) = (globalcs_HIT(:,12+k,i) + globalcs_HIT(:,27+k,i))/2; 
    end
end
for i = 1:length(HRT)  % global CS with midline shoulder point added
    for k=1:3
    globalcs_HRT(:,48+k,i) = (globalcs_HRT(:,12+k,i) + globalcs_HRT(:,27+k,i))/2; 
    end
end
numpts = 17; %update numpts to reflect added point


%% Look at 3d plot to see if global coordinate system seems correct
        figure, grid on, axis equal, hold on
        d=globalcs_HINT;
        o=d_sm;
            %for i=1:numpts
            i = 1
            j = 1
                hold on
                title('left wingtip, global coordinate system vs old smoothed');
                plot3(d(:,i*3-2,j),d(:,i*3-1,j),d(:,i*3,j),'-*')
                hold on
                plot3(o(:,i*3-2,j),o(:,i*3-1,j),o(:,i*3,j),'-*')   
                xlabel('X-Axis'), ylabel('Y-Axis'), zlabel('Z-Axis')
            %end
            

%% Calculate body pitch (tilt angle of line from midline shoulder to tail base)

body_HINT = NaN(size(globalcs_HINT,1), 3, size(globalcs_HINT,3)); 
body_HIT = NaN(size(globalcs_HIT,1), 3, size(globalcs_HIT,3)); 
body_HRT = NaN(size(globalcs_HRT,1), 3, size(globalcs_HRT,3)); 
bodytilt_HINT = NaN(size(globalcs_HINT,1), 1, size(globalcs_HINT,3));
bodytilt_HIT = NaN(size(globalcs_HIT,1), 1, size(globalcs_HIT,3)); 
bodytilt_HRT = NaN(size(globalcs_HRT,1), 1, size(globalcs_HRT,3)); 
z = [0,1];

    % BQ: calculate body angle tilt over whole trial with body defined as
    % line between shoulder midline and tail base; will redo with body defined
    % as line between halfway point between shoulders on midline and
    % tailbase
    for i=1:length(HINT)
        body_HINT(:,:,i) = globalcs_HINT(:,49:51,i) - globalcs_HINT(:,25:27,i);
        for j=1:length(globalcs_HINT)
            bodytilt_HINT(j,:,i) = 90 - acosd( (dot(body_HINT(j,[2,3],i), z)) / (norm(body_HINT(j,[2,3],i)) * norm(z)));
        end
    end
    for i=1:length(HIT)
        body_HIT(:,:,i) = globalcs_HIT(:,49:51,i) - globalcs_HIT(:,25:27,i);
        for j=1:length(globalcs_HIT)
            bodytilt_HIT(j,:,i) = 90 - acosd( (dot(body_HIT(j,[2,3],i), z)) / (norm(body_HIT(j,[2,3],i)) * norm(z)));
        end
    end
    for i=1:length(HRT)
        body_HRT(:,:,i) = globalcs_HRT(:,49:51,i) - globalcs_HRT(:,25:27,i);
        for j=1:length(globalcs_HRT)
            bodytilt_HRT(j,:,i) = 90 - acosd( (dot(body_HRT(j,[2,3],i), z)) / (norm(body_HRT(j,[2,3],i)) * norm(z)));
        end
    end
    
% quick plot
figure
hold on
plot(bodytilt_HRT(:,1,6));

% create new cell array with only digitized data, all leading and trailing
% NaNs removed
A2_HINT = cell(1,length(numtrls));
A2_HIT = cell(1,length(numtrls));
A2_HRT = cell(1,length(numtrls));

for i = 1:length(HINT)
    tilt_temp = bodytilt_HINT(:,:,i);
    A2_HINT{i} = {tilt_temp(f1(HINT(i)):f2(HINT(i))) transpose(1:length(tilt_temp(f1(HINT(i)):f2(HINT(i)))))};
    clear tilt_temp
end
for i = 1:length(HIT)
    tilt_temp = bodytilt_HIT(:,:,i);
    A2_HIT{i} = {tilt_temp(f1(HIT(i)):f2(HIT(i))) transpose(1:length(tilt_temp(f1(HIT(i)):f2(HIT(i)))))};
    clear tilt_temp
end
for i = 1:length(HRT) 
    tilt_temp = bodytilt_HRT(:,:,i);
    A2_HRT{i} = {tilt_temp(f1(HRT(i)):f2(HRT(i))) transpose(1:length(tilt_temp(f1(HRT(i)):f2(HRT(i)))))};
    clear tilt_temp
end

% normalize that digitized data section from 0 to 1 (0 = start of
% downstroke, 1 = start of 3rd downstroke, so this is full two wingbeats
% and 0.5 may not be exactly start of 2nd downstroke, can make this more
% precise later)
normbodytilt_HINT = cell(1,length(numtrls));
normbodytilt_HIT = cell(1,length(numtrls));
normbodytilt_HRT = cell(1,length(numtrls));

for i = 1:length(HINT)
    normbodytilt_HINT{i} = {mat2gray(A2_HINT{1,i}{1,2}) A2_HINT{1,i}{1,1}};
end
for i = 1:length(HIT)
    normbodytilt_HIT{i} = {mat2gray(A2_HIT{1,i}{1,2}) A2_HIT{1,i}{1,1}};
end
for i = 1:length(HRT)
    normbodytilt_HRT{i} = {mat2gray(A2_HRT{1,i}{1,2}) A2_HRT{1,i}{1,1}};
end

% split up body tilt first and second wingbeats so I can overlay them


% make some plots of body tilt (pitch)
figure
hold on
for i=1:length(HINT) % or HINT or HIT or HRT
    plot(normbodytilt_HINT{1,i}{1,1}, normbodytilt_HINT{1,i}{1,2},'color',col(HINT(i)))
    hold on
end
for i=1:length(HIT)
    plot(normbodytilt_HIT{1,i}{1,1}, normbodytilt_HIT{1,i}{1,2},'color',col(HIT(i)))
    hold on
end
for i=1:length(HRT)
    plot(normbodytilt_HRT{1,i}{1,1}, normbodytilt_HRT{1,i}{1,2},'color',col(HRT(i)))
    hold on
end
%legend
ylabel('Angle')
xlabel('Wingbeat cycle (0, 0.5, and 1 = start of downstrokes)')
title('Body tilt angle: pink=hair intact no turb; light blue=hair intact turb; dark blue=hair removed turb')
hold off 


% hair intact no turbulence
for i=1:length(HINT)
    varbodytilt_HINT(i,1) = var(normbodytilt_HINT{1,i}{1,2}, 'omitnan');
end
avar_bodytilt_HINT = mean(varbodytilt_HINT);
iqr_bodytilt_HINT = iqr(varbodytilt_HINT);
maxtomin_tilt_HINT = NaN(length(numtrls),1);
for i=1:length(HINT)
    maxtomin_tilt_HINT(i,1) = max(normbodytilt_HINT{1,i}{1,2})-min(normbodytilt_HINT{1,i}{1,2});
end
mean_tilt_diff_HINT = mean(maxtomin_tilt_HINT,'omitnan');

% hair intact turbulence
for i=1:length(HIT)
    varbodytilt_HIT(i,1) = var(normbodytilt_HIT{1,i}{1,2}, 'omitnan');
end
avar_bodytilt_HIT = mean(varbodytilt_HIT);
iqr_bodytilt_HIT = iqr(varbodytilt_HIT);
maxtomin_tilt_HIT = NaN(length(numtrls),1);
for i=1:length(HIT)
    maxtomin_tilt_HIT(i,1) = max(normbodytilt_HIT{1,i}{1,2})-min(normbodytilt_HIT{1,i}{1,2});
end
mean_tilt_diff_HIT = mean(maxtomin_tilt_HIT,'omitnan');

% hair removed turbulence
for i=1:length(HRT)
    varbodytilt_HRT(i,1) = var(normbodytilt_HRT{1,i}{1,2}, 'omitnan');
end
avar_bodytilt_HRT = mean(varbodytilt_HRT);
iqr_bodytilt_HRT = iqr(varbodytilt_HRT);
maxtomin_tilt_HRT = NaN(length(numtrls),1);
for i=1:length(HRT)
    maxtomin_tilt_HRT(i,1) = max(normbodytilt_HRT{1,i}{1,2})-min(normbodytilt_HRT{1,i}{1,2});
end
mean_tilt_diff_HRT = mean(maxtomin_tilt_HRT,'omitnan');

% statistics on body tilt between conditions

% should probably remove maxtomin_tilt honestly, too easily skewed by outliers
[h,p,ci,stats] = ttest2(maxtomin_tilt_HINT(:,1), maxtomin_tilt_HIT(:,1))
[h,p,ci,stats] = ttest2(maxtomin_tilt_HINT(:,1), maxtomin_tilt_HRT(:,1))
[h,p,ci,stats] = ttest2(maxtomin_tilt_HIT(:,1), maxtomin_tilt_HRT(:,1))

[h,p,ci,stats] = ttest2(varbodytilt_HINT(:,1), varbodytilt_HIT(:,1))
[h,p,ci,stats] = ttest2(varbodytilt_HINT(:,1), varbodytilt_HRT(:,1))
[h,p,ci,stats] = ttest2(varbodytilt_HIT(:,1), varbodytilt_HRT(:,1))

% Plot max - min diffs in body tilt angle per trial
% grp = [ones(size(maxtomin_tilt_HINT)); 2*ones(length(maxtomin_tilt_HIT(HIT,1)),1); 3*ones(length(maxtomin_tilt_HRT(HRT,1)),1)];
% figure, boxplot([maxtomin_tilt_HINT(HINT,1);maxtomin_tilt_HIT(HIT,1); maxtomin_tilt_HRT(HRT,1)],grp)
% hold on
% xlabel('Condition: HINT   HIT   HRT')
% ylabel('Range in body tilt')

figure
plot(maxtomin_tilt_HINT)


%% Body pitch/tilt over one wingbeat only (stacked two wingbeats digitized per trial)
% body_HINT / HIT / HRT is defined in previous section of body tilt
% calculations

% split up body tilt first and second wingbeats so I can overlay them
first_A2_HINT = cell(1,length(numtrls));
first_A2_HIT = cell(1,length(numtrls));
first_A2_HRT = cell(1,length(numtrls));
first_A2_HINT = cell(1,length(numtrls));
first_A2_HIT = cell(1,length(numtrls));
first_A2_HRT = cell(1,length(numtrls));
for i=1:length(HINT)
    first_A2_HINT{1,i}{1,1} = A2_HINT{1,i}{1,1}(1:ssd(HINT(i))-f1(HINT(i)),1);
    first_A2_HINT{1,i}{1,2} = A2_HINT{1,i}{1,2}(1:ssd(HINT(i))-f1(HINT(i)),1); %7 is problem
    second_A2_HINT{1,i}{1,1} = A2_HINT{1,i}{1,1}(ssd(HINT(i))-f1(HINT(i)):f2(HINT(i))-f1(HINT(i)),1);
    second_A2_HINT{1,i}{1,2} = A2_HINT{1,i}{1,2}(ssd(HINT(i))-f1(HINT(i)):f2(HINT(i))-f1(HINT(i)),1);  
end
for i=1:length(HIT)
    first_A2_HIT{1,i}{1,1} = A2_HIT{1,i}{1,1}(1:ssd(HIT(i))-f1(HIT(i)),1);
    first_A2_HIT{1,i}{1,2} = A2_HIT{1,i}{1,2}(1:ssd(HIT(i))-f1(HIT(i)),1);
    second_A2_HIT{1,i}{1,1} = A2_HIT{1,i}{1,1}(ssd(HIT(i))-f1(HIT(i)):f2(HIT(i))-f1(HIT(i)),1);
    second_A2_HIT{1,i}{1,2} = A2_HIT{1,i}{1,2}(ssd(HIT(i))-f1(HIT(i)):f2(HIT(i))-f1(HIT(i)),1);
end
for i=1:length(HRT)
    first_A2_HRT{1,i}{1,1} = A2_HRT{1,i}{1,1}(1:ssd(HRT(i))-f1(HRT(i)),1);
    first_A2_HRT{1,i}{1,2} = A2_HRT{1,i}{1,2}(1:ssd(HRT(i))-f1(HRT(i)),1);
    second_A2_HRT{1,i}{1,1} = A2_HRT{1,i}{1,1}(ssd(HRT(i))-f1(HRT(i)):f2(HRT(i))-f1(HRT(i)),1);
    second_A2_HRT{1,i}{1,2} = A2_HRT{1,i}{1,2}(ssd(HRT(i))-f1(HRT(i)):f2(HRT(i))-f1(HRT(i)),1);  
end

% Normalize
first_normbodytilt_HINT = cell(1,length(HINT));
for i = 1:length(HINT)
    first_normbodytilt_HINT{i} = {mat2gray(first_A2_HINT{1,i}{1,2}) first_A2_HINT{1,i}{1,1}};
end
first_normbodytilt_HIT = cell(1,length(HIT));
for i = 1:length(HIT)
    first_normbodytilt_HIT{i} = {mat2gray(first_A2_HIT{1,i}{1,2}) first_A2_HIT{1,i}{1,1}};
end
first_normbodytilt_HRT = cell(1,length(HRT));
for i = 1:length(HRT)
    first_normbodytilt_HRT{i} = {mat2gray(first_A2_HRT{1,i}{1,2}) first_A2_HRT{1,i}{1,1}};
end

second_normbodytilt_HINT = cell(1,length(HINT));
for i = 1:length(HINT)
    second_normbodytilt_HINT{i} = {mat2gray(second_A2_HINT{1,i}{1,2}) second_A2_HINT{1,i}{1,1}};
end
second_normbodytilt_HIT = cell(1,length(HIT));
for i = 1:length(HIT)
    second_normbodytilt_HIT{i} = {mat2gray(second_A2_HIT{1,i}{1,2}) second_A2_HIT{1,i}{1,1}};
end
second_normbodytilt_HRT = cell(1,length(HRT));
for i = 1:length(HRT)
    second_normbodytilt_HRT{i} = {mat2gray(second_A2_HRT{1,i}{1,2}) second_A2_HRT{1,i}{1,1}};
end

% make some plots of body tilt (pitch)
figure
hold on
for i=1:length(HINT) 
    plot(first_normbodytilt_HINT{1,i}{1,1}, first_normbodytilt_HINT{1,i}{1,2},'color',col(i))
    hold on
end
for i=1:length(HINT)
    plot(second_normbodytilt_HINT{1,i}{1,1}, second_normbodytilt_HINT{1,i}{1,2},'color',col(i))
    hold on
end
for i=1:length(HIT)
    plot(first_normbodytilt_HIT{1,i}{1,1}, first_normbodytilt_HIT{1,i}{1,2},'color',col(length(HINT)+i))
    hold on
end
for i=1:length(HIT)
    plot(second_normbodytilt_HIT{1,i}{1,1}, second_normbodytilt_HIT{1,i}{1,2},'color',col(length(HINT)+i))
    hold on
end
for i=1:length(HRT)
    plot(first_normbodytilt_HRT{1,i}{1,1}, first_normbodytilt_HRT{1,i}{1,2},'color',col(length(HINT)+length(HIT)+i))
    hold on
end
for i=1:length(HRT)
    plot(second_normbodytilt_HRT{1,i}{1,1}, second_normbodytilt_HRT{1,i}{1,2},'color',col(length(HINT)+length(HIT)+i))
    hold on
end
%legend
ylabel('Pitch angle (higher = pitched up')
xlabel('Wingbeat cycle (0 = start of downstroke, 0.5 = start of upstroke')
title('Body tilt angle: pink=hair intact no turb; light blue=hair intact turb; dark blue=hair removed turb')
hold off 

% check out means, iqrs, variance, etc
% hair intact no turbulence
for i=1:length(HINT)
    first_meanbodytilt_HINT(i,1) = mean(first_normbodytilt_HINT{1,i}{1,2},'omitnan');
    second_meanbodytilt_HINT(i,1) = mean(second_normbodytilt_HINT{1,i}{1,2},'omitnan');
    first_varbodytilt_HINT(i,1) = var(first_normbodytilt_HINT{1,i}{1,2}, 'omitnan');
    second_varbodytilt_HINT(i,1) = var(second_normbodytilt_HINT{1,i}{1,2}, 'omitnan');
    first_maxtomin_tilt_HINT(i,1) = max(first_normbodytilt_HINT{1,i}{1,2})-min(first_normbodytilt_HINT{1,i}{1,2});
    second_maxtomin_tilt_HINT(i,1) = max(second_normbodytilt_HINT{1,i}{1,2})-min(second_normbodytilt_HINT{1,i}{1,2});
end
comb_meanbodytilt_HINT = [first_meanbodytilt_HINT;second_meanbodytilt_HINT];
trtmnt_avg_bodytilt_HINT = mean(comb_meanbodytilt_HINT,'omitnan');
comb_varbodytilt_HINT = [first_varbodytilt_HINT;second_varbodytilt_HINT];
comb_avar_bodytilt_HINT = mean(comb_varbodytilt_HINT);
comb_iqr_bodytilt_HINT = iqr(comb_varbodytilt_HINT);
comb_maxtomin_tilt_HINT = NaN(length(numtrls),1);
comb_maxtomin_tilt_HINT = [first_maxtomin_tilt_HINT;second_maxtomin_tilt_HINT];
comb_mean_tilt_diff_HINT = mean(comb_maxtomin_tilt_HINT,'omitnan');



% hair intact turbulence
for i=1:length(HIT)
    first_meanbodytilt_HIT(i,1) = mean(first_normbodytilt_HIT{1,i}{1,2},'omitnan');
    second_meanbodytilt_HIT(i,1) = mean(second_normbodytilt_HIT{1,i}{1,2},'omitnan');
    first_varbodytilt_HIT(i,1) = var(first_normbodytilt_HIT{1,i}{1,2}, 'omitnan');
    second_varbodytilt_HIT(i,1) = var(second_normbodytilt_HIT{1,i}{1,2}, 'omitnan');
    first_maxtomin_tilt_HIT(i,1) = max(first_normbodytilt_HIT{1,i}{1,2})-min(first_normbodytilt_HIT{1,i}{1,2});
    second_maxtomin_tilt_HIT(i,1) = max(second_normbodytilt_HIT{1,i}{1,2})-min(second_normbodytilt_HIT{1,i}{1,2});
end
comb_meanbodytilt_HIT = [first_meanbodytilt_HIT;second_meanbodytilt_HIT];
trtmnt_avg_bodytilt_HIT = mean(comb_meanbodytilt_HIT, 'omitnan');
comb_varbodytilt_HIT = [first_varbodytilt_HIT;second_varbodytilt_HIT];
comb_avar_bodytilt_HIT = mean(comb_varbodytilt_HIT);
comb_iqr_bodytilt_HIT = iqr(comb_varbodytilt_HIT);
comb_maxtomin_tilt_HIT = NaN(length(numtrls),1);
comb_maxtomin_tilt_HIT = [first_maxtomin_tilt_HIT;second_maxtomin_tilt_HIT];
comb_mean_tilt_diff_HIT = mean(comb_maxtomin_tilt_HIT,'omitnan');


% hair removed turbulence
for i=1:length(HRT)
    first_meanbodytilt_HRT(i,1) = mean(first_normbodytilt_HRT{1,i}{1,2},'omitnan');
    second_meanbodytilt_HRT(i,1) = mean(second_normbodytilt_HRT{1,i}{1,2},'omitnan');
    first_varbodytilt_HRT(i,1) = var(first_normbodytilt_HRT{1,i}{1,2}, 'omitnan');
    second_varbodytilt_HRT(i,1) = var(second_normbodytilt_HRT{1,i}{1,2}, 'omitnan');
    first_maxtomin_tilt_HRT(i,1) = max(first_normbodytilt_HRT{1,i}{1,2})-min(first_normbodytilt_HRT{1,i}{1,2});
    second_maxtomin_tilt_HRT(i,1) = max(second_normbodytilt_HRT{1,i}{1,2})-min(second_normbodytilt_HRT{1,i}{1,2});
end
comb_meanbodytilt_HRT = [first_meanbodytilt_HRT;second_meanbodytilt_HRT];
trtmnt_avg_bodytilt_HRT = mean(comb_meanbodytilt_HRT,'omitnan');
comb_varbodytilt_HRT = [first_varbodytilt_HRT;second_varbodytilt_HRT];
comb_avar_bodytilt_HRT = mean(comb_varbodytilt_HRT);
comb_iqr_bodytilt_HRT = iqr(comb_varbodytilt_HRT);
comb_maxtomin_tilt_HRT = NaN(length(numtrls),1);
comb_maxtomin_tilt_HRT = [first_maxtomin_tilt_HRT;second_maxtomin_tilt_HRT];
comb_mean_tilt_diff_HRT = mean(comb_maxtomin_tilt_HRT,'omitnan');


% Stats here with combined first and second
% remove outliers first
comb_meanbodytilt_HIT(3) = [];
[h,p,ci,stats] = ttest2(comb_meanbodytilt_HINT(:,1), comb_meanbodytilt_HIT(:,1)) % statistical 0.06
[h,p,ci,stats] = ttest2(comb_meanbodytilt_HIT(:,1), comb_meanbodytilt_HRT(:,1)) % not statistical, 0.99
[h,p,ci,stats] = ttest2(comb_meanbodytilt_HINT(:,1), comb_meanbodytilt_HRT(:,1)) % statistical 0.06



% look at timing of minimum body tilt over wingbeat cycle
for i = 1:length(HINT)
    [c,d] = min(first_normbodytilt_HINT{1,i}{1,2});
    mintilt_pos_HINT(i) = first_normbodytilt_HINT{1,i}{1,1}(d);
    clear c d
end
for i = 1:length(HIT)
    [c,d] = min(first_normbodytilt_HIT{1,i}{1,2});
    mintilt_pos_HIT(i) = first_normbodytilt_HIT{1,i}{1,1}(d);
    clear c d
end
for i = 1:length(HRT)
    [c,d] = min(first_normbodytilt_HRT{1,i}{1,2});
    mintilt_pos_HRT(i) = first_normbodytilt_HRT{1,i}{1,1}(d);
    clear c d
end

% look at timing of maximum body tilt over wingbeat
for i = 1:length(HINT)
    [c,d] = min(first_normbodytilt_HINT{1,i}{1,2});
    mintilt_pos_HINT(i) = first_normbodytilt_HINT{1,i}{1,1}(d);
    clear c d
end
for i = 1:length(HIT)
    [c,d] = min(first_normbodytilt_HIT{1,i}{1,2});
    mintilt_pos_HIT(i) = first_normbodytilt_HIT{1,i}{1,1}(d);
    clear c d
end
for i = 1:length(HRT)
    [c,d] = min(first_normbodytilt_HRT{1,i}{1,2});
    mintilt_pos_HRT(i) = first_normbodytilt_HRT{1,i}{1,1}(d);
    clear c d
end

%% Calculate wingbeat frequency

% Find frame number of maxima and minima in z (vertical) axis
% insert 3 or 45 into wingtipz_1 and wingtipz_2 lines for column to extract
% from bodycs for left (3) or right (45) wingtip


for t = 1:size(trls,2)
    wingtipz_1(:,t) = bodycs(f1(t):f1(t)+25,3,t); % extract z component of left wingtip motion in BCS around start of downstroke digitization
    wingtipz_2(:,t) = bodycs(f2(t)-25:f2(t),3,t); % extract z component of left wingtip motion in BCS around start of third downstroke digitization (aka end of two full wingbeat cycles)
    if isnan(wingtipz_1(10,t)) 
        clear wingtipz_1(:,t)
        wingtipz_1(:,t) = bodycs(f1(t):f1(t)+25,45,t); % use right wingtip instead
        if isnan(wingtipz_1(1,t))
            clear wingtipz_1(:,t)
            wingtipz_1(:,t) = bodycs(f1(t):f1(t)+25,6,t); % use left wrist instead
            if isnan(wingtipz_1(1,t))
                wingtipz_1(:,t) = ones(26,1);
            end
        end
    end
    if isnan(wingtipz_2(10,t))
        clear wingtipz_1(:,t)
        wingtipz_2(:,t) = bodycs(f2(t)-25:f2(t),45,t);
        if isnan(wingtipz_1(1,t))
            clear wingtipz_1(:,t)
            wingtipz_2(:,t) = bodycs(f2(t)-25:f2(t),6,t); % use left wrist instead
            if isnan(wingtipz_2(1,t))
                wingtipz_2(:,t) = ones(26,1);
            end
        end
    end
    [maxi(t),maxframe(t)] = max(wingtipz_1(:,t));
    [nextmax(t),nextmaxframe(t)] = max(wingtipz_2(:,t));
    numfrm(t) = (f2(t)-25+nextmaxframe(t))-(f1(t)+maxframe(t));
    freq(t) = 1400 / numfrm(t); % (2 wingbeats / x frames) * (700 frames / 1 sec)
end

%%% testing ways to calculate WF from here to ...
% use frame numbers from original digitized data only
for t = 1:size(trls,2)
    fre(t) = 1400 / (f2(t)-f1(t))
end
% break up into conditions
fre_HINT = fre(HINT);
mean_fre_HINT = mean(fre_HINT);
fre_HIT = fre(HIT);
mean_fre_HIT = mean(fre_HIT);
fre_HRT = fre(HRT);
mean_fre_HRT = mean(fre_HRT);
% statistics for wingbeat frequency differences
[h,p,ci,stats] = ttest2(fre_HINT, fre_HIT)
[h,p,ci,stats] = ttest2(fre_HINT, fre_HRT) % now it's statistically significant
[h,p,ci,stats] = ttest2(fre_HIT, fre_HRT) % here too

for t = 1:size(trls,2)
    freq1(t) = 700 / (ssd(t)-f1(t));
    freq2(t) = 700 / (f2(t)-ssd(t));
    frequ(t) = (freq1(t) + freq2(t))/2;
end

freq1_HINT = freq1(HINT);
freq2_HINT = freq2(HINT);
% break up into conditions
frequ_HINT = frequ(HINT);
mean_frequ_HINT = mean(frequ_HINT);
frequ_HIT = frequ(HIT);
mean_frequ_HIT = mean(frequ_HIT);
frequ_HRT = frequ(HRT);
mean_frequ_HRT = mean(frequ_HRT);
% statistics for wingbeat frequency differences
[h,p,ci,stats] = ttest2(frequ_HINT, frequ_HIT)
[h,p,ci,stats] = ttest2(frequ_HINT, frequ_HRT) % now it's statistically significant
[h,p,ci,stats] = ttest2(frequ_HIT, frequ_HRT) % here too
%%% ...here


% break up into conditions
freq_HINT = freq(HINT);
mean_freq_HINT = mean(freq_HINT);
freq_HIT = freq(HIT);
mean_freq_HIT = mean(freq_HIT);
freq_HRT = freq(HRT);
mean_freq_HRT = mean(freq_HRT);

% statistics for wingbeat frequency differences
[h,p,ci,stats] = ttest2(freq_HINT, freq_HIT)
[h,p,ci,stats] = ttest2(freq_HINT, freq_HRT) 
[h,p,ci,stats] = ttest2(freq_HIT, freq_HRT)

% look at variation?
var(freq_HINT);
var(freq_HIT);
var(freq_HRT);

[H,P] = vartest2(freq_HINT, freq_HIT);
[H,P] = vartest2(freq_HIT, freq_HRT);
[H,P] = vartest2(freq_HRT, freq_HINT)

%mytest = nbintest(freq_HINT,freq_HIT);

% repeat with wrist (left) instead of wingtip for frequency calculation
%39 = z column for right wrist, 6 = z column for left wrist

for t = 1:size(trls,2)
    wristz_1(:,t) = bodycs(f1(t):f1(t)+25,6,t); % extract z component of left wingtip motion in BCS around start of downstroke digitization
    wristz_2(:,t) = bodycs(f2(t)-25:f2(t),6,t); % extract z component of left wingtip motion in BCS around start of third downstroke digitization (aka end of two full wingbeat cycles)
    [maxiw(t),maxframew(t)] = max(wristz_1(:,t));
    [nextmaxw(t),nextmaxframew(t)] = max(wristz_2(:,t));
    numfrmw(t) = (f2(t)-25+nextmaxframew(t))-(f1(t)+maxframew(t));
    freqw(t) = 1400 / numfrmw(t); % (2 wingbeats / x frames) * (700 frames / 1 sec)
end


% break up into conditions
freq_HINTw = freqw(HINT);
mean_freq_HINTw = mean(freq_HINTw);
freq_HITw = freqw(HIT);
mean_freq_HITw = mean(freq_HITw);
freq_HRTw = freqw(HRT);
mean_freq_HRTw = mean(freq_HRTw);

% statistics for wingbeat frequency differences
[h,p,ci,stats] = ttest2(freq_HINTw, freq_HITw)
[h,p,ci,stats] = ttest2(freq_HINTw, freq_HRTw) % still the same result (as of 12.14.22) with wrist
[h,p,ci,stats] = ttest2(freq_HITw, freq_HRTw)






%% Calculate basic velocity
% get velocity by dividing position / time for sternum marker
vel_HINT = NaN(400,length(HINT));
for i=HINT
    f1 = find(~isnan(d_unfilt(:,1,i))==1,1,'first'); 
    f2 = find(~isnan(d_unfilt(:,1,i))==1,1,'last');
    dist_str = squareform(pdist(globalcs_HINT(f1:f2,22:24,i)));
    dist_str_1 = NaN(length(f1:f2),1);
    for j=1:(length(dist_str)-1)
        dist_str_1(j,1) = dist_str(j+1,j);
    end
    vel_HINT(1:length(f1:f2),i) = dist_str_1*700;
    clear f1 f2 dist_str dist_str_1
end
vel_HIT = NaN(400,length(date_11_12));
for i=HIT
    f1 = find(~isnan(d_unfilt(:,1,i))==1,1,'first'); 
    f2 = find(~isnan(d_unfilt(:,1,i))==1,1,'last');
    dist_str = squareform(pdist(globalcs_HIT(f1:f2,22:24,i)));
    dist_str_1 = NaN(length(f1:f2),1);
    for j=1:(length(dist_str)-1)
        dist_str_1(j,1) = dist_str(j+1,j);
    end
    vel_HIT(1:length(f1:f2),i) = dist_str_1*700;
    clear f1 f2 dist_str dist_str_1
end
vel_HRT = NaN(400,length(date_11_12));
for i=HRT
    f1 = find(~isnan(d_unfilt(:,1,i))==1,1,'first'); 
    f2 = find(~isnan(d_unfilt(:,1,i))==1,1,'last');
    dist_str = squareform(pdist(globalcs_HRT(f1:f2,22:24,i)));
    dist_str_1 = NaN(length(f1:f2),1);
    for j=1:(length(dist_str)-1)
        dist_str_1(j,1) = dist_str(j+1,j);
    end
    vel_HRT(1:length(f1:f2),i) = dist_str_1*700;
    clear f1 f2 dist_str dist_str_1
end



%% Calculate (proxy for) camber
% left: 4:6 (left wrist) and 7:9 (left 5th digit)
% right: 37:39 (right wrist)and 40:42  (right 5th digit)
%mat2gray function used in this case to normalize such that max value from column
%equals 1, so every row value is the length of the chord at that time point
%(distance from wrist to tip of fifth digit) divided by the maximum chord
%length for that trial. So, higher value = more flat wing, less cambered.
%Lower value = shorter chord at that moment compared to max, so more
%cambered.
for i = HINT
    chordL_HINT(:,i) = sqrt((bodycs(:,4,i)-bodycs(:,7,i)).^2+(bodycs(:,5,i)-bodycs(:,8,i)).^2+(bodycs(:,6,i)-bodycs(:,9,i)).^2);
    chordR_HINT(:,i) = sqrt((bodycs(:,37,i)-bodycs(:,40,i)).^2+(bodycs(:,38,i)-bodycs(:,41,i)).^2+(bodycs(:,39,i)-bodycs(:,42,i)).^2);
    cambL_HINT(:,i) = mat2gray(chordL_HINT(:,i));     
    cambR_HINT(:,i) = mat2gray(chordR_HINT(:,i));
end
cambL_HINT = cambL_HINT(:,HINT);
cambR_HINT = cambR_HINT(:,HINT);

for i = HIT
    chordL_HIT(:,i) = sqrt((bodycs(:,4,i)-bodycs(:,7,i)).^2+(bodycs(:,5,i)-bodycs(:,8,i)).^2+(bodycs(:,6,i)-bodycs(:,9,i)).^2);
    chordR_HIT(:,i) = sqrt((bodycs(:,37,i)-bodycs(:,40,i)).^2+(bodycs(:,38,i)-bodycs(:,41,i)).^2+(bodycs(:,39,i)-bodycs(:,42,i)).^2);
    cambL_HIT(:,i) = mat2gray(chordL_HIT(:,i));     
    cambR_HIT(:,i) = mat2gray(chordR_HIT(:,i));
end
cambL_HIT = cambL_HIT(:,HIT);
cambR_HIT = cambR_HIT(:,HIT);

for i = HRT
    chordL_HRT(:,i) = sqrt((bodycs(:,4,i)-bodycs(:,7,i)).^2+(bodycs(:,5,i)-bodycs(:,8,i)).^2+(bodycs(:,6,i)-bodycs(:,9,i)).^2);
    chordR_HRT(:,i) = sqrt((bodycs(:,37,i)-bodycs(:,40,i)).^2+(bodycs(:,38,i)-bodycs(:,41,i)).^2+(bodycs(:,39,i)-bodycs(:,42,i)).^2);
    cambL_HRT(:,i) = mat2gray(chordL_HRT(:,i));     
    cambR_HRT(:,i) = mat2gray(chordR_HRT(:,i));
end
cambL_HRT = cambL_HRT(:,HRT);
cambR_HRT = cambR_HRT(:,HRT);

% normalize camber by wingbeat cycle (2 wingbeats, so 50% is end of a full
% wingbeat) for plotting. code adds the frame numbers as a percentage of
% the maximum frame number for that digitized trial into added columns, so
% you can plot that as the x axis and the camber values (which are already
% normalized to 1 = minimum camber for that trial) on the y axis. 
for i = 1:length(HINT)
    j = HINT(i);
    normcambL_HINT(1:length(f1(j):f2(j)),i) = cambL_HINT(f1(j):f2(j),i);
    normcambL_HINT(1:f2(j)-f1(j),length(HINT)+i) = 1:(f2(j)-f1(j));
    normcambL_HINT(1:f2(j)-f1(j),length(HINT)*2+i) = normcambL_HINT(1:f2(j)-f1(j),length(HINT)+i)/max(normcambL_HINT(1:f2(j)-f1(j),length(HINT)+i));
end
for i = 1:length(HIT)
    j = HIT(i);
    normcambL_HIT(1:length(f1(j):f2(j)),i) = cambL_HIT(f1(j):f2(j),i);
    normcambL_HIT(1:f2(j)-f1(j),length(HIT)+i) = 1:(f2(j)-f1(j));
    normcambL_HIT(1:f2(j)-f1(j),length(HIT)*2+i) = normcambL_HIT(1:f2(j)-f1(j),length(HIT)+i)/max(normcambL_HIT(1:f2(j)-f1(j),length(HIT)+i));
end
for i = 1:length(HRT)
    j = HRT(i);
    normcambL_HRT(1:length(f1(j):f2(j)),i) = cambL_HRT(f1(j):f2(j),i);
    normcambL_HRT(1:f2(j)-f1(j),length(HRT)+i) = 1:(f2(j)-f1(j));
    normcambL_HRT(1:f2(j)-f1(j),length(HRT)*2+i) = normcambL_HRT(1:f2(j)-f1(j),length(HRT)+i)/max(normcambL_HRT(1:f2(j)-f1(j),length(HRT)+i));
end
%repeat for right wing camber
for i = 1:length(HINT)
    j = HINT(i);
    normcambR_HINT(1:length(f1(j):f2(j)),i) = cambR_HINT(f1(j):f2(j),i);
    normcambR_HINT(1:f2(j)-f1(j),length(HINT)+i) = 1:(f2(j)-f1(j));
    normcambR_HINT(1:f2(j)-f1(j),length(HINT)*2+i) = normcambR_HINT(1:f2(j)-f1(j),length(HINT)+i)/max(normcambR_HINT(1:f2(j)-f1(j),length(HINT)+i));
end
for i = 1:length(HIT)
    j = HIT(i);
    normcambR_HIT(1:length(f1(j):f2(j)),i) = cambR_HIT(f1(j):f2(j),i);
    normcambR_HIT(1:f2(j)-f1(j),length(HIT)+i) = 1:(f2(j)-f1(j));
    normcambR_HIT(1:f2(j)-f1(j),length(HIT)*2+i) = normcambR_HIT(1:f2(j)-f1(j),length(HIT)+i)/max(normcambR_HIT(1:f2(j)-f1(j),length(HIT)+i));
end
for i = 1:length(HRT)
    j = HRT(i);
    normcambR_HRT(1:length(f1(j):f2(j)),i) = cambR_HRT(f1(j):f2(j),i);
    normcambR_HRT(1:f2(j)-f1(j),length(HRT)+i) = 1:(f2(j)-f1(j));
    normcambR_HRT(1:f2(j)-f1(j),length(HRT)*2+i) = normcambR_HRT(1:f2(j)-f1(j),length(HRT)+i)/max(normcambR_HRT(1:f2(j)-f1(j),length(HRT)+i));
end

% calculate position of maximum camber for each trial, compare between
% conditions (e.g. does max camber occur during downstroke in HINT and
% upstroke in HIT?)
% using min function to find max camber because we're looking for lowest
% value from mat2gray output, which makes max = 1, and the values we're
% looking at are chord length at time point compared to max chord length.
% more camber = smaller ratio
for i = 1:length(HINT)
    j = HINT(i);
    [a,b] = min(normcambL_HINT(1:(f2(j)-f1(j)),i));
    maxcambL_pos_HINT(i) = normcambL_HINT(b,length(HINT)*2+i);
    clear a b
end

for i = 1:length(HIT)
    j = HIT(i);
    [a,b] = min(normcambL_HIT(1:(f2(j)-f1(j)),i));
    maxcambL_pos_HIT(i) = normcambL_HIT(b,length(HIT)*2+i);
    clear a b
end

for i = 1:length(HRT)
    j = HRT(i);
    [a,b] = min(normcambL_HRT(1:(f2(j)-f1(j)),i));
    maxcambL_pos_HRT(i) = normcambL_HRT(b,length(HRT)*2+i);
    clear a b
end

% Repeat for right wing
for i = 1:length(HINT)
    j = HINT(i);
    [a,b] = min(normcambR_HINT(1:(f2(j)-f1(j)),i));
    maxcambR_pos_HINT(i) = normcambR_HINT(b,length(HINT)*2+i);
    clear a b
end

for i = 1:length(HIT)
    j = HIT(i);
    [a,b] = min(normcambR_HIT(1:(f2(j)-f1(j)),i));
    maxcambR_pos_HIT(i) = normcambR_HIT(b,length(HIT)*2+i);
    clear a b
end

for i = 1:length(HRT)
    j = HRT(i);
    [a,b] = min(normcambR_HRT(1:(f2(j)-f1(j)),i));
    maxcambR_pos_HRT(i) = normcambR_HRT(b,length(HRT)*2+i);
    clear a b
end

%% Plot camber two wingbeats
% plot camber over normalized frame range so 0% is start of downstroke and 100% is
% start of third downstroke (so end of two full wingbeats)
figure
for i = 1:length(HRT)
    j = HRT(i);
    f1 = find(~isnan(d_sm(:,1,j))==1,1,'first');
    f2 = find(~isnan(d_sm(:,1,j))==1,1,'last'); 
    plot(normcambR_HRT(1:f2-f1,length(HRT)*2+i),normcambR_HRT(1:f2-f1,i))
            hold on
            ylabel('Camber')
            xlabel('Percentage of wingbeat cycle (0,0.5,and 1=start of downstrokes)')
            title('Camber (right) Hair Removed Turbulence')
            hold on
end
hold off

figure
for i = HIT
    f1 = find(~isnan(d_sm(:,1,i))==1,1,'first');
    f2 = find(~isnan(d_sm(:,1,i))==1,1,'last'); 
    plot(normcambR_HIT(1:f2-f1,max(HIT)*2+i),normcambR_HIT(1:f2-f1,i))
            hold on
            ylabel('Camber')
            xlabel('Percentage of wingbeat cycle (0,0.5,and 1=start of downstrokes)')
            title('Camber (right) Hair Intact Turbulence')
            hold on
end
hold off

figure
for i = HINT
    f1 = find(~isnan(d_sm(:,1,i))==1,1,'first');
    f2 = find(~isnan(d_sm(:,1,i))==1,1,'last'); 
    plot(normcambR_HINT(1:f2-f1,max(HINT)*2+i),normcambR_HINT(1:f2-f1,i))
            hold on
            ylabel('Camber')
            xlabel('Percentage of wingbeat cycle (0,0.5,and 1=start of downstrokes)')
            title('Camber (right) Hair Intact No Turbulence')
            hold on
end
hold off

figure
hold on
for i=HINT
    plot(normcambL_HINT(1:f2-f1,max(HINT)*2+i),normcambL_HINT(1:f2-f1,i),'color',col(i))
    hold on
end
for i=HIT
    plot(normcambL_HIT(1:f2-f1,max(HIT)*2+i),normcambL_HIT(1:f2-f1,i),'color',col(i))
    hold on
end
for i=HRT % or HINT or HIT or HRT
    plot(normcambL_HRT(1:f2-f1,max(HRT)*2+i),normcambL_HRT(1:f2-f1,i),'color',col(i))
    hold on
end
ylabel('Normalized proxy for camber (1=max length of left wrist to 5th digit)')
xlabel('Wingbeat cycle (0, 0.5, and 1 = start of downstrokes')
title('Camber: pink=hair intact no turb; light blue=hair intact turb; dark blue=hair removed turb')
hold off 


[h,p,ci,stats] = ttest2(normcambL_HINT(:,1), normcambL_HINT(:,2));


%% Calculate camber over one wingbeat (break up two wingbeat digitizations and stack)
% normalize camber by wingbeat cycle after breaking up 2 wingbeats (so
% looking at one wingbeat, 50% is start of upstroke). 
% code adds the frame numbers as a percentage of the maximum frame number for that digitized trial into added columns, so
% you can plot that as the x axis and the camber values (which are already
% normalized to 1 = minimum cambering for that trial) on the y axis. 

for i = 1:length(HINT)
    j = HINT(i);
    first_normcambL_HINT(1:length(f1(j):ssd(j)),i) = cambL_HINT(f1(j):ssd(j),i);
    first_normcambL_HINT(1:ssd(j)-f1(j),length(HINT)+i) = 1:(ssd(j)-f1(j));
    first_normcambL_HINT(1:ssd(j)-f1(j),length(HINT)*2+i) = first_normcambL_HINT(1:ssd(j)-f1(j),length(HINT)+i)/max(first_normcambL_HINT(1:ssd(j)-f1(j),length(HINT)+i));
end
for i = 1:length(HINT)
    j = HINT(i);
    second_normcambL_HINT(1:length(ssd(j):f2(j)),i) = cambL_HINT(ssd(j):f2(j),i);
    second_normcambL_HINT(1:f2(j)-ssd(j),length(HINT)+i) = 1:(f2(j)-ssd(j));
    second_normcambL_HINT(1:f2(j)-ssd(j),length(HINT)*2+i) = second_normcambL_HINT(1:f2(j)-ssd(j),length(HINT)+i)/max(second_normcambL_HINT(1:f2(j)-ssd(j),length(HINT)+i));
end
comb_camber_HINT = {first_normcambL_HINT(:,1:length(HINT)*3), second_normcambL_HINT(:,1:length(HINT)*3)};

for i = 1:length(HIT)
    j = HIT(i);
    first_normcambL_HIT(1:length(f1(j):ssd(j)),i) = cambL_HIT(f1(j):ssd(j),i);
    first_normcambL_HIT(1:ssd(j)-f1(j),length(HIT)+i) = 1:(ssd(j)-f1(j));
    first_normcambL_HIT(1:ssd(j)-f1(j),length(HIT)*2+i) = first_normcambL_HIT(1:ssd(j)-f1(j),length(HIT)+i)/max(first_normcambL_HIT(1:ssd(j)-f1(j),length(HIT)+i));
end
for i = 1:length(HIT)
    j = HIT(i);
    second_normcambL_HIT(1:length(ssd(j):f2(j)),i) = cambL_HIT(ssd(j):f2(j),i);
    second_normcambL_HIT(1:f2(j)-ssd(j),length(HIT)+i) = 1:(f2(j)-ssd(j));
    second_normcambL_HIT(1:f2(j)-ssd(j),length(HIT)*2+i) = second_normcambL_HIT(1:f2(j)-ssd(j),length(HIT)+i)/max(second_normcambL_HIT(1:f2(j)-ssd(j),length(HIT)+i));
end
comb_camber_HIT = {first_normcambL_HIT(:,1:length(HIT)*3), second_normcambL_HIT(:,1:length(HIT)*3)}; 


for i = 1:length(HRT)
    j = HRT(i);
    first_normcambL_HRT(1:length(f1(j):ssd(j)),i) = cambL_HRT(f1(j):ssd(j),i);
    first_normcambL_HRT(1:ssd(j)-f1(j),length(HRT)+i) = 1:(ssd(j)-f1(j));
    first_normcambL_HRT(1:ssd(j)-f1(j),length(HRT)*2+i) = first_normcambL_HRT(1:ssd(j)-f1(j),length(HRT)+i)/max(first_normcambL_HRT(1:ssd(j)-f1(j),length(HRT)+i));
end
for i = 1:length(HRT)
    j = HRT(i);
    second_normcambL_HRT(1:length(ssd(j):f2(j)),i) = cambL_HRT(ssd(j):f2(j),i);
    second_normcambL_HRT(1:f2(j)-ssd(j),length(HRT)+i) = 1:(f2(j)-ssd(j));
    second_normcambL_HRT(1:f2(j)-ssd(j),length(HRT)*2+i) = second_normcambL_HRT(1:f2(j)-ssd(j),length(HRT)+i)/max(second_normcambL_HRT(1:f2(j)-ssd(j),length(HRT)+i));
end
comb_camber_HRT = {first_normcambL_HRT(:,1:length(HRT)*3), second_normcambL_HRT(:,1:length(HRT)*3)};


% Find position of maximum camber (which is minimum value in this calc)
for i = 1:length(HINT)
    j = HINT(i);
    [e,f] = min(comb_camber_HINT{1,1}(1:ssd(j)-f1(j),i));
    [g,h] = min(comb_camber_HINT{1,2}(1:f2(j)-ssd(j),i));
    if e < g 
        comb_maxcambL_pos_HINT(1,i) = e;
    else
        comb_maxcambL_pos_HINT(1,i) = g;        
    end
    clear e f g h
end

for i = 1:length(HIT)
    j = HIT(i);
    [e,f] = min(comb_camber_HIT{1,1}(1:ssd(j)-f1(j),i));
    [g,h] = min(comb_camber_HIT{1,2}(1:f2(j)-ssd(j),i));
    if e < g 
        comb_maxcambL_pos_HIT(1,i) = e;
    else
        comb_maxcambL_pos_HIT(1,i) = g;        
    end
    clear e f g h
end

for i = 1:length(HRT)
    j = HRT(i);
    [e,f] = min(comb_camber_HRT{1,1}(1:ssd(j)-f1(j),i));
    [g,h] = min(comb_camber_HRT{1,2}(1:f2(j)-ssd(j),i));
    if e < g 
        comb_maxcambL_pos_HRT(1,i) = e;
    else
        comb_maxcambL_pos_HRT(1,i) = g;        
    end
    clear e f g h
end

%% Write Output Data to Directory
answer = questdlg('Would you like to write output files?','Output Data','Yes','No','Yes');
     
 switch answer
     case 'Yes'
     currDate = strrep(datestr(datetime), ':', '_');
    
         if comp==1
             outputDir = [dir1 'Output_' currDate '\'];
             mkdir([outputDir])
         elseif comp==2
             outputDir = [dir2 'Output_' currDate '/'];
             mkdir([outputDir])
         end
         
         for i=1:size(trls,2)
             %Write table for Raw Data in Flower-Coordinate-System
             %(filename example = output_FCS_20200129_1_1_control.csv)
                 tFCS = array2table(d_fcs_unfilt(:,:,i), 'VariableNames', ... 
                     {'Flower 12 oclock_X' 'Flower 12 oclock_Y' 'Flower 12 oclock_Z' ...
                      'Flower 10 oclock_X' 'Flower 10 oclock_Y' 'Flower 10 oclock_Z' ...
                      'Flower 9 oclock_X' 'Flower 9 oclock_Y' 'Flower 9 oclock_Z' ...
                      'Flower 7 oclock_X' 'Flower 7 oclock_Y' 'Flower 7 oclock_Z' ...
                      'Upper lip_X' 'Upper lip_Y' 'Upper lip_Z' 'Eye_X' 'Eye_Y' ...
                      'Eye_Z' 'Nose bridge_X' 'Nose bridge_Y' 'Nose bridge_Z' ...
                      'Noseleaf tip_X' 'Noseleaf tip_Y' 'Noseleaf tip_Z' ...
                      'Tip of tongue_X' 'Tip of tongue_Y' 'Tip of tongue_Z' ...
                      'Surface top_X' 'Surface top_Y' 'Surface top_Z' 'Surface bottom_X' ...
                      'Surface bottom_Y' 'Surface bottom_Z'});
                  writetable(tFCS, [outputDir 'output_FCS_' trls{i}]);
                  
             %Write table for Approach-Phase Raw Data in Flower-Coordinate-System
             %(filename example = output_Approach_FCS_20200129_1_1_control.csv)
                 tFCS = array2table(d_fcs_unfilt_Approach(:,:,i), 'VariableNames', ... 
                     {'Flower 12 oclock_X' 'Flower 12 oclock_Y' 'Flower 12 oclock_Z' ...
                      'Flower 10 oclock_X' 'Flower 10 oclock_Y' 'Flower 10 oclock_Z' ...
                      'Flower 9 oclock_X' 'Flower 9 oclock_Y' 'Flower 9 oclock_Z' ...
                      'Flower 7 oclock_X' 'Flower 7 oclock_Y' 'Flower 7 oclock_Z' ...
                      'Upper lip_X' 'Upper lip_Y' 'Upper lip_Z' 'Eye_X' 'Eye_Y' ...
                      'Eye_Z' 'Nose bridge_X' 'Nose bridge_Y' 'Nose bridge_Z' ...
                      'Noseleaf tip_X' 'Noseleaf tip_Y' 'Noseleaf tip_Z' ...
                      'Tip of tongue_X' 'Tip of tongue_Y' 'Tip of tongue_Z' ...
                      'Surface top_X' 'Surface top_Y' 'Surface top_Z' 'Surface bottom_X' ...
                      'Surface bottom_Y' 'Surface bottom_Z'});
                  writetable(tFCS, [outputDir 'output_Approach_FCS_' trls{i}]);    
              
             %Write table for in-FLower-Phase Raw Data in Flower-Coordinate-System
             %(filename example = output_inFlower_FCS_20200129_1_1_control.csv)
                 tFCS = array2table(d_fcs_unfilt_inFlower(:,:,i), 'VariableNames', ... 
                     {'Flower 12 oclock_X' 'Flower 12 oclock_Y' 'Flower 12 oclock_Z' ...
                      'Flower 10 oclock_X' 'Flower 10 oclock_Y' 'Flower 10 oclock_Z' ...
                      'Flower 9 oclock_X' 'Flower 9 oclock_Y' 'Flower 9 oclock_Z' ...
                      'Flower 7 oclock_X' 'Flower 7 oclock_Y' 'Flower 7 oclock_Z' ...
                      'Upper lip_X' 'Upper lip_Y' 'Upper lip_Z' 'Eye_X' 'Eye_Y' ...
                      'Eye_Z' 'Nose bridge_X' 'Nose bridge_Y' 'Nose bridge_Z' ...
                      'Noseleaf tip_X' 'Noseleaf tip_Y' 'Noseleaf tip_Z' ...
                      'Tip of tongue_X' 'Tip of tongue_Y' 'Tip of tongue_Z' ...
                      'Surface top_X' 'Surface top_Y' 'Surface top_Z' 'Surface bottom_X' ...
                      'Surface bottom_Y' 'Surface bottom_Z'});
                  writetable(tFCS, [outputDir 'output_inFlower_FCS_' trls{i}]);  
                  
              %Write table for Head Tilt Angle during Approach (when lip_x crosses ≤ 0, i.e. when the lip crosses into the flower)
              % (filename example = output_Approach_HeadTilt_20200129_1_1_control.csv)
                  tHTapp = array2table(headtilt_approach(:,:,i), 'VariableNames',{'HeadTilt_Approach'});
                  writetable(tHTapp, [outputDir 'output_Approach_HeadTilt_' trls{i}]);
  
              %Write table for Head Tilt Angle within the Flower (when lip_x > 0, i.e. within the flower)        
              % (filename example = output_inFlower_HeadTilt_20200129_1_1_control.csv)
              tHTinFl = array2table(headtilt_inFlower(:,:,i), 'VariableNames',{'HeadTilt_inFlower'});
              writetable(tHTinFl, [outputDir 'output_inFlower_HeadTilt_' trls{i}]);
              
              %Write table for 3D distance from origin (flower opening center) to all points at each frame
              % (filename example = output_3DptDistance_20200129_1_1_control.csv)
              tptDist = array2table(ptDist(:,:,i), 'VariableNames',{'Flower 12 oclock' ...
                   'Flower 10 oclock' 'Flower 9 oclock' 'Flower 7 oclock' ...
                   'Upper lip' 'Eye' 'Nose bridge' 'Noseleaf tip' ...
                   'Tip of tongue' 'Surface top' 'Surface bottom'});
              writetable(tptDist, [outputDir 'output_3DptDistance_' trls{i}]);
              
              %Write table for Approach-phase 3D distance from origin (flower opening center) to all points at each frame
              % (filename example = output_Approach_3DptDistance_20200129_1_1_control.csv)
              tptDist = array2table(ptDist_Approach(:,:,i), 'VariableNames',{'Flower 12 oclock' ...
                   'Flower 10 oclock' 'Flower 9 oclock' 'Flower 7 oclock' ...
                   'Upper lip' 'Eye' 'Nose bridge' 'Noseleaf tip' ...
                   'Tip of tongue' 'Surface top' 'Surface bottom'});
              writetable(tptDist, [outputDir 'output_Approach_3DptDistance_' trls{i}]);
              
              %Write table for in-FLower-phase 3D distance from origin (flower opening center) to all points at each frame
              % (filename example = output_inFlower_3DptDistance_20200129_1_1_control.csv)
              tptDist = array2table(ptDist_inFlower(:,:,i), 'VariableNames',{'Flower 12 oclock' ...
                   'Flower 10 oclock' 'Flower 9 oclock' 'Flower 7 oclock' ...
                   'Upper lip' 'Eye' 'Nose bridge' 'Noseleaf tip' ...
                   'Tip of tongue' 'Surface top' 'Surface bottom'});
              writetable(tptDist, [outputDir 'output_inFlower_3DptDistance_' trls{i}]);
              
              clear tFCS tHTapp tHTinFl tptDist
         end
         
         

     case 'No'
         disp('Ouput data not written to file')
 end
 