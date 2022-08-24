%%% Freely Moving DLC Output Processing Script %%%
%%% Dewan Lab %%%
%%% Austin Pauley & Sam Caton %%%
%%% 8-16-2022 %%%

%% ===== Read Data File and Adjustments  ===== %%
clear; %%Refresh Everything

framerate = 30;
confidenceThreshold = 0.9;

%Import tracking data
[file, path] = uigetfile( '*.csv','Select Tracking Data: ', 'D:\');
opts = detectImportOptions(strcat(path,file));
opts = setvartype(opts,1,'double');
PoseData = readtable(strcat(path,file),opts);
PoseData.Properties.VariableNames = ["index","noseX", "noseY", "noseLike", "headX", "headY",  "headLike",...
     "LHeadbarX", "LHeadbarY", "LHeadbarLike", "RHeadbarX", "RHeadbarY", "RHeadbarLike", "bodyX", "bodyY", "bodyLike",...
     "rearX","rearY","rearLike","LEDX","LEDY","LEDLike"]; %%Resets all of the headers to make grabbing data easier; head/cannula/endoscope will all be called 'head' here for simplicity
PoseData = PoseData((3:end),:);

%Import Experimental H5
[file, path] = uigetfile('*.h5', 'Select Eperiment Output:', 'D:\');
ExperimentData = h5read(strcat(path,file),"/Trials");

%Extract odors from the H5
ImportOdors = ExperimentData.odor';
Odors = cell(length(ImportOdors),1);

for i = 1: length(Odors)
    Odors{i} = [ImportOdors(i,find(~isspace(ImportOdors(i,:))))];
end

%[file, path] = uigetfile('*.mp4', 'Select Eperiment Video:', 'D:\');

totalFrames = length(PoseData.index);

%% ======= Initial Settings ==== %%
initialSettings = inputApp;
initialSettings.Settings.Visible = 1;
uiwait(initialSettings.UIFigure);

%% ======= X/Y Coordinates ===== %%
noseCoords = [PoseData.noseX, PoseData.noseY];
headCoords = [PoseData.headX, PoseData.headY];
LHeadbarCoords = [PoseData.LHeadbarX, PoseData.LHeadbarY];
RHeadbarCoords = [PoseData.RHeadbarX, PoseData.RHeadbarY];
bodyCoords = [PoseData.bodyX, PoseData.bodyY];
rearCoords = [PoseData.rearX, PoseData.rearY];
LEDCoords = [PoseData.LEDX, PoseData.LEDY];

%% ======== Certainty Values ====== %

noseLike = [PoseData.noseLike];
headLike = [PoseData.headLike];
LHeadbarLike = [PoseData.LHeadbarLike];
RHeadbarLike = [PoseData.RHeadbarLike];
bodyLike = [PoseData.bodyLike];
rearLike = [PoseData.rearLike];
LEDLike = [PoseData.LEDLike];

%% ======== Bad Frames ======= %
%%Gather all of the indicies (frame number + 1) where the certainty is below our
%%defined threshold
badNoseFrames = find(noseLike < confidenceThreshold);
badHeadFrames = find(headLike < confidenceThreshold);
badLHFrames = find (LHeadbarLike < confidenceThreshold);
badRHFrames = find(RHeadbarLike < confidenceThreshold);
badBodyFrames = find(bodyLike < confidenceThreshold);
badRearFrames = find(rearLike < confidenceThreshold);

%% ======== LED Stats ======== %%
ledState = zeros(size(LEDCoords));
for i = 1:length(LEDCoords)
    if(LEDLike(i) >= 0.98)
        ledState(i,1) = 1; %% If led is on set the column 1 digit to 1
        if(LEDCoords(i,1) <= ledLeftX)
            ledState(i,2) = 0; %% If led is on the left, set column 2 digit to 0
        elseif(LEDCoords(i,1) >= ledRightX)
            ledState(i,2) = 1; %% If led is on the right, set column 2 digit to 1
        end
    else
        ledState(i,1) = 0; %% if led is off set column 1 digit to 0
        ledState(i,2) = -1; %% Since led is off set side to -1
    end
end

%% ======== Trial Stats ====== %%
trials = table;

trials.StartFrame = find(ledState(:,1), 1, 'first'); %%Find the start of trial one
trialCounter = 1;                               %% Initialize to trial one
firstFound = true;                              %% Start frame has been found
  
for i = (trials.StartFrame(1)+1):length(ledState)        %% Start one frame after trial 1 start and iterate to end
    if(firstFound)                              %% Search for the first instance of 0 (led OFF) after the start of a trial
        if(ledState(i,1) == 0)                  %% This will be the end of the respective trial
            trials.EndFrame(trialCounter) = i;        %% Set end frame to current index
            trialCounter = trialCounter + 1;    %% Advance trial counter
            firstFound = false;                 %% Reset firstFound so we can search for the first index of the next trial
        end
    else
        if(ledState(i,1) == 1)                  %% If we have not found the first frame of the trial, look for it
            firstFound = true;
            trials.StartFrame(trialCounter) = i;         %% Set first frame of this new trial to current index once found
        end
    end
    
    if(i == length(ledState))                   %% Edge case for if the last frame is during a trial it will set the end frame to the last index
        if(firstFound)                          %% Typically not the case, but incase of a crash or video failure 
            trials.EndFrame(trialCounter) = i;
        end
    end
end

%numTrials = length(ExperimentData.trialNumber);   
numTrials = length(trials);

for j = 1:numTrials
    trials.TrialType(j) = ledState(trials.StartFrame(j),2);      %% Set whether a specific trial is an L (0) or R (1) trial in col 3   
end

trials.odor = Odors(1:length(trials.StartFrame));

%% ======= Total Distance Calculations===== %%
noseDistance = zeros(totalFrames, 1);
bodyDistance = zeros(totalFrames, 1);


for index = 1:totalFrames-1 
%%=======Nose Distance=======%%
       pair1 = [noseCoords(index,1) noseCoords(index,2)];       %% Get the first ordered pair
       pair2 = [noseCoords(index+1,1) noseCoords(index+1,2)];   %% Get the next ordered pair in line
       coordinate = [pair1; pair2];                             %% Vertically concatenate the two pairs into a 2x2 matrix
       noseDistance(index) = pdist(coordinate);                 %% Get the euclidean distance between the two points
%%=======Body Distance=======%%
       pair1 = [bodyCoords(index,1) bodyCoords(index,2)];       %% Get the first ordered pair
       pair2 = [bodyCoords(index+1,1) bodyCoords(index+1,2)];   %% Get the next ordered pair in line
       coordinate = [pair1; pair2];                             %% Vertically concatenate the two pairs into a 2x2 matrix
       bodyDistance(index) = pdist(coordinate);                 %% Get the euclidean distance between the two points
  
end


totalDistance_nose = sum(noseDistance); %% Total distance the nose traveled 
totalDistance_body = sum(bodyDistance); %% Total distance the body traveled

%% ======= Total Distance Per Trial ======= %%

for i = 1:length(trials.StartFrame)
   tempArray = [];
   for j = trials.StartFrame(i):(trials.EndFrame(i)-1)
       pair1 = [noseCoords(j,1) noseCoords(j,2)];       %% Get the first ordered pair
       pair2 = [noseCoords(j+1,1) noseCoords(j+1,2)];   %% Get the next ordered pair in line
       coordinate = [pair1; pair2];                             %% Vertically concatenate the two pairs into a 2x2 matrix
       tempArray(j) = pdist(coordinate);
   end
       trials.trialNoseDistance(i) = sum(tempArray);
   for k = trials.StartFrame(i):(trials.EndFrame(i)-1)
       pair1 = [bodyCoords(k,1) bodyCoords(k,2)];       %% Get the first ordered pair
       pair2 = [bodyCoords(k+1,1) bodyCoords(k+1,2)];   %% Get the next ordered pair in line
       coordinate = [pair1; pair2];                             %% Vertically concatenate the two pairs into a 2x2 matrix
       bodyDistance(k) = pdist(coordinate);                 %% Get the euclidean distance between the two points
   end
    trials.trialBodyDistance(i) = sum(tempArray);
end

%% ======= Velocities (Gotta go fast boi) ==== % 


