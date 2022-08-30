%%% Freely Moving DLC Output Processing Script %%%
%%% Dewan Lab %%%
%%% Austin Pauley & Sam Caton %%%
%%% 8-16-2022 %%%

%% ===== Read Data File and Adjustments  ===== %%
clear; %%Refresh Everything

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

[file, path] = uigetfile('*.mp4', 'Select Eperiment Video:', 'D:\');

video = VideoReader(strcat(path,file));

framerate = video.framerate;
timePerFrame = 1/framerate;
totalFrames = length(PoseData.index);
state = -1;
%% ======= Initial Settings ==== %%
initialSettings = inputApp;
initialSettings.Settings.Visible = 1;
uiwait(initialSettings.UIFigure);

%% ======= ROI Determination ==== %%
odorLROI = zeros(4,2); 
odorRROI = zeros(4,2);

L_Port1 = zeros(1,2);
L_Port3 = zeros(1,2);
R_Port1 = zeros(1,2);
R_Port3 = zeros(1,2);
       
%%Input ROI if needed
if inputROI == 1
        state = 1;
        inputGUI = inputApp;
        inputGUI.ROI.Visible = 1;
        inputGUI.Settings.Visible = 0;
        uiwait(inputGUI.UIFigure);
end

initFrame = readFrame(video);
%%Draw ROIs
while state ~= 1
    display = imshow(initFrame);
    
%     msg = msgbox('Please label the three ports for the RIGHT bank', 'Bank 1', 'modal');
%     uiwait(msg);
%     B1P1 = drawpoint('Color', 'Cyan');
%     B1P2 = drawpoint('Color', 'Green');
%     B1P3 = drawpoint('Color', 'Red');
%     
%     msg = msgbox('Please label the three ports for the LEFT bank', 'Bank 2','modal');
%     uiwait(msg);
%     B2P1 = drawpoint('Color', 'Cyan');
%     B2P2 = drawpoint('Color', 'Green');
%     B3P3 = drawpoint('Color', 'Red');
    
    msg = msgbox('Please outline the ROI for ODOR on the Right', 'Odor','modal');
    uiwait(msg);
    odorRight = drawrectangle('LineWidth', 7, 'Color', 'Yellow');
    
    msg = msgbox('Please outline the ROI for ODOR on the Left', 'Odor','modal');
    uiwait(msg);
    odorLeft = drawrectangle('LineWidth', 7, 'Color', 'Magenta');
    
    confirmation = questdlg('Would you like to redraw the ROIs?' ,'Finished?', 'No, Process Video', 'Yes', 'No, Process Video');

    switch confirmation
        case 'No, Process Video'
            state = 1;  
%             bank1Locs(1,:) = B1P1.Position;
%             bank1Locs(2,:) = B1P2.Position;
%             bank1Locs(3,:) = B1P3.Position;
%             bank2Locs(1,:) = B2P1.Position;
%             bank2Locs(2,:) = B2P2.Position;
%             bank2Locs(3,:) = B2P2.Position;

            odorLROI = odorLeft.Vertices;
            odorRROI = odorRight.Vertices;
            
%             L_Port1 = bank1Locs(1,:);
%             L_Port3 = bank1Locs(3,:);
%             R_Port1 = bank2Locs(1,:);
%             R_Port3 = bank2Locs(3,:);
            
            close(display.Parent.Parent);
        case 'Yes, Restart'
            state = 0;
            close(display.Parent.Parent);
    end
end

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
        ledState(i,1) = 1;                                                  %% If led is on set the column 1 digit to 1
        if(LEDCoords(i,1) <= ledLeftX)
            ledState(i,2) = 0;                                              %% If led is on the left, set column 2 digit to 0
        elseif(LEDCoords(i,1) >= ledRightX)
            ledState(i,2) = 1;                                              %% If led is on the right, set column 2 digit to 1
        end
    else
        ledState(i,1) = 0;                                                  %% if led is off set column 1 digit to 0
        ledState(i,2) = -1;                                                 %% Since led is off set side to -1
    end
end

%% ======== Trial Stats ====== %%
trialStats = table;

trialStats.StartFrame = find(ledState(:,1), 1, 'first');                    %%Find the start of trial one
trialCounter = 1;                                                           %% Initialize to trial one
firstFound = true;                                                          %% Start frame has been found
  
for i = (trialStats.StartFrame(1)+1):length(ledState)                       %% Start one frame after trial 1 start and iterate to end
    if(firstFound)                                                          %% Search for the first instance of 0 (led OFF) after the start of a trial
        if(ledState(i,1) == 0)                                              %% This will be the end of the respective trial
            trialStats.EndFrame(trialCounter) = i;                          %% Set end frame to current index
            trialCounter = trialCounter + 1;                                %% Advance trial counter
            firstFound = false;                                             %% Reset firstFound so we can search for the first index of the next trial
        end
    else
        if(ledState(i,1) == 1)                                              %% If we have not found the first frame of the trial, look for it
            firstFound = true;
            trialStats.StartFrame(trialCounter) = i;                        %% Set first frame of this new trial to current index once found
        end
    end
    
    if(i == length(ledState))                                               %% Edge case for if the last frame is during a trial it will set the end frame to the last index
        if(firstFound)                                                      %% Typically not the case, but incase of a crash or video failure 
            trialStats.EndFrame(trialCounter) = i;
        end
    end
end

numTrials = length(ExperimentData.trialNumber);                             %%This is the real way to do it just does not work with the test data
%numTrials = length(trialStats.StartFrame);

for j = 1:numTrials
    trialStats.TrialType(j) = ledState(trialStats.StartFrame(j),2);         %% Set whether a specific trial is an L (0) or R (1) trial in col 3   
end

trialStats.odor = Odors(1:length(trialStats.StartFrame));

%% ======= Total Distance Calculations===== %%
noseDistance = zeros(totalFrames, 1);
bodyDistance = zeros(totalFrames, 1);


for index = 1:totalFrames-1 
%%=======Nose Distance=======%%
       pair1 = [noseCoords(index,1) noseCoords(index,2)];                   %% Get the first ordered pair
       pair2 = [noseCoords(index+1,1) noseCoords(index+1,2)];               %% Get the next ordered pair in line
       coordinate = [pair1; pair2];                                         %% Vertically concatenate the two pairs into a 2x2 matrix
       noseDistance(index) = pdist(coordinate);                             %% Get the euclidean distance between the two points
%%=======Body Distance=======%%
       pair1 = [bodyCoords(index,1) bodyCoords(index,2)];                   %% Get the first ordered pair
       pair2 = [bodyCoords(index+1,1) bodyCoords(index+1,2)];               %% Get the next ordered pair in line
       coordinate = [pair1; pair2];                                         %% Vertically concatenate the two pairs into a 2x2 matrix
       bodyDistance(index) = pdist(coordinate);                             %% Get the euclidean distance between the two points 
end


totalDistance_nose = sum(noseDistance);                                     %% Total distance the nose traveled 
totalDistance_body = sum(bodyDistance);                                     %% Total distance the body traveled

%%======= % Distance Change / Time ==========%%
%percentChange = (bodyDistance / totalDistance_body) *100;                  %% Divide our small changes by total change to get % change

for i = 0:floor(totalFrames/framerate)-1                                    %% Bin data by adding distance changes over *framerate* length bins
    changeSum = 0;                                                          %% disregards last few frames due to rounding totalFrames/framerate down

    for j = 1:framerate
       changeSum = changeSum + bodyDistance(int32(j+(i*framerate)));
    end
    
    binnedPercentChange(i+1,1) = changeSum;
    
end

percentChange = binnedPercentChange / totalDistance_body;

plot(1:floor(totalFrames/framerate), percentChange);                        %% Graph these relative to our frames
title('% Change in Distance over Time');
xlabel('Time (s)');
ylabel('Percent total movement');



%% ======= Per Trial Distance  ======= %%

for i = 1:length(trialStats.StartFrame)

    trialStats.trialNoseDistance(i) = sum(noseDistance(trialStats.StartFrame(i):trialStats.EndFrame(i))); %% Total Nose Distance per trial
    trialStats.trialBodyDistance(i) = sum(bodyDistance(trialStats.StartFrame(i):trialStats.EndFrame(i))); %% Total Body Distance per trial
    
end

%% ======= Speeds (Gotta go fast boi) ==== %%
noseSpeed(:,1) = movmean(noseDistance/timePerFrame, frames_to_average);
bodySpeed(:,1) = movmean(bodyDistance/timePerFrame, frames_to_average);

averageNoseSpeed = mean(noseSpeed(:,1));                                    %% PX/Second
averageBodySpeed = mean(bodySpeed(:,1));                                    %% PX/Second

for i = 1:length(trialStats.StartFrame)
    trialStats.averageTrialNoseSpeed(i) = mean(noseSpeed(trialStats.StartFrame(i):trialStats.EndFrame(i)));
    trialStats.maxTrialNoseSpeed(i) = max(noseSpeed(trialStats.StartFrame(i):trialStats.EndFrame(i)));
    trialStats.minTrialNoseSpeed(i) = min(noseSpeed(trialStats.StartFrame(i):trialStats.EndFrame(i)));
    trialStats.averageTrialBodySpeed(i) = mean(bodySpeed(trialStats.StartFrame(i):trialStats.EndFrame(i)));
    trialStats.maxTrialBodySpeed(i) = max(bodySpeed(trialStats.StartFrame(i):trialStats.EndFrame(i)));
    trialStats.minTrialBodySpeed(i) = min(bodySpeed(trialStats.StartFrame(i):trialStats.EndFrame(i)));
end



%%======= Speed / Time ==========%
plot(1:totalFrames, bodySpeed);                                             %% Graph these relative to our frames
title('Body Speed over Time');
xlabel('Time (frames)');
ylabel('Body Speed (px/s)');

%% ======= Time in ROI ======= %%
for i = 1:numTrials
    trialStats.LROI_Nose_Time(i) = 0;
    trialStats.RROI_Nose_Time(i) = 0;
    trialStats.LROI_Body_Time(i) = 0;
    trialStats.RROI_Body_Time(i) = 0;
    startFrame = trialStats.StartFrame(i);
    endFrame = trialStats.EndFrame(i);
    
    [inL] = inpolygon(noseCoords(startFrame:endFrame,1),noseCoords(startFrame:endFrame,2),odorLROI(1,:),odorLROI(2,:));
    trialStats.LROI_Nose_Time(i) = sum(inL)*timePerFrame;
    [inR] =  inpolygon(noseCoords(startFrame:endFrame,1),noseCoords(startFrame:endFrame,2),odorRROI(1,:),odorRROI(2,:));
    trialStats.RROI_Nose_Time(i) = sum(inR)*timePerFrame;
    [inL] = inpolygon(bodyCoords(startFrame:endFrame,1),bodyCoords(startFrame:endFrame,2),odorLROI(1,:),odorLROI(2,:));
    trialStats.LROI_Body_Time(i) = sum(inL)*timePerFrame;
    [inR] =  inpolygon(bodyCoords(startFrame:endFrame,1),bodyCoords(startFrame:endFrame,2),odorRROI(1,:),odorRROI(2,:));
    trialStats.RROI_Body_Time(i) = sum(inR)*timePerFrame;

end