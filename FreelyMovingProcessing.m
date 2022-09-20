%%% Freely Moving DLC Output Processing Script %%%
%%% Dewan Lab %%%
%%% Austin Pauley & Sam Caton %%%
%%% 8-16-2022 %%%

%% ===== Read Data File and Adjustments  ===== %%
clear; %%Refresh Everything


%% ======= Initial Settings ======= %%
initialSettings = inputApp;
initialSettings.Settings.Visible = 1;
uiwait(initialSettings.UIFigure);


%Import tracking data
if batchProcess == 1
    poseDataDir = uigetdir('D:\','Select Pose Data Folder:');               % If we're batch processing, get all the files in the folder
    poseDataFiles = dir(fullfile(poseDataDir,'*.csv'));
else
    [file, directory] = uigetfile('*.csv','Select Experimental File: ');    % If we're single file processing, get the single file
    poseDataFiles = dir(fullfile(directory,file));
end

experimentDataDir = uigetdir('D:\', 'Select Experimental Data Folder:');
videoDir = uigetdir('D:\', 'Select Videos Folder: ');

for z = 1:length(poseDataFiles)

clearvars -except inputROI confidenceThreshold createVideo frames_to_average...     % Global variables we don't want to clear between runs
    ledLeftX ledRightX max_bad_frames ...
    z experimentDataDir poseDataDir videoDir poseDataFiles;

poseDataPath = strcat(poseDataFiles(z).folder,'\', poseDataFiles(z).name);          % Folder with the PoseData files
fileStem = poseDataFiles(z).name(1:end-4);                                          % Get the AnimalName-Treatment prefix
experiementDataPath = strcat(experimentDataDir, '\', fileStem, '.h5');              % Get H5 File Folder
videoPath = strcat(videoDir, '\', fileStem,'.mp4');                                 % Get Video Folders


opts = detectImportOptions(poseDataPath);
opts = setvartype(opts,1,'double');

try
    PoseData = readtable(poseDataPath,opts);                                        % Try to open the PoseData file
catch                                                                               % If it fails, skip to the next file
    fprintf(2,strcat('Unable to open PoseData file! Trying next file!\n'));
    continue;
end

PoseData.Properties.VariableNames = ["index","noseX", "noseY", "noseLike", "headX", "headY",  "headLike",...
     "LHeadbarX", "LHeadbarY", "LHeadbarLike", "RHeadbarX", "RHeadbarY", "RHeadbarLike", "bodyX", "bodyY", "bodyLike",...
     "rearX","rearY","rearLike","LEDX","LEDY","LEDLike"];                           %%Resets all of the headers to make grabbing data easier; head/cannula/endoscope will all be called 'head' here for simplicity
%PoseData = PoseData((1:end),:);
prefix = fileStem(1:7);

totalFrames = length(PoseData.index);

if(PoseData.noseX(totalFrames) == 0)                                                %% Occasionally there will be a data file with a row of zeros at the end that breaks the code
    PoseData(totalFrames,:) = [];                                                   %% If this exists, delete it
    totalFrames = totalFrames - 1;
end

%Import Experimental H5

try
    ExperimentData = h5read(experiementDataPath,"/Trials");                         %% Try to import the H5 File, if it doesn't work or doesn't exist, skip to the next PoseData file
catch
    fprintf(2, strcat('Unable to open Experimental Data File: ', fileStem,'.h5! Trying next file!\n'));
    continue;
end

numTrials = length(ExperimentData.trialNumber);                                     %% Get number of Trials from Experimental Data

%Extract odors from the H5
ImportOdors = ExperimentData.odor';                                                 %% The Odor data is stored awkwardly in the H5 (sideways and with lots of spaces
Odors = cell(length(ImportOdors),1);                                                %% Transform the odor data and parse out all of the empty space

for i = 1:length(ImportOdors)
    Odors{i} = deblank(ImportOdors(i,:));
end

%Import Video
try
    video = VideoReader(videoPath);                                                 %% Try to import the video, if it doesn't work or doesn't exist, skip to the next PoseData file!
catch
    fprintf(2,strcat('Cannot open video file: ', fileStem,'.mp4! Trying next file!'));
    continue;
end

framerate = video.framerate;                                                        %% Set up some variables from the video file
timePerFrame = 1/framerate;
state = -1;

%% ======= ROI Determination ======= %%
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

%% ======= X/Y Coordinates ======= %%
noseCoords = [PoseData.noseX, PoseData.noseY];
headCoords = [PoseData.headX, PoseData.headY];
LHeadbarCoords = [PoseData.LHeadbarX, PoseData.LHeadbarY];
RHeadbarCoords = [PoseData.RHeadbarX, PoseData.RHeadbarY];
bodyCoords = [PoseData.bodyX, PoseData.bodyY];
rearCoords = [PoseData.rearX, PoseData.rearY];
LEDCoords = [PoseData.LEDX, PoseData.LEDY];

bodyCoords(:,2) = abs(headCoords(:,2)-video.Height);


%% ======== Certainty Values ======== %%

noseLike = [PoseData.noseLike];
headLike = [PoseData.headLike];
LHeadbarLike = [PoseData.LHeadbarLike];
RHeadbarLike = [PoseData.RHeadbarLike];
bodyLike = [PoseData.bodyLike];
rearLike = [PoseData.rearLike];
LEDLike = [PoseData.LEDLike];

%% ======== Bad Frames ======= %%
%Gather all of the indicies (frame number + 1) where the certainty is below our defined threshold
badNoseFrames = find(noseLike < confidenceThreshold);
badHeadFrames = find(headLike < confidenceThreshold);
badLHFrames = find (LHeadbarLike < confidenceThreshold);
badRHFrames = find(RHeadbarLike < confidenceThreshold);
badBodyFrames = find(bodyLike < confidenceThreshold);
badRearFrames = find(rearLike < confidenceThreshold);

%% ======== LED Stats ======== %%
ledState = zeros(size(LEDCoords));
for i = 1:length(LEDCoords)
    if(LEDLike(i) >= 0.95)
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

trialStats.StartFrame = find(ledState(:,1), 1, 'first');                    %% Find the start of trial one
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

if(length(trialStats.StartFrame) ~= numTrials)                              %% If the video cuts off early, use the number of trials in video as numTrials
    numTrials = min(numTrials, length(trialStats.StartFrame));
end

for j = 1:numTrials
    trialStats.TrialType(j) = ledState(trialStats.StartFrame(j),2);         %% Set whether a specific trial is an L (0) or R (1) trial in col 3   
end

%trialStats.odor = Odors(1:length(trialStats.StartFrame));
trialStats.odor = Odors(1:numTrials);

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


totalDistance_nose = sum(noseDistance(trialStats.StartFrame(1):trialStats.EndFrame(numTrials)));                                     %% Total distance the nose traveled between the first trial and last trial
totalDistance_body = sum(bodyDistance(trialStats.StartFrame(1):trialStats.EndFrame(numTrials)));                                     %% Total distance the body traveled between the first trial and last trial

%%======= % Distance Change / Time ==========%%
experimentStartFrame = trialStats.StartFrame(1);
experimentEndFrame = trialStats.EndFrame(numTrials);

for i = 0:floor((experimentEndFrame-experimentStartFrame)/framerate)-1      %% Bin data by adding distance changes over *framerate* length bins minus the endFrame
    changeSum = 0;                                                          %% May discard last few frames due to rounding totalFrames/framerate down

    for j = 1:framerate
       changeSum = changeSum + bodyDistance(int32(j+((i*framerate)+(experimentStartFrame-1))));  %% Due to rounding we ended up with a float, convert back to int for indexing
    end
    
    binnedPercentChange(i+1,1) = changeSum;                                 %% Add 1 since index i starts at 0 for multiplying i*j to get our frame
    
end

percentChange = binnedPercentChange / totalDistance_body;                   %% Divide by total distance traveled to get percent change

%% ======= Per Trial Distance  ======= %%

for i = 1:numTrials
    trialStats.trialNoseDistance(i) = sum(noseDistance(trialStats.StartFrame(i):trialStats.EndFrame(i))); %% Total Nose Distance per trial
    trialStats.trialBodyDistance(i) = sum(bodyDistance(trialStats.StartFrame(i):trialStats.EndFrame(i))); %% Total Body Distance per trial    
end

%% ======= Speeds (Gotta go fast boi) ==== %%
noseSpeed(:,1) = movmean(noseDistance/timePerFrame, frames_to_average);
bodySpeed(:,1) = movmean(bodyDistance/timePerFrame, frames_to_average);

averageNoseSpeed = mean(noseSpeed((experimentStartFrame:experimentEndFrame),1));                                    %% PX/Second
averageBodySpeed = mean(bodySpeed((experimentStartFrame:experimentEndFrame),1));                                    %% PX/Second

for i = 1:numTrials
    trialStats.averageTrialNoseSpeed(i) = mean(noseSpeed(trialStats.StartFrame(i):trialStats.EndFrame(i)));
    trialStats.maxTrialNoseSpeed(i) = max(noseSpeed(trialStats.StartFrame(i):trialStats.EndFrame(i)));
    trialStats.minTrialNoseSpeed(i) = min(noseSpeed(trialStats.StartFrame(i):trialStats.EndFrame(i)));
    trialStats.averageTrialBodySpeed(i) = mean(bodySpeed(trialStats.StartFrame(i):trialStats.EndFrame(i)));
    trialStats.maxTrialBodySpeed(i) = max(bodySpeed(trialStats.StartFrame(i):trialStats.EndFrame(i)));
    trialStats.minTrialBodySpeed(i) = min(bodySpeed(trialStats.StartFrame(i):trialStats.EndFrame(i)));
end

%% ====== Total Stats ======== %
totalStats = table;

totalStats.TotalBodyDistance = totalDistance_body;
totalStats.AverageBodySpeed = averageBodySpeed;
totalStats.MaxBodySpeed = max(bodySpeed(experimentStartFrame:experimentEndFrame));

%% ======= Time in ROI ======= %%
ROIStats = table;
centralROI = [700,475;700,625;1100,475;1100,625];
for i = 1:numTrials
%   trialStats.LROI_Nose_Time(i) = 0;
%   trialStats.RROI_Nose_Time(i) = 0;
%   trialStats.LROI_Body_Time(i) = 0;
%   trialStats.RROI_Body_Time(i) = 0;
    ROIStats.PreTimeBodyLROI(i) = 0;
    ROIStats.PreTimeBodyRROI(i) = 0;
    ROIStats.PreTimeBodyCROI(i) = 0;
    ROIStats.TimeBodyLROI(i) = 0;
    ROIStats.TimeBodyRROI(i) = 0;
    ROIStats.TimeBodyRROI(i) = 0;

    startFrame = trialStats.StartFrame(i);
    endFrame = trialStats.EndFrame(i);
 
    
%   [inL] = inpolygon(noseCoords(startFrame:endFrame,1),noseCoords(startFrame:endFrame,2),odorLROI(1,:),odorLROI(2,:));
%   trialStats.LROI_Nose_Time(i) = sum(inL)*timePerFrame;
%   [inR] =  inpolygon(noseCoords(startFrame:endFrame,1),noseCoords(startFrame:endFrame,2),odorRROI(1,:),odorRROI(2,:));
%   trialStats.RROI_Nose_Time(i) = sum(inR)*timePerFrame;
    [inL] = inpolygon(bodyCoords(startFrame:endFrame,1),bodyCoords(startFrame:endFrame,2),odorLROI(:,1),odorLROI(:,2));         
    ROIStats.TimeBodyLROI(i) = sum(inL) * timePerFrame;
    [inR] =  inpolygon(bodyCoords(startFrame:endFrame,1),bodyCoords(startFrame:endFrame,2),odorRROI(:,1),odorRROI(:,2));
    ROIStats.TimeBodyRROI(i) = sum(inR) * timePerFrame;
    [inC] = inpolygon(bodyCoords(startFrame:endFrame, 1), bodyCoords(startFrame:endFrame,2),centralROI(:,1),centralROI(:,2));
    ROIStats.TimeBodyCROI(i) = sum(inC) * timePerFrame;
    if i ~= 1                                                               % Trial one does not always have enough time beforehand; either ignore or find workaround
        preTimeStartFrame = startFrame - (framerate*10-1);
        preTimeEndFrame = startFrame - 1;
        [inL] = inpolygon(bodyCoords(preTimeStartFrame:preTimeEndFrame,1),bodyCoords(preTimeStartFrame:preTimeEndFrame,2),odorLROI(:,1),odorLROI(:,2));
        ROIStats.PreTimeBodyLROI(i) = sum(inL) * timePerFrame;
        [inR] = inpolygon(bodyCoords(preTimeStartFrame:preTimeEndFrame,1),bodyCoords(preTimeStartFrame:preTimeEndFrame,2),odorRROI(:,1),odorRROI(:,2));
        ROIStats.PreTimeBodyRROI(i) = sum(inR) * timePerFrame;
        [inC] = inpolygon(bodyCoords(preTimeStartFrame:preTimeEndFrame, 1), bodyCoords(preTimeStartFrame:preTimeEndFrame,2),centralROI(:,1),centralROI(:,2));
        ROIStats.PreTimeBodyCROI(i) = sum(inC) * timePerFrame;
    end
end

totalStats.TotalPreTimeBodyLROI(1) = sum(ROIStats.PreTimeBodyLROI);
totalStats.TotalTimeBodyLROI(1) = sum(ROIStats.TimeBodyLROI);
totalStats.TotalPreTimeBodyRROI(1) = sum(ROIStats.PreTimeBodyRROI);
totalStats.TotalTimeBodyRROI(1) = sum(ROIStats.TimeBodyRROI);
totalStats.TotalPreTimeBodyCROI(1) = sum(ROIStats.PreTimeBodyCROI);
totalStats.TotalTimeBodyCROI(1) = sum(ROIStats.TimeBodyCROI);

%% ====== Aversion Index ====== %%
% Defined as preRoiTime - ROITime for the active ROI during an odor presentation
ROIStats.trialType(1) = trialStats.TrialType(1);                            % Sets trial type for trial one since its skipped in the loop
ROIStats.aversionIndex(1) = 0;
ROIStats.odorCheck(1) = 0;                                                  % 1 if animal enters ROI with odor, 0 if animal does not enter ROI with odor
for trial = 2:numTrials                                                     % Ignore trial one due to no pretrial ROI time
    if trialStats.TrialType(trial) == 0
        ROIStats.trialType(trial) = 0;
        if(ROIStats.TimeBodyLROI(trial) ~= 0)
            ROIStats.odorCheck(trial) = 1; 
            ROIStats.aversionIndex(trial) =  ROIStats.TimeBodyLROI(trial) - ROIStats.PreTimeBodyLROI(trial);
        end     
    elseif trialStats.TrialType(trial) == 1
        ROIStats.trialType(trial) = 1;
        if(ROIStats.TimeBodyRROI(trial) ~= 0)
            ROIStats.odorCheck(trial) = 1;
            ROIStats.aversionIndex(trial) = ROIStats.TimeBodyLROI(trial) - ROIStats.PreTimeBodyLROI(trial);
        end
    end
end


%% ====== Per-Odor ROI Stats ====== %%
individualOdors = unique(Odors);

OdorStats = table();

for i = 1:length(individualOdors)
    Odor = cell2mat(individualOdors(i));
    col1 = strcat("Odor_", Odor,"_pre");
    col2 = strcat("Odor_", Odor, "_post");
    col3 = strcat("Odor_", Odor, "_AI");
    OdorStats.(col1) = zeros(numTrials,1);
    OdorStats.(col2) = zeros(numTrials,1);
    odorTrials = find(strcmp(trialStats.odor, Odor));
    
    for trials = 1:length(odorTrials)
        currentTrial = odorTrials(trials);
        if(trialStats.TrialType(currentTrial) == 0)
            OdorStats.(col1)(currentTrial) = ROIStats.PreTimeBodyLROI(currentTrial);
            OdorStats.(col2)(currentTrial) = ROIStats.TimeBodyLROI(currentTrial);
        else
            OdorStats.(col1)(currentTrial) = ROIStats.PreTimeBodyRROI(currentTrial);
            OdorStats.(col2)(currentTrial) = ROIStats.TimeBodyRROI(currentTrial);
        end
    end

    totalStats.(col3) = sum(OdorStats.(col2)) - sum(OdorStats.(col1))  ;       % Per Odor Aversion Index
end






%% ======= EXPORT DATA ======= %
path = strcat('.\\Output\', prefix);
mkdir(path);
cd(path);
writetable(totalStats,strcat(prefix,'-totalStats.xlsx'), 'Sheet','Data');
writetable(table(percentChange),strcat(prefix,'-PercentChangePos.xlsx'),'Sheet','Data','WriteVariableNames',0);
writetable(ROIStats,strcat(prefix,'-ROIStats.xlsx'), 'Sheet','Data');
writetable(OdorStats,strcat(prefix,'-OdorStats.xlsx'), 'Sheet','Data');
%Positional Heatmap
colormap('hot');
histogram2(bodyCoords(:,1),bodyCoords(:,2),[(floor(max(bodyCoords(:,1)))),floor(max(bodyCoords(:,2)))],'DisplayStyle','tile')
savefig(strcat(prefix,'-PosHeatmap.fig'));
hold on;
colormap('parula')
imagesc(readFrame(video));
histogram2(bodyCoords(:,1),bodyCoords(:,2),[(floor(max(bodyCoords(:,1)))),floor(max(bodyCoords(:,2)))],'DisplayStyle','tile')
savefig(strcat(prefix,'-PosHeatmap-Img.fig'));
hold off;


%Percent Change Plot
plot(1:floor((experimentEndFrame-experimentStartFrame)/framerate), percentChange);    %% Graph these relative to our frames
title('% Change in Distance over Time');
xlabel('Time (s)');
ylabel('Percent total movement');
graph = gcf;
exportgraphics(graph,strcat(prefix,'-percentChange.tiff'));
cd ../..

end