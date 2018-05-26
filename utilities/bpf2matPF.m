function retVal = bpf2matPF(filename)
%converts BPF into Matlab file
% (C) Dino Dvorak 2015
% dino@indus3.net

doTetrode = 1;
doPosition = 1;
doEEG = 0;
doSync = 0;
doInput = 0;

D = dir(filename);
sizeoffile=D.bytes;

fid = fopen(filename, 'r');

%record lengts in bytes
eegRecordLength = 0;
syncRecordLength = 0;
singleRecordLength = 0;
stereoRecordLength = 0;
tetrodeRecordLength = 0;
keyRecordLength = 0;
inputRecordLength = 0;
outputRecordLength = 0;
roomRecordLength = 0;
arenaRecordLength = 0;

%EEG
eegChannelsPerBlock = 0;
eegSamplesPerBlock = 0;
eegBlocks = 0;
eegFS = 0;

%tetrode (multiple samples per block of data)
tetrodeChannelsPerBlock = 0;
tetrodeSamplesPerBlock = 0;
tetrodeBlocks = 0;
tetrodeFS = 0; %sampling frequency

%sync (multiple samples per block of data)
syncSamplesPerBlock = 0;
syncBlocks = 0;

%length of the header in bytes
headerLength = 0;
headerLen = 0;

%input
inputBlocks = 0;

%output
outputBlocks = 0;

%room arena
roomBlocks = 0;
arenaBlocks = 0;

bytesRead = 0; %byte counter

%% Read Header
for i = 1:500  %first 500 lines of file, break at header
   tline = fgetl(fid);
   
   %add length of the line to the counter
   headerLength = headerLength + length(tline);
   
   %EegSamplingFrequency.0 (2000)
   stringToFind = '%EegSamplingFrequency';
   if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        eegFS = str2double(tline(26:length(tline)-1));
        continue;
   end
   
   %SpikeSamplingFrequency.0 (2000)
   stringToFind = '%SpikeSamplingFrequency';
   if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        tetrodeFS = str2double(tline(28:length(tline)-1));
        continue;
   end
  
   %EegRecordFormat.10005 Identifier.1 100usTimeStamp.4 Voltage.2[5][1000]
   stringToFind = '%EegRecordFormat';
   if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %parts = strsplit(' ', tline, 'omit'); 
        parts = strsplit(tline,' ');
        part = parts{1};
        part = part(18:end);
        eegRecordLength = str2double(part);
        
        part = parts{4}; %Voltage.2[5][1000]
        part = part(11:end-1); % 5][1000
        %parts = strsplit(']',part,'omit');
        parts = strsplit(part,']');
        eegChannelsPerBlock = str2double(parts{1});
        part = parts{2};
        eegSamplesPerBlock = str2double(part(2:end));
        continue;
   end
   
   %SyncRecordFormat.2005 Identifier.1 100usTimeStamp.4 Voltage.2[1000]
   stringToFind = '%SyncRecordFormat';
   if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %parts = strsplit(' ', tline, 'omit'); 
        parts = strsplit(tline,' ');
        part = parts{1};
        part = part(19:end);
        syncRecordLength = str2double(part);
        
        part = parts{4};
        %parts = strsplit('[',part,'omit');
        parts = strsplit(part,'[');
        part = parts{2};
        part = part(1:end-1);
        syncSamplesPerBlock = str2double(part);
        continue;
   end
  
   %SingleRecordFormat.71 Identifier.1 100usTimeStamp.4 Channel.1 
   stringToFind = '%SingleRecordFormat';
   if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %parts = strsplit(' ', tline, 'omit');  
        parts = strsplit(tline,' ');
        part = parts{1};
        part = part(21:end);
        singleRecordLength = str2double(part);
        continue;
   end
   
   %StereoRecordFormat.135 Identifier.1 100usTimeStamp.4 StereoChannel.1
   stringToFind = '%StereoRecordFormat';
   if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %parts = strsplit(' ', tline, 'omit');  
        parts = strsplit(tline,' ');
        part = parts{1};
        part = part(21:end);
        stereoRecordLength = str2double(part);
        continue;
   end
   
   %TetrodeRecordFormat.775 Identifier.1 100usTimeStamp.4 TetrodeChannel.1
   %DiscriminatedUnit.1 Voltage.2[4][96]
   stringToFind = '%TetrodeRecordFormat';
   if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %parts = strsplit(' ', tline, 'omit'); 
        parts = strsplit(tline,' ');
        part = parts{1};
        part = part(22:end);
        tetrodeRecordLength = str2double(part);
                
        part = parts{6}; %Voltage.2[4][96]
        part = part(11:end-1); % 4][96
        %parts = strsplit(']',part,'omit');
        parts = strsplit(part,']');
        tetrodeChannelsPerBlock = str2double(parts{1});
        part = parts{2};
        tetrodeSamplesPerBlock = str2double(part(2:end));
        continue;
   end
   
   %KeyEventRecordFormat.6 Identifier.1 100usTimestamp.4 ASCII.1
   stringToFind = '%KeyEventRecordFormat';
   if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %parts = strsplit(' ', tline, 'omit'); 
        parts = strsplit(tline,' ');
        part = parts{1};
        part = part(23:end);
        keyRecordLength = str2double(part);
        continue;
   end
   
  %InputEventIdentifier.1 ('I')
  %InputEventRecordFormat.7 Identifier.1 100usTimeStamp.4 Value.2
  stringToFind = '%InputEventRecordFormat';
  if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %parts = strsplit(' ', tline, 'omit'); 
        parts = strsplit(tline,' ');
        part = parts{1};
        part = part(25:end);
        inputRecordLength = str2double(part);
        continue;
  end   
  
  %OutputEventIdentifier.1 ('O')
  %OutputEventRecordFormat.7 Identifier.1 100usTimeStamp.4 Value.2
  stringToFind = '%OutputEventRecordFormat';
  if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %parts = strsplit(' ', tline, 'omit'); 
        parts = strsplit(tline,' ');
        part = parts{1};
        part = part(26:end);
        outputRecordLength = str2double(part);
        continue;
  end   

  %RoomPositionIdentifier.1 ('R')
  %RoomPositionRecordFormat.9 Identifier.1 100usTimestamp.4 RoomX1.1 RoomY1.1 RoomAng.2
  stringToFind = '%RoomPositionRecordFormat';
  if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %parts = strsplit(' ', tline, 'omit'); 
        parts = strsplit(tline,' ');
        part = parts{1};
        part = part(27:end);
        roomRecordLength = str2double(part);
        continue;
  end   
  
 %ArenaPositionIdentifier.1 ('A')
 %ArenaPositionRecordFormat.9 Identifier.1 100usTimestamp.4 ArenaX1.1 ArenaY1.1 ArenaAng.2 
  stringToFind = '%ArenaPositionRecordFormat';
  if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %parts = strsplit(' ', tline, 'omit'); 
        parts = strsplit(tline,' ');
        part = parts{1};
        part = part(28:end);
        arenaRecordLength = str2double(part);
        continue;
  end  
  
  %ListOfEegChannels.0 ( 1 2 3 4 5 6 7 8 28)
  stringToFind = '%ListOfEegChannels';
  if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %get nubers between ( and )
        stx = strfind(tline,'(');
        edx = strfind(tline,')');
        stx = stx + 2;
        edx = edx - 1; %also exclude space
        parts = tline(stx:edx);
        %parts = strsplit(' ', parts, 'omit'); 
        parts = strsplit(parts,' ');

        eegChannels = zeros(1,length(parts));

        for gI = 1:length(parts);
            eegChannels(gI) = str2double(parts{gI});
        end
        continue;
  end
  
  %ListOfGains.0 (8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 )
  stringToFind = '%ListOfGains';
   if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %get nubers between ( and )
        stx = strfind(tline,'(');
        edx = strfind(tline,')');
        stx = stx + 1;
        edx = edx - 2; %also exclude space
        parts = tline(stx:edx);
        %parts = strsplit(' ', parts, 'omit'); 
        parts = strsplit(parts,' ');
        
        gains = zeros(1,length(parts));
        
        for gI = 1:length(parts);
            gains(gI) = str2double(parts{gI});
        end
        continue;
   end
   
  %ListOfCalibrations.0 ( 500.000000 500.000000 500.000000 500.000000 500.000000 500.000000 500.000000 500.000000 1.0)
  stringToFind = '%ListOfCalibrations';
   if length(tline) > length(stringToFind) && strcmp(tline(1:length(stringToFind)),stringToFind)
        %get nubers between ( and )
        stx = strfind(tline,'(');
        edx = strfind(tline,')');
        stx = stx + 2;
        edx = edx - 2; %also exclude space
        parts = tline(stx:edx);
        %parts = strsplit(' ', parts, 'omit'); 
        parts = strsplit(parts,' ');
        
        calibrations = zeros(1,length(parts));
        
        for gI = 1:length(parts);
            calibrations(gI) = str2double(parts{gI});
        end
        continue;
   end
  
  %12345678901234567890123456789012345678901234567890
   %%END_HEADER
   stringToFind = '%%END_HEADER';
   if strcmp(tline,stringToFind)
       headerLength = headerLength + 2*i; %need to add newline characters
       headerLen = i;
       break;
   end
   
   %matches = findstr(tline, literal);
   %fprintf(1,'%d:%s\n',num,tline);

end

disp('header ok...');

%% Get sample count
%back to beginning
frewind(fid);

%skip header
fread(fid,headerLength,'uchar');

bytesRead = bytesRead + headerLength;

%get gains for EEG channels
gainsEEG = gains(eegChannels);

%get data length (reads file to the end)
while bytesRead < sizeoffile
    
    ch = fread(fid,1,'uint8=>char');
    
    if (ch == 'E')%eng
        fread(fid,eegRecordLength-1,'uchar');
        bytesRead = bytesRead + eegRecordLength;
        eegBlocks = eegBlocks + 1;
    elseif (ch == 'Y')%sync
        fread(fid,syncRecordLength-1,'uchar');
        bytesRead = bytesRead + syncRecordLength;
        syncBlocks = syncBlocks + 1;
    elseif (ch == 'U')%single
        fread(fid,singleRecordLength-1,'uchar');
        bytesRead = bytesRead + singleRecordLength;
    elseif (ch == 'S')%stereo
        fread(fid,stereoRecordLength-1,'uchar');
        bytesRead = bytesRead + stereoRecordLength;
    elseif (ch == 'T')%tetrode
        fread(fid,tetrodeRecordLength-1,'uchar');
        bytesRead = bytesRead + tetrodeRecordLength;
        tetrodeBlocks = tetrodeBlocks + 1;
    elseif (ch == 'K')%key
        fread(fid,keyRecordLength-1,'uchar');
        bytesRead = bytesRead + keyRecordLength;
    elseif (ch == 'I')%position
        fread(fid,inputRecordLength-1,'uchar');
        bytesRead = bytesRead + inputRecordLength;
        inputBlocks = inputBlocks + 1;
    elseif (ch == 'O')%position
        fread(fid,outputRecordLength-1,'uchar');
        bytesRead = bytesRead + outputRecordLength;
        outputBlocks = outputBlocks + 1;
    elseif (ch == 'R')%room
        fread(fid,roomRecordLength-1,'uchar');
        bytesRead = bytesRead + roomRecordLength;
        roomBlocks = roomBlocks + 1;
    elseif (ch == 'A')%arena
        fread(fid,arenaRecordLength-1,'uchar');
        bytesRead = bytesRead + arenaRecordLength;
        arenaBlocks = arenaBlocks + 1;
    else 
        disp(['Unknown record format ' ch ]);
        fclose(fid);
        retVal = 0;
        return;
    end
end

disp('sample count ok...');

%% read data
if doEEG == 1
    eegData = zeros(eegChannelsPerBlock,eegSamplesPerBlock*eegBlocks,'single');
    eegTimeStamp = zeros(eegBlocks,1,'uint32');
end

%these arrays hold the data
if doTetrode == 1
    tetrodeTimestamp = zeros(tetrodeBlocks,1,'uint32');
    tetrodeData = zeros(tetrodeChannelsPerBlock,tetrodeSamplesPerBlock*tetrodeBlocks, 'single');
    tetrodeChannel = zeros(tetrodeBlocks,1,'uint8');
    tetrodeUnit = zeros(tetrodeBlocks,1,'uint8');
end

if doSync == 1
    syncData = zeros(syncSamplesPerBlock * syncBlocks,1);
    syncTimeStamp = zeros(syncBlocks,1);
end

if doInput == 1
    inputTimeStamps = zeros(inputBlocks,1,'uint32');
    inputData = zeros(inputBlocks,1);
    
    outputTimeStamps = zeros(outputBlocks,1,'uint32');
    outputData = zeros(outputBlocks,1);
end

if doPosition == 1
    roomTimeStamps = zeros(roomBlocks,1,'uint32');
    roomXY = zeros(roomBlocks,2,'uint8');
    %roomAngle = zeros(roomBlocks,1,'single');

    %arenaTimeStamps = zeros(arenaBlocks,1,'uint32');
    %arenaXY = zeros(arenaBlocks,2,'uint8');
    %arenaAngle = zeros(arenaBlocks,1,'single');
end

%back to beginning
frewind(fid);

%skip header
fread(fid,headerLength,'uchar');
bytesRead = 0;
bytesRead = bytesRead + headerLength;

%sample counters
sampleCounterTetrode = 1;       %counts tetrode samples
sampleCounerTetrodeBlocks = 1;  %counts tetrode blocks (for timestamp...)
sampleCounterSync = 1;
sampleCounterSyncBlocks = 1;
sampleCounterEeg = 1;
sampleCounterEegBlocks = 1;
sampleCounterInput = 1;
sampleCounterOutput = 1;
sampleCounterRoom = 1;
sampleCounterArena = 1;

while bytesRead < sizeoffile
    
    ch = fread(fid,1,'uint8=>char');
    
    if (ch == 'E')%eng
        if doEEG == 1
            %read 4 bytes timestamp
            eegTimeStamp(sampleCounterEegBlocks) = fread(fid,1,'uint32');

            %read data
            for chI = 1:eegChannelsPerBlock
                
                sample = fread(fid,eegSamplesPerBlock,'int16');
                sample = single(sample);
                
                %sample = single(sample) / 15.753365; //orig
                
                sample = sample * 1.5; %1.5V is the maximal range of ADC
                sample = sample / 32768; %make it -1:1
                sample = sample / gainsEEG(chI); %divide by gain
                
                eegData(chI,sampleCounterEeg:sampleCounterEeg+eegSamplesPerBlock-1) = sample;
                
                %alternative way - reading individual bytes
                %byteArray = fread(fid,2,'uchar');
                %byteArray = uint8(byteArray);
                %value = typecast(byteArray, 'int16');
                
            end

            sampleCounterEeg = sampleCounterEeg + eegSamplesPerBlock;
            sampleCounterEegBlocks = sampleCounterEegBlocks + 1;
        else
            fread(fid,eegRecordLength-1,'uchar');
        end
        bytesRead = bytesRead + eegRecordLength;
    elseif (ch == 'Y')%sync
        if doSync == 1
            %read 4 bytes timestamp
            syncTimeStamp(sampleCounterSyncBlocks) = fread(fid,1,'uint32');

            %read data
            sample = fread(fid,syncSamplesPerBlock,'int16');
            sample = single(sample);
            sample = sample * 1.5; %1.5V is the maximal range of ADC
            sample = sample / 32768; %make it -1:1
            %sample = sample / gainsEEG(chI); %divide by gain
            syncData(sampleCounterSync:sampleCounterSync+syncSamplesPerBlock-1) = sample;

            sampleCounterSync = sampleCounterSync + syncSamplesPerBlock;
            sampleCounterSyncBlocks = sampleCounterSyncBlocks + 1;
        else
            fread(fid,syncRecordLength-1,'uchar');
        end
        bytesRead = bytesRead + syncRecordLength;
    elseif (ch == 'U')%single
        fread(fid,singleRecordLength-1,'uchar');
        bytesRead = bytesRead + singleRecordLength;
    elseif (ch == 'S')%stereo
        fread(fid,stereoRecordLength-1,'uchar');
        bytesRead = bytesRead + stereoRecordLength;
    elseif (ch == 'T')%tetrode
        if doTetrode == 1
            %read 4 bytes timestamp
            tetrodeTimestamp(sampleCounerTetrodeBlocks) = fread(fid,1,'uint32');

            %read 1 byte tetrode channel
            tch = fread(fid,1,'uint8');
            tetrodeChannel(sampleCounerTetrodeBlocks) = tch;

            %read discriminated unit
            tetrodeUnit(sampleCounerTetrodeBlocks) = fread(fid,1,'uint8');

            %read data
            for chI = 1:tetrodeChannelsPerBlock
                sample = fread(fid,tetrodeSamplesPerBlock,'int16');
                sampleS = single(sample);
                sampleS = sampleS / 32768; %make it -1:1
                sampleS = sampleS * 1.5;
                
                %get actual channel
                chx = (tch*4) + chI;
                gain = gains(chx);
                sampleS = sampleS / gain;
                
                tetrodeData(chI,sampleCounterTetrode:sampleCounterTetrode+tetrodeSamplesPerBlock-1) = sampleS;
            end

            sampleCounterTetrode = sampleCounterTetrode + tetrodeSamplesPerBlock;
            sampleCounerTetrodeBlocks = sampleCounerTetrodeBlocks + 1;
        else
                fread(fid,tetrodeRecordLength-1,'uchar');
        end
        bytesRead = bytesRead + tetrodeRecordLength;
    elseif (ch == 'K')%key
        fread(fid,keyRecordLength-1,'uchar');
        bytesRead = bytesRead + keyRecordLength;
    elseif (ch == 'I')%input
        if doInput == 1
            %read 4 bytes timestamp
            inputTimeStamps(sampleCounterInput) =  fread(fid,1,'uint32');
            %read 2 bytes value
            inputData(sampleCounterInput) = fread(fid,1,'int16');

            sampleCounterInput = sampleCounterInput + 1;
        else
            fread(fid,inputRecordLength-1,'uchar');
        end
        bytesRead = bytesRead + inputRecordLength;
    elseif (ch == 'O')%output
        if doInput == 1
            %read 4 bytes timestamp
            outputTimeStamps(sampleCounterOutput) =  fread(fid,1,'uint32');
            %read 2 bytes value
            outputData(sampleCounterOutput) = fread(fid,1,'int16');

            sampleCounterOutput = sampleCounterOutput + 1;
        else
            fread(fid,outputRecordLength-1,'uchar');
        end
        bytesRead = bytesRead + outputRecordLength;
    elseif (ch == 'R')%room
         if doPosition == 1
            %RoomPositionRecordFormat.14 
            %Identifier.1 
            %100usTimestamp.4 
            %RoomX1.1 RoomY1.1 RoomZ1.1 RoomAng.2 
            %RoomAngX.2 RoomAngY.2
             
            %read 4 bytes timestamp
            roomTimeStamps(sampleCounterRoom) =  fread(fid,1,'uint32');

            %read 4 bytes X Y
            temp = fread(fid,roomRecordLength-5,'uint8');
            
            %xyz
            roomXY(sampleCounterRoom,1) = temp(1); %X
            roomXY(sampleCounterRoom,2) = temp(2); %Y
            %roomXY(sampleCounterRoom,2) = temp(3); %Z

            %angle
            %sample = typecast([xyaa(3) xyaa(4)],'int16');
            %sampleS = single(sample) / 15.753365;
            %positionAngle(sampleCounterPosition) = sampleS;

            sampleCounterRoom = sampleCounterRoom + 1;
        else
            fread(fid,roomRecordLength-1,'uchar');
        end
        bytesRead = bytesRead + roomRecordLength;
    elseif (ch == 'A')%arena
       if doPosition == 2
            %read 4 bytes timestamp
            arenaTimeStamps(sampleCounterArena) =  fread(fid,1,'uint32');

            temp = fread(fid,9,'uint8');
            
            %xyz
            arenaXY(sampleCounterArena,1) = temp(1); %X
            arenaXY(sampleCounterArena,2) = temp(2); %Y
            %roomXY(sampleCounterRoom,2) = temp(3); %Z
          
            sampleCounterArena = sampleCounterArena + 1;
        else
            fread(fid,arenaRecordLength-1,'uchar');
        end
        bytesRead = bytesRead + arenaRecordLength;
    else 
        disp(['Unknown record format ' ch ]);
        fclose(fid);
        retVal = 0;
        return;
    end
   
end

fclose(fid);

%% save all data
disp('saving data...');
[~,name,ext] = fileparts(filename);
fnOut = [name '.mat'];
dr = '/Users/hupeiyao/Movies/PH/DATA/';
if exist(dr,'dir') == 0; mkdir(dr); end

save ([dr fnOut], 'tetrodeData', 'tetrodeTimestamp','tetrodeChannel','tetrodeUnit',...
    'tetrodeFS', 'tetrodeSamplesPerBlock', 'roomTimeStamps', 'roomXY');

retVal = 1;

