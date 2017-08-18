FileExtn = '.txt';
ExcludeFiles = {'ProcSettings';
                'coords';
                'Coords';
                'Summary';
                'WithDistances'};

xCoordCol = 1;	// x coordinates column
yCoordCol = 2;  // y coordinates column
areasCol = 3;	// area represented by each xy pair (e.g. area of a cluster if xy is the centroid)
densityCol = 4; // intensity represented by each xy pair (e.g. staining brightness within the area given above).

xMin = 0;
xMax = 1500;
yMin = 0;
yMax = 1500;

FurthestFriendID = 1;   % Measure relative to this nth nearest neighbour

FileList = GetListOfData(FileExtn, ExcludeFiles);

% ask some questions

prompt = {'Nth Neighbour:','xCoord Column:','yCoord Column:','xMin','xMax:','yMin:','yMax:'};
dlg_title = 'Enter the NEW settings to apply...';
num_lines = 1;
defaults = {num2str(FurthestFriendID),num2str(xCoordCol),num2str(yCoordCol),num2str(xMin),num2str(xMax),num2str(yMin),num2str(yMax)};
answer = inputdlg(prompt,dlg_title,num_lines,defaults);

if ~isempty(answer)
    FurthestFriendID = str2double(answer(1,1));
    xCoordCol = str2double(answer(2,1));
    yCoordCol = str2double(answer(3,1));
    xMin = str2double(answer(4,1));
    xMax = str2num([answer{5,1}]);
    yMin = str2num([answer{6,1}]);
    yMax = str2num([answer{7,1}]);
else
    error('Cancelled?! So rude.');
end

% begin processing

ProcMessage = ['Processing ',num2str(size(FileList,1)),' files...'];
disp('=====================');
disp(ProcMessage);
disp('=====================');

% split file list into Ch1 and Ch2 files
FileList_Ch1 = FileList(find(~cellfun(@isempty,strfind(FileList,'Ch1'))),1);
FileList_Ch2 = FileList(find(~cellfun(@isempty,strfind(FileList,'Ch2'))),1);

data_stats = {...
                'FileID',...
                'Nth',...
                'Avg dist to Nth (Self)',...
                'StdDev(Avg,Self)',...
                'Avg sum of dists to Nth (Self)',...
                'StdDev(Sum, Self)',...
                'Avg dist to Nth (other ch)',...
                'StdDev(Avg) (Other Ch)',...
                'Avg sum of dist to Nth',...
                'StdDev(Sum) (Other Ch)',...
                'Normalised Distances',...
                'Distance Polarisation',...
                'Total Points',...
                'Total Points  (Other Ch)'...
                'Average Area'...
				'Average Intensity'...
                };
            
data_nb_self_stats = data_stats;
data_nb_other_stats = data_stats;

data_DN = struct;
data_DP = struct;

data_DN_nonboundary_self = struct;
data_DP_nonboundary_self = struct;

data_DN_nonboundary_other = struct;
data_DP_nonboundary_other = struct;


for t = 1:size(FileList,1)
    
    %Load the file
    CurrentFile  = FileList{t,1};
    FileNameOnly = strsplit(CurrentFile,FileExtn);
    SaveFileName = strcat(FileNameOnly{1},' - NN-',num2str(FurthestFriendID),FileExtn);

    data = importdata(CurrentFile);

    % Set the columns to save calcs to
    MainDistToNCol = size(data.data,2) + 1; % Distance to nth in this column
    MainSumToNCol = size(data.data,2) + 2;  % Sum of distances to nearest 10 in this column
    MainBoundaryFlagCol = size(data.data,2) + 3;  % Sum of distances to nearest 10 in this column
    
    OtherDistToNCol = size(data.data,2) + 4; % Distance to nth in this column
    OtherSumToNCol = size(data.data,2) + 5;  % Sum of distances to nearest 10 in this column
    OtherBoundaryFlagCol = size(data.data,2) + 6;  % Sum of distances to nearest 10 in this column
    
    NormCol = size(data.data,2) + 7; 
    PolarisationCol = size(data.data,2) + 8; 
    
    % Find and import it's matching channel
    FileList_Ch_line = find(~cellfun(@isempty,strfind(FileList_Ch1,CurrentFile)));
    
    if isempty(FileList_Ch_line)
        FileList_Ch_line = find(~cellfun(@isempty,strfind(FileList_Ch2,CurrentFile)));
        data_otherCh = importdata(FileList_Ch1{FileList_Ch_line,1});
    else
        data_otherCh = importdata(FileList_Ch2{FileList_Ch_line,1});
    end
    
    % measure the distance to nth nearest event, for all events
    for eventID1 = 1:size(data.data,1)
        x1 = data.data(eventID1,xCoordCol);
        y1 = data.data(eventID1,yCoordCol);

        % Relative to self
        Geometry_Self_Temp = zeros(size(data.data,1),1);

        for eventID2 = 1:size(data.data,1)
            x2 = data.data(eventID2,xCoordCol);
            y2 = data.data(eventID2,yCoordCol);

            deltaX = x1-x2;
            deltaY = y1-y2;

            EucDist = sqrt((deltaX*deltaX) + (deltaY*deltaY));
            Geometry_Self_Temp(eventID2,1) = EucDist;
        end

        % Sort the distances by increasing value and delete the
        % self-to-self distance (always zero)
        Geometry_Self_Temp2 = sort(Geometry_Self_Temp(:,1));
        Geometry_Self_Temp2(1,:) = [];

        FurthestFriend = Geometry_Self_Temp2(FurthestFriendID,1);
        AllFriendDistance = sum(Geometry_Self_Temp2(1:FurthestFriendID,1));
        
        % test if the distance to the nth neighbour is greater than the
        % distance to the nearest boundary.
        dxMinEdge = x1 - xMin;
        dxMaxEdge = xMax - x1;
        dyMinEdge = y1 - yMin;
        dyMaxEdge = yMax - y1;

        dxbound = min(dxMinEdge,dxMaxEdge);
        dybound = min(dyMinEdge,dyMaxEdge);
        
        dAnyEdge = min(dxbound,dybound);
               
        if dAnyEdge < FurthestFriend
            IsBoundaryPoint = true;
        else
            IsBoundaryPoint = false;
        end
        
        data.data(eventID1,MainDistToNCol) = FurthestFriend;
        data.data(eventID1,MainSumToNCol) = AllFriendDistance;
        data.data(eventID1,MainBoundaryFlagCol) = IsBoundaryPoint;
              
        % Relative to other channel
        Geometry_Other_Temp = zeros(size(data_otherCh.data,1),1);

        for eventID2 = 1:size(data_otherCh.data,1)
            x2 = data_otherCh.data(eventID2,xCoordCol);
            y2 = data_otherCh.data(eventID2,yCoordCol);

            deltaX = x1-x2;
            deltaY = y1-y2;

            EucDist = sqrt((deltaX*deltaX) + (deltaY*deltaY));
            Geometry_Other_Temp(eventID2,1) = EucDist;
        end

        % Sort the distances by increasing value
        Geometry_Self_Temp2 = sort(Geometry_Other_Temp(:,1));

        FurthestFriendOther = Geometry_Self_Temp2(FurthestFriendID,1);
        AllFriendDistanceOther = sum(Geometry_Self_Temp2(1:FurthestFriendID,1));
        
        % test if the distance to the nth neighbour is greater than the
        % distance to the nearest boundary.
        dxMinEdge = x1 - xMin;
        dxMaxEdge = xMax - x1;
        dyMinEdge = y1 - yMin;
        dyMaxEdge = yMax - y1;

        dxbound = min(dxMinEdge,dxMaxEdge);
        dybound = min(dyMinEdge,dyMaxEdge);
        
        dAnyEdge = min(dxbound,dybound);
              
        if dAnyEdge < FurthestFriendOther
            IsBoundaryPoint = true;
        else
            IsBoundaryPoint = false;
        end
        
        data.data(eventID1,OtherDistToNCol) = FurthestFriendOther;
        data.data(eventID1,OtherSumToNCol) = AllFriendDistanceOther;
        data.data(eventID1,OtherBoundaryFlagCol) = IsBoundaryPoint;  
        
%         % Normalise the distance to other Ch by distance to self
%         data.data(eventID1,SumToNCol+5) = FurthestFriendOther/FurthestFriend;
    end
    
    data.data(:,NormCol) = data.data(:,OtherDistToNCol) ./ data.data(:,MainDistToNCol);
    data.data(:,PolarisationCol) = (data.data(:,MainDistToNCol) - data.data(:,OtherDistToNCol)) ./ (data.data(:,MainDistToNCol) + data.data(:,OtherDistToNCol));

    % add DN and DP data to the holder structure
    FileNameStripped = regexprep(FileNameOnly{1,1},'[^\w'']','');
    % all points
    data_DN.(['DN_',FileNameStripped]) = data.data(:,NormCol);
    data_DP.(['DP_',FileNameStripped]) = data.data(:,PolarisationCol);
    % non-boundary points - self
    data_DN_nonboundary_self.(['DN_',FileNameStripped]) = data.data(~data.data(:,MainBoundaryFlagCol),NormCol);
    data_DP_nonboundary_self.(['DP_',FileNameStripped]) = data.data(~data.data(:,MainBoundaryFlagCol),PolarisationCol);
    % non-boundary points - other
    data_DN_nonboundary_other.(['DN_',FileNameStripped]) = data.data(~data.data(:,OtherBoundaryFlagCol),NormCol);
    data_DP_nonboundary_other.(['DP_',FileNameStripped]) = data.data(~data.data(:,OtherBoundaryFlagCol),PolarisationCol);

    
    
    %calc stats for all points in self
    data_stats{t+1,1} = FileNameOnly{1};
    data_stats{t+1,2} = FurthestFriendID;
    
    data_stats{t+1,3} = sum(data.data(:,MainDistToNCol))/size(data.data,1);
    data_stats{t+1,4} = std(data.data(:,MainDistToNCol));
    data_stats{t+1,5} = sum(data.data(:,MainSumToNCol))/size(data.data,1);
    data_stats{t+1,6} = std(data.data(:,MainSumToNCol));
    
    data_stats{t+1,7} = sum(data.data(:,OtherDistToNCol))/size(data.data,1);
    data_stats{t+1,8} = std(data.data(:,OtherDistToNCol));
    data_stats{t+1,9} = sum(data.data(:,OtherSumToNCol))/size(data.data,1);
    data_stats{t+1,10} = std(data.data(:,OtherSumToNCol));
    
    data_stats{t+1,11} = sum(data.data(:,NormCol)/size(data.data,1));
    data_stats{t+1,12} = sum(data.data(:,PolarisationCol))/size(data.data,1);

    data_stats{t+1,13} = size(data.data,1);
    data_stats{t+1,14} = size(data_otherCh.data,1);
    
    data_stats{t+1,15} = sum(data.data(:,areasCol))/size(data.data,1);
	data_stats{t+1,16} = sum(data.data(:,densityCol))/size(data.data,1);

    
    % calc stats for non-boundary (self) points
    data_nonboundary_self = data.data(~data.data(:,MainBoundaryFlagCol),:);
    
    data_nb_self_stats{t+1,1} = FileNameOnly{1};
    data_nb_self_stats{t+1,2} = FurthestFriendID;
    
    data_nb_self_stats{t+1,3} = sum(data_nonboundary_self(:,MainDistToNCol))/size(data_nonboundary_self,1);
    data_nb_self_stats{t+1,4} = std(data_nonboundary_self(:,MainDistToNCol));
    data_nb_self_stats{t+1,5} = sum(data_nonboundary_self(:,MainSumToNCol))/size(data_nonboundary_self,1);
    data_nb_self_stats{t+1,6} = std(data_nonboundary_self(:,MainSumToNCol));
    
    data_nb_self_stats{t+1,7} = sum(data_nonboundary_self(:,OtherDistToNCol))/size(data_nonboundary_self,1);
    data_nb_self_stats{t+1,8} = std(data_nonboundary_self(:,OtherDistToNCol));
    data_nb_self_stats{t+1,9} = sum(data_nonboundary_self(:,OtherSumToNCol))/size(data_nonboundary_self,1);
    data_nb_self_stats{t+1,10} = std(data_nonboundary_self(:,OtherSumToNCol));
    
    data_nb_self_stats{t+1,11} = sum(data_nonboundary_self(:,NormCol)/size(data_nonboundary_self,1));
    data_nb_self_stats{t+1,12} = sum(data_nonboundary_self(:,PolarisationCol))/size(data_nonboundary_self,1);

    data_nb_self_stats{t+1,13} = size(data_nonboundary_self,1);
    data_nb_self_stats{t+1,14} = size(data_otherCh.data,1);
    
    data_nb_self_stats{t+1,15} = sum(data_nonboundary_self(:,areasCol))/size(data_nonboundary_self,1);
	data_nb_self_stats{t+1,16} = sum(data_nonboundary_self(:,densityCol))/size(data_nonboundary_self,1);
    
    
    % calc stats for non-boundary (other channel) points
    data_nonboundary_other = data.data(~data.data(:,MainBoundaryFlagCol),:);
    
    data_nb_other_stats{t+1,1} = FileNameOnly{1};
    data_nb_other_stats{t+1,2} = FurthestFriendID;
    
    data_nb_other_stats{t+1,3} = sum(data_nonboundary_other(:,MainDistToNCol))/size(data_nonboundary_other,1);
    data_nb_other_stats{t+1,4} = std(data_nonboundary_other(:,MainDistToNCol));
    data_nb_other_stats{t+1,5} = sum(data_nonboundary_other(:,MainSumToNCol))/size(data_nonboundary_other,1);
    data_nb_other_stats{t+1,6} = std(data_nonboundary_other(:,MainSumToNCol));
    
    data_nb_other_stats{t+1,7} = sum(data_nonboundary_other(:,OtherDistToNCol))/size(data_nonboundary_other,1);
    data_nb_other_stats{t+1,8} = std(data_nonboundary_other(:,OtherDistToNCol));
    data_nb_other_stats{t+1,9} = sum(data_nonboundary_other(:,OtherSumToNCol))/size(data_nonboundary_other,1);
    data_nb_other_stats{t+1,10} = std(data_nonboundary_other(:,OtherSumToNCol));
    
    data_nb_other_stats{t+1,11} = sum(data_nonboundary_other(:,NormCol)/size(data_nonboundary_other,1));
    data_nb_other_stats{t+1,12} = sum(data_nonboundary_other(:,PolarisationCol))/size(data_nonboundary_other,1);

    data_nb_other_stats{t+1,13} = size(data_nonboundary_other,1);
    data_nb_other_stats{t+1,14} = size(data_otherCh.data,1);
    
    data_nb_other_stats{t+1,15} = sum(data_nonboundary_other(:,areasCol))/size(data_nonboundary_other,1);
	data_nb_other_stats{t+1,16} = sum(data_nonboundary_other(:,densityCol))/size(data_nonboundary_other,1);
   
    %write new centroids file
    dlmwrite(SaveFileName,data.data,'delimiter','\t');
    
    ProgMessage = ['Completed processing: ',CurrentFile];
    disp(ProgMessage);
end

% Write file containing DP data
disp('Saving file of pooled Distance Polarisation data');
DP_elements = fieldnames(data_DP);
maximum_size = max(structfun(@length, data_DP));
DP_tmp = NaN(maximum_size,length(DP_elements));
for i = 1:length(DP_elements)
    current_length = length(data_DP.(DP_elements{i}));
    DP_tmp(1:current_length,i) = data_DP.(DP_elements{i});
end
format_str_headers = ['%s',repmat('\t%s',1,size(DP_elements,1)-1),'\r\n'];
format_str_data = ['%f',repmat('\t%f',1,size(DP_elements,1)-1),'\r\n'];
SummaryFileName = strcat('Distance Polarisation Values (All).txt');
fid = fopen(SummaryFileName,'w');
fprintf(fid,format_str_headers,DP_elements{:,1});
for datarows = 1:size(DP_tmp,1)
    fprintf(fid,format_str_data,DP_tmp(datarows,:));
end
fid = fclose(fid);

% data_DN_nonboundary_self
% data_DN_nonboundary_other

% Write file containing DN data
disp('Saving file of pooled Distance-Normalised data');
DN_elements = fieldnames(data_DN);
maximum_size = max(structfun(@length, data_DN));
DN_tmp = NaN(maximum_size,length(DN_elements));
for i = 1:length(DN_elements)
    current_length = length(data_DN.(DN_elements{i}));
    DN_tmp(1:current_length,i) = data_DN.(DN_elements{i});
end
format_str_headers = ['%s',repmat('\t%s',1,size(DN_elements,1)-1),'\r\n'];
format_str_data = ['%f',repmat('\t%f',1,size(DN_elements,1)-1),'\r\n'];
SummaryFileName = strcat('Distance Polarisation Values (All).txt');
fid = fopen(SummaryFileName,'w');
fprintf(fid,format_str_headers,DN_elements{:,1});
for datarows = 1:size(DN_tmp,1)
    fprintf(fid,format_str_data,DN_tmp(datarows,:));
end
fid = fclose(fid);


%write all stats file
disp('Saving summary file for all clusters...');
% write headers
headers = data_stats(1,:);
format_str_headers = ['%s',repmat('\t%s',1,size(headers,2)-1),'\r\n'];
format_str_data = ['%s\t%d',repmat('\t%f',1,size(headers,2)-4),'\t%d\t%d\r\n'];
SummaryFileName = strcat('Summary - All points to NN-',num2str(FurthestFriendID),'.txt');
fid = fopen(SummaryFileName,'w');
fprintf(fid,format_str_headers,headers{1,:});
% write data
for datarows = 2:size(data_stats,1)
    fprintf(fid,format_str_data,data_stats{datarows,:});
end
fid = fclose(fid);


%write summary stats file
disp('Saving summary file for non-boundary (self) clusters...');
% write headers
headers2 = data_nb_self_stats(1,:);
format_str_headers = ['%s',repmat('\t%s',1,size(headers2,2)-1),'\r\n'];
format_str_data = ['%s\t%d',repmat('\t%f',1,size(headers2,2)-4),'\t%d\t%d\r\n'];

SummaryFileName = strcat('Summary - Non-boundary (Self) to NN-',num2str(FurthestFriendID),'.txt');
fid = fopen(SummaryFileName,'w');
fprintf(fid,format_str_headers,headers2{1,:});
% write data
for datarows = 2:size(data_stats,1)
    fprintf(fid,format_str_data,data_nb_self_stats{datarows,:});
end
fid = fclose(fid);

%write summary stats file
disp('Saving summary file for non-boundary (other) clusters...');
% write headers
headers2 = data_nb_other_stats(1,:);
format_str_headers = ['%s',repmat('\t%s',1,size(headers2,2)-1),'\r\n'];
format_str_data = ['%s\t%d',repmat('\t%f',1,size(headers2,2)-4),'\t%d\t%d\r\n'];

SummaryFileName = strcat('Summary - Non-boundary (Other Ch) to NN-',num2str(FurthestFriendID),'.txt');
fid = fopen(SummaryFileName,'w');
fprintf(fid,format_str_headers,headers2{1,:});
% write data
for datarows = 2:size(data_stats,1)
    fprintf(fid,format_str_data,data_nb_other_stats{datarows,:});
end
fid = fclose(fid);


disp('=====================');
disp('Processing Complete!');
disp('=====================');
