function [detHistory, trackHistory] = getDetectionHistoryGNN(infoLog, detLog, trackLog, tID)
%% Using Hampton University Skyler Radar Sensor's dataset to track drones
% Work done at TAN's LAB (Time-sensitive networking (T), 
%                           AI-driven cybersecurity (A), 
%                    NextG communication networking (N), 
%         Time-series Analysis via Network science (TAN))

detHistory = objectDetection.empty;
for i=1:numel(infoLog)

    curInfo = infoLog{i};
    existTrack = any(curInfo.TrackIDsAtStepBeginning == tID);
    if ~existTrack
        % check if track was created
        if any(curInfo.InitiatedTrackIDs == tID)
            % Add initial detection
            detHistory(end+1) = detLog{i}{find(curInfo.InitiatedTrackIDs == tID)};
        end
        continue
    end
    
    if any(curInfo.DeletedTrackIDs == tID)
        %track was deleted
        break
    end

    trackMatches = find(curInfo.Assignments(:,1) == tID );
    assignedDetectionIndices = curInfo.Assignments(trackMatches, 2);
    assignedDets = [detLog{i}{assignedDetectionIndices}];
    for j=1:numel(assignedDets)
        detHistory(end+1) =  assignedDets(j);
    end
end

trackarray = [trackLog{:}];
alltrackids = [trackarray.TrackID];
trackHistory = trackarray(alltrackids == tID);