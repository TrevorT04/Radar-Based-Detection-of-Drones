function [detHistory, trackHistory] = getDetectionHistoryJPDA(infoLog, detLog, trackLog, tID, probThreshold)
%% Using Hampton University Skyler Radar Sensor's dataset to track drones
% Work done at TAN's LAB (Time-sensitive networking (T), 
%                           AI-driven cybersecurity (A), 
%                    NextG communication networking (N), 
%         Time-series Analysis via Network science (TAN))

detHistory = objectDetection.empty;
for i=1:numel(infoLog)

    curInfo = infoLog{i};
    existTrack = any(curInfo.TrackIDsAtStepBeginning == tID);
    isUnassigned = any(curInfo.UnassignedTracks == tID);
    if ~existTrack || isUnassigned
        % check if track was created
        if any(curInfo.InitializedTrackIDs == tID)
            % Add initial detection
            detHistory(end+1) = detLog{i}{find(curInfo.InitializedTrackIDs == tID)};
        end
        continue
    end
    
    if any(curInfo.DeletedTrackIDs == tID)
        %track was deleted
        break
    end

    % Find the cluster with tID
    hasTID = cellfun(@(x) any(x.TrackIDs == tID), curInfo.Clusters);
    cluster = curInfo.Clusters{hasTID};
    trackIndexInCluster = find(cluster.TrackIDs == tID );
    detIndexInCluster = find(cluster.MarginalProbabilities(1:end-1,trackIndexInCluster) > probThreshold);
    assignedDetectionIndices = cluster.DetectionIndices(detIndexInCluster);
    assignedDets = [detLog{i}{assignedDetectionIndices}];
    for j=1:numel(assignedDets)
        detHistory(end+1) =  assignedDets(j);
    end
end
trackLog = [trackLog{:}];
alltrackids = [trackLog.TrackID];
trackHistory = trackLog(alltrackids == tID);