function [numClusters, clusterTracks, clusterDetections, clusterProbabilities] = getClusterData(infolog, detlog, tracklog, step)
%% Using Hampton University Skyler Radar Sensor's dataset to track drones
% Work done at TAN's LAB (Time-sensitive networking (T), 
%                           AI-driven cybersecurity (A), 
%                    NextG communication networking (N), 
%         Time-series Analysis via Network science (TAN))

info = infolog{step};
detections = detlog{step};
numClusters = numel(info.Clusters);

% Retrieve tracks in cluster
initialTracks = tracklog{step -1};
clusterTracks = cell(1,numClusters);
clusterDetections = cell(1,numClusters);
clusterProbabilities = cell(1,numClusters);

for c=1:numClusters
    clusterTrackIDs = info.Clusters{c}.TrackIDs;
    clusterTracks{c} = initialTracks(ismember([initialTracks.TrackID],clusterTrackIDs));
    clusterDetections{c} = [detections{info.Clusters{c}.DetectionIndices}];
    clusterProbabilities{c} = info.Clusters{c}.MarginalProbabilities;
end
