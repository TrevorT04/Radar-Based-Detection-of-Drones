function plotTrackAndDets(fnew,detectionHistory, trackHistory)
%% Using Hampton University Skyler Radar Sensor's dataset to track drones
% Work done at TAN's LAB (Time-sensitive networking (T), 
%                           AI-driven cybersecurity (A), 
%                    NextG communication networking (N), 
%         Time-series Analysis via Network science (TAN))
% Create a new figure, copy hax, and plot all of detectionHistory and
% trackHistory in the same color


if isempty(trackHistory)
return
end

colors = lines(10);
tid = trackHistory(1).TrackID;
thisColor = colors(mod(tid,7)+1,:);



% Add a new theaterPlot to add the detections and track
parent = findobj(fnew,'type','axes');
tp = theaterPlot('Parent',parent(1),'XLim', [-1000, 4000], 'YLim', [-3000 1500], 'ZLim', [-100 100]);
set(tp.Parent,'YDir','reverse', 'ZDir','reverse');
detp = detectionPlotter(tp,'DisplayName',['Detections for Track' num2str(tid)],'MarkerFaceColor',thisColor);
trackp = trajectoryPlotter(tp, 'DisplayName', ['Track', num2str(tid)],'Color',thisColor, 'LineStyle','-','LineWidth',3);

% Plot
trackpositions = getTrackPositions(trackHistory,'constvel');
plotTrajectory(trackp,{trackpositions});

meas = cat(2,detectionHistory.Measurement)';
measCov = cat(3,detectionHistory.MeasurementNoise);
plotDetection(detp,meas,measCov);
% warning off
% s =struct(detp);
% m = s.PositionsLine.MarkerHandle;
% m.FaceColorData = 0*m.FaceColorData;

end
