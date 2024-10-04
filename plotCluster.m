function plotCluster(f, tracks, detections, probabilities, color)
%% Using Hampton University Skyler Radar Sensor's dataset to track drones
% Work done at TAN's LAB (Time-sensitive networking (T), 
%                           AI-driven cybersecurity (A), 
%                    NextG communication networking (N), 
%         Time-series Analysis via Network science (TAN))

hax = findobj(f,'Type','Axes');
tp = theaterPlot("Parent",hax(1));

detp = detectionPlotter(tp,'DisplayName','Cluster Detections','MarkerFaceColor', color);
trackp = trackPlotter(tp, 'DisplayName', 'Cluster Tracks','FontSize',12,'MarkerSize',12,'MarkerFaceColor',color);

[trackpos,trackcov] = getTrackPositions(tracks,'constvel');
plotTrack(trackp,trackpos,trackcov, arrayfun(@(x) string(x.TrackID), tracks'));

meas = cat(2,detections.Measurement)';
measCov = cat(3,detections.MeasurementNoise);
plotDetection(detp,meas,measCov);

yoffset = -15; % for offsetting text from line
for t=1:numel(tracks)
% line(hax(1),[meas(:,1)';repmat(trackpos(t,1),1,numel(detections))],[meas(:,2)';repmat(trackpos(t,2),1,numel(detections))],'Color',color,'HandleVisibility','off','LineWidth', 2);
% add probability
for d=1:numel(detections)
    line(hax(1),[meas(d,1)';trackpos(t,1)],[meas(d,2)'; trackpos(t,2)],'Color',[color,0.3+0.7* probabilities(d,t)],'HandleVisibility','off','LineWidth', 2);
    text(0.5*(meas(d,1)+trackpos(t,1)),yoffset+0.5*(meas(d,2)+trackpos(t,2)),[num2str(round(100*probabilities(d,t))) '%'],"FontSize",12);
end
end