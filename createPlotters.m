function [platp, detp, trackp] = createPlotters(scenario)
%% Using Hampton University Skyler Radar Sensor's dataset to track drones
% Work done at TAN's LAB (Time-sensitive networking (T), 
%                           AI-driven cybersecurity (A), 
%                    NextG communication networking (N), 
%         Time-series Analysis via Network science (TAN))
% Create plotters
f = figure("Units","normalized",OuterPosition=[0.2 0.2 0.45 0.6]);
hax = axes(f);
tp = theaterPlot('Parent',hax,'XLim', [-1000, 4000], 'YLim', [-3000 1500], 'ZLim', [-100 100]);
set(tp.Parent,'YDir','reverse', 'ZDir','reverse');
platp = platformPlotter(tp,'DisplayName','Platforms','MarkerFaceColor','k');
detp = detectionPlotter(tp,'DisplayName','Detections','MarkerSize',6,'MarkerFaceColor',[0.85 0.325 0.098],'MarkerEdgeColor','k','History',10000);
trackp = trackPlotter(tp, 'DisplayName', 'tracks','ConnectHistory','on');

%
trajp = trajectoryPlotter(tp, 'DisplayName','Target Trajectories','LineWidth',1, 'LineStyle','--');
sampTimes = linspace(0,18,20);
for i=1:3
    trajpose{i} = lookupPose(scenario.Platforms{i}.Trajectory,sampTimes );

end

plotTrajectory(trajp, trajpose);
end