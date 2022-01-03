clear all;
close all;
clc;

for flowField=1:3
    markers={'*','o','s'};
    maxpercenterror=[];
    meanpercenterror=[];
    maxnormdiff=[];
    meannormdiff=[];
    X=[];
    for numIterations=0:10
        X=[X,numIterations];
        [metrics,~,~,~,~,~,~]=OccMotionTomography(numIterations,flowField);
        meanpercenterror=[meanpercenterror, metrics(2)];
        %maxpercenterror=[maxpercenterror, metrics(1)];
        %maxnormdiff=[maxnormdiff, metrics(3)];
        %meannormdiff=[meannormdiff, metrics(4)];
    end
    plot(X,meanpercenterror,markers{flowField});
    hold on;
end
legend('Flow field 1','Flow field 2','Flow field 3')
title('Number of iterations vs Mean Error');
xlabel('Number of Iterations') 
ylabel('Mean Error') 
hold off