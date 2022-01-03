clear all;
close all;
clc;

flowField = 1;
numIterations = 5;
[metrics,X,Y,VectorField,true_paths,approximate_paths,VectorFieldApprox]=OccMotionTomography(numIterations,flowField);

figure
quiver(X,Y,VectorField(1,:)',VectorField(2,:)', 'Color','black');
hold on;

for i=1:size(true_paths,3)
    plot(true_paths(1,:,i),true_paths(2,:,i));
    hold on;
end
xlim([-0.2,1.2]);
ylim([-0.2,1.2]);
title('True Field and Paths');
hold off

figure
quiver(X,Y,VectorFieldApprox(1,:)',VectorFieldApprox(2,:)', 'Color','black');
hold on

for i=1:size(approximate_paths,3)
    plot(approximate_paths(1,:,i),approximate_paths(2,:,i));
end
xlim([-0.2,1.2]);
ylim([-0.2,1.2]);
if(numIterations == 1)
    iterationstring = " Iteration";
else
    iterationstring = " Iterations";
end
    title("Approximate Field and Paths with " + numIterations + iterationstring);
hold off