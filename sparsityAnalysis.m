%% Visualize different sparsity constraint
clc
clear all;

% Hyper-prameters
size_ = 500;

% Generate L1 constraint
xL1 = linspace(0, 1, size_);
yL1 = 1 - xL1;

% Generate L1/2 constraint
xL1_2 = linspace(0, 1, size_);
yL1_2 = 1 - xL1_2;

% Generate L2 constraint
xL2 = linspace(0, 1, size_);
yL2 = 1 - xL2;

% Visualize these sparsity constraint
figure;
hold on
axis square
plot(xL1, yL1, 'r')
plot(xL1_2.^2, yL1_2.^2, 'b')
plot(xL2.^0.5, yL2.^0.5, 'k')
legend('L1','L2','L1/2')