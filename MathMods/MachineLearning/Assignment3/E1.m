%%	Exercise 1 - Linear mapping

%%	Initialization
clear all; clc;
load('Adot.mat');

theta = pi/3;
V = [cos(theta) -sin(theta);sin(theta) cos(theta)];

Y = V * X;
Z = V' * Y;
