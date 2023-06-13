clc;
clear all
close all;
warning off;
addpath(genpath('Utils'));

%% Puppets (Figure 1 and Figure 2)
Puppets_ReproduceFigures;

%% Torus (Figure 3 and Figure 4)
CommonPolodial_GetParameters;
Simulations_GetData;
Simulations_AnalyzePair;

%% Condition monitoring (Figure 5 and Table 1)
CM_ReproduceFigures;
CM_ReproduceTable;

%% E-Nose (Figure 6 and Table 2)
ENose_ReproduceFigures;
ENose_ReproduceTable;

%% Noisy MNIST (Figure 8 and Table 3)
NoisyMnist_ReproduceTable;
NoisyMnist_ReproduceFigures;

