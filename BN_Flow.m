function [end_clock] = BN_Flow(ns,nc)
%=========================================================================%
% 2d coupled flow and biomass-nutrient dynamics simulator
%-------------------------------------------------------------------------%
% Inputs:
% ns: number for scenario
% nc: number for case
% user must first complete the input data script in the +Scenario folder.
% Example:
% To import input data from +Scenario/S1C1.m, 
% run BN_Flow(1,1) with ns = 1 and nc = 1.
% Simulation results are saved in Simulation/Scenario1/Case1 folder 
%            under different version number starting V1. 
%            If we run same example for multiple times, we will have 
%            multiple folders with names V1,V2,....           
%-------------------------------------------------------------------------%
% See, e.g., an inputcard +Scenario/S1C1.m for details about input data
%-------------------------------------------------------------------------%
% Outputs:
% end_clock: the total computational time
% This simulator also saves the plots of numerical solutions 
%      in the appropriate folder.
%-------------------------------------------------------------------------%
% supp.input_data reads the chosen example input data and translate into 
%                       appropriate units and format
% start_clock & end_clock measures the total computational time
% T is the final time constraint
% flag.restart is an input flag in +Scenario/S1C1.m which allows 
%                 to restart the simulation from the last recorded data
%                 when flag.restart = 1
% [B0,N0,L0] are the initial biomass, nutrient, and Lagrange multiplier
% G.b0 is the indicator for biofilm phase where B_*<= B <= B^*
% flag.Gb is empty until we first see the biofilm formation 
%            which can skip calculating the flow 
%            while there is no change in permeability
% G.kb is the permeability of biofilm defined in +Scenario/S1C1.m
% supp.BN_Brinkman solves the coupled flow and biomass-nutrient dynamics
%                  when G.kb > 0: biofilm is (partially) permeable. 
%                  We use the heterogeneous Brinkman flow.
% supp.BN_Stokes solves the coupled flow and biomass-nutrient dynamics 
%                when G.kb = 0: biofilm is impermeable. 
%                We use the Stokes flow.
% supp.record records data: [B,N,L] and [U,V,P,K] if flow is enabled
% supp.figGen saves figures: B, N, and velocity, P if flow is enabled 
%=========================================================================%
    %% import data
    supp.input_data;

    start_clock = tic;    
    %% initialize
    T = num.T;
    if flag.restart == 1
    tn = 0;
    while isfolder(sprintf('%s/plt%d',dir.version,tn)); tn = tn + 1; end
    tn = tn - 1;
    dir.tn = sprintf('%s/plt%d',dir.version,tn);
    B0 = load(sprintf('%s/B.txt',dir.tn));
    N0 = load(sprintf('%s/N.txt',dir.tn));
    L0 = load(sprintf('%s/L.txt',dir.tn));
    t = load(sprintf('%s/T.txt',dir.tn));
    else
    B0 = b.IC;
    N0 = n.IC;
    L0 = 0*B0;
    tn = 0;
    t = 0;
    end

    G.b0 = sparse(y.n,x.n);
    G.b0(B0 >= b.star2 & B0 <= b.star) = 1;
    G.bn = G.b0;
    flag.Gb = [];

    fprintf('tn = %d, t = %g [%s]\n', tn, t, unit.str_time);

    %% iteration
    if G.kb > 0
    supp.BN_Brinkman;
    else
    supp.BN_Stokes;
    end

    %% record at final time T

    end_clock = toc(start_clock);
    close('all','hidden');

    diary(sprintf('%s/record.txt',dir.version));
    diary on
    fprintf('Total clock time = %.4f minutes\n',end_clock./60);
    fprintf('========================================\n\n');
    diary off

    supp.record;
    supp.figGen;
end
