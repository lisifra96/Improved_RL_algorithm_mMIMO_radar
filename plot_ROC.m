%%  Copyright 2021 Francesco Lisi
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

% This script plots the ROC curves for the different proposed algorithms.
% The PFA used in the x axis are the theoretical ones and not the effective
% ones computed via Monte Carlo simulations. Thus the results are valid
% only if the CFAR property is verified by the system.

clear all;
% close all;
clc;

SimulationId={'Results_2021_09_03_13_23_PFA_0001_Scenario_2_smoothed'};
Nsim=length(SimulationId);

% Retrieve Path of the folder where the results will be stored
PathName_ResultsNew=what('Results');
PathName_ResultsNew=PathName_ResultsNew.path;

load([PathName_ResultsNew '/' SimulationId{1}])
NumberOfTargets=length(NuTargetIndex{1});

TargetSNRdb_Mat=zeros(Nsim,NumberOfTargets);

% These flags are 1 if the corresponding algorithm was computed in the
% results that are being merged, otherwise it is 0
if exist('DetectionFrequency_Ort','var')
    OrtFlag=1;
else
    OrtFlag=0;
end

if exist('DetectionFrequency_Opt','var')
    OptFlag=1;
else
    OptFlag=0;
end

if exist('DetectionFrequency_Adaptive','var')
    AdptFlag=1;
else
    AdptFlag=0;
end

if exist('DetectionFrequency_SARSA','var')
    SARSAFlag=1;
else
    SARSAFlag=0;
end

if OrtFlag
    % Orthogonal merged results variable initialization
    DetectionFrequencyTarget_Ort=zeros(Nsim,Nloop,NumberOfTargets);
end

if OptFlag
    % Optimal merged results variable initialization
    DetectionFrequencyTarget_Opt=zeros(Nsim,Nloop,NumberOfTargets);
end

if AdptFlag
    % Adaptive merged results variable initialization
    DetectionFrequencyTarget_Adaptive=zeros(Nsim,Nloop,NumberOfTargets);
end

if SARSAFlag
    % SARSA merged results variable initialization
    DetectionFrequencyTarget_SARSA=zeros(Nsim,Nloop,NumberOfTargets);
end

for s=1:Nsim
   
    load([PathName_ResultsNew '/' SimulationId{s}])

    TargetSNRdb_Mat(s,:)=TargetSNRdB{1};
    
    for nn=1:Nloop
        
    StartIndex=strfind(FileName_Results,LoopVar)+length(LoopVar);
    StopIndex=strfind(FileName_Results,'_Scenario');
        if isnumeric(loop_array{nn})
            FileName_Results=[FileName_Results(1:StartIndex), sprintf('%s',num2str(loop_array{nn})) , FileName_Results(StopIndex:end)];
        else
            if ischar(loop_array{nn})
                FileName_Results=[FileName_Results(1:StartIndex), sprintf('%s',loop_array{nn}) , FileName_Results(StopIndex:end)];
            else
                FileName_Results=[FileName_Results(1:StartIndex), sprintf('%s',func2str(loop_array{nn})) , FileName_Results(StopIndex:end)];
            end
        end
        FileName_Results(find(FileName_Results=='.'))='';

        load([PathName_ResultsNew '/' FileName_Results])
        
%         Tmax=50;
        
        for m=1:NumberOfTargets
            if OrtFlag
                % Orthogonal merged results variable initialization
                DetectionFrequencyTarget_Ort(s,nn,m)=DetectionFrequency_Ort(NuTargetIndex{1}(m),Tmax);
            end

            if OptFlag
                % Optimal merged results variable initialization
                DetectionFrequencyTarget_Opt(s,nn,m)=DetectionFrequency_Opt(NuTargetIndex{1}(m),Tmax);
            end

            if AdptFlag
                % Adaptive merged results variable initialization
                DetectionFrequencyTarget_Adaptive(s,nn,m)=DetectionFrequency_Adaptive(NuTargetIndex{1}(m),Tmax);
            end

            if SARSAFlag
                % SARSA merged results variable initialization
                DetectionFrequencyTarget_SARSA(s,nn,m)=DetectionFrequency_SARSA(NuTargetIndex{1}(m),Tmax);
            end

        end
    end
    
end

%% Figures

% Figure parameters
line_width=2;
FontSize=18;
Marker='o';

% Starting figure counters
nFigSingle=101;
LineStyle_array={'-','--',':','-.'};

loop_array=cell2mat(loop_array);

for s=1:Nsim
    
    LineStyle=LineStyle_array{s};
    
    for m=1:NumberOfTargets
        
        figure(m)
        
        if OrtFlag
            PdLegend=sprintf('Orthogonal, SNR=%d dB',TargetSNRdb_Mat(s,m));
            semilogx(loop_array,DetectionFrequencyTarget_Ort(s,:,m),'Color','b','LineStyle',LineStyle,'Marker',Marker,'LineWidth',line_width,'DisplayName',PdLegend)
            legend('off')
            legend('show')
            legend({},'FontSize',FontSize)
            hold on;
        end
        
        if OptFlag
            PdLegend=sprintf('Optimal, SNR=%d dB',TargetSNRdb_Mat(s,m));
            semilogx(loop_array,DetectionFrequencyTarget_Opt(s,:,m),'Color','k','LineStyle',LineStyle,'Marker',Marker,'LineWidth',line_width,'DisplayName',PdLegend)
            legend('off')
            legend('show')
            legend({},'FontSize',FontSize)
            hold on;
        end
        
        if AdptFlag
            PdLegend=sprintf('Adaptive, SNR=%d dB',TargetSNRdb_Mat(s,m));
            semilogx(loop_array,DetectionFrequencyTarget_Adaptive(s,:,m),'Color','g','LineStyle',LineStyle,'Marker',Marker,'LineWidth',line_width,'DisplayName',PdLegend)
            legend('off')
            legend('show')
            legend({},'FontSize',FontSize)
            hold on;
        end
        
        if SARSAFlag
            PdLegend=sprintf('SARSA, SNR=%d dB',TargetSNRdb_Mat(s,m));
            semilogx(loop_array,DetectionFrequencyTarget_SARSA(s,:,m),'Color','r','LineStyle',LineStyle,'Marker',Marker,'LineWidth',line_width,'DisplayName',PdLegend)
            legend('off')
            legend('show')
            legend({},'FontSize',FontSize)
            hold on;
        end
        
        axis([loop_array(1) loop_array(end) 0 1]);
        xlabel('P_{FA}','FontSize',FontSize);
        ylabel('P_D','FontSize',FontSize);
        title(sprintf('ROC curve of the target corresponding to \\nu=%.2f',nu_array(NuTargetIndex{u}(m))),'FontSize',FontSize);
        set(gca,'FontSize',FontSize)
        grid on;

    end
end
