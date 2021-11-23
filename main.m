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

%% Parameters

% compute current date and time. These strings are added to the saved
% results file name
c=clock;
year=num2str(c(1));
month=num2str(c(2));
if length(month)==1
    month=['0',month];
end
day=num2str(c(3));
if length(day)==1
    day=['0',day];
end
hour=num2str(c(4));
if length(hour)==1
    hour=['0',hour];
end
minute=num2str(c(5));
if length(minute)==1
    minute=['0',minute];
end
second=num2str(c(6));
if c(6)<10
    second=['0',second];
end

% Adaptive epsilon parameters
EpsilonMin=0.1;                         % Minimum epsilon value
EpsilonMax=0.8;                         % Maximum epsilon value
EpsilonIncreaseFactor1=2;               % Increase multiplicative coefficient 
EpsilonDecreaseFactor=0.8;              % Decrease multiplicative coefficient
VaryingEpsilonThresh1=0.5;              % 1st varying epsilon threshold
VaryingEpsilonThresh2=1.8;              % 2nd varying epsilon threshold
EpsilonInitial=EpsilonMax;              % Initial epsilon value

% Adaptive alpha parameters
AlphaMin=0.2;                           % Minimum alpha value
AlphaMax=0.6;                           % Maximum alpha value
AlphaIncreaseFactor1=2.5;               % Increase multiplicative coefficient 
AlphaDecreaseFactor=0.9;                % Decrease multiplicative coefficient
VaryingAlphaThresh1=0.5;                % 1st varying alpha threshold
VaryingAlphaThresh2=1.8;                % 2nd varying alpha threshold
AlphaInitial=AlphaMax;                  % Initial epsilon value

% In the following all the quantities that are formatted as cell are
% variables that you can loop on to compare the performances over different
% values

% SARSA parameters
epsilon={EpsilonMin,EpsilonMax,'VaryingEpsilon'};               % epsilon value
alpha={AlphaMin,AlphaMax,'VaryingAlpha'};                       % alpha value
SARSAparam.gamma = 0.8;                                         % gamma value of the total reward (controls the weight of early reward) [0,1]
Qinitialization={'zero','identity'};                            % Initial Q matrix

% function selection
GetStateReward=@GetStateReward_Pd;                                            % Selection of the function to use to compute the state and the reward from the collected data
nGetStateReward=length(GetStateReward);
epsilonGreedy_array={@epsilonGreedy,@quasi_epsilonGreedy,@quasi_epsilonGreedy_wTargetRecovery};          % Selection of the epsilon greedy policy type

% Retrieve Path of the folder where the presaved Wcube are (offline mode)
PathName_Wcubes=what('Wcubes');
PathName_Wcubes=PathName_Wcubes.path;

% Retrieve Path of the folder where the results will be stored
PathName_Results=what('Results');
PathName_Results=PathName_Results.path;

% Simulation parameters         
P_FA_array=num2cell(10.^(-5:1:-1));                             % Array containing different PFA values
Nt_array= num2cell(floor(logspace(1,2,1)));                     % Array containing different number of transmit antennas
nuBinsNumber=20;                                                % Number of ν bins 
nu_array=-0.5+(0:nuBinsNumber-1)'/nuBinsNumber;                 % Array defining the grid of ν=d/λ*sin(θ) values with d=λ/2
Ptot=1;                                                         % Total transmit power
MaxDetectableTargetsOffline=5;                                  
Processing='online';
if isequal(Processing,'offline')
    MaxDetectableTargets=MaxDetectableTargetsOffline;           % Maximum number of target that can be detected in offline mode (To maintain the Q matrix small)
else
    MaxDetectableTargets=5;                                     % Maximum number of target that can be detected in online mode (To maintain the Q matrix small)
end
GetStateRewardInputStruct.MaxTargetNumber=MaxDetectableTargets; % GetStateRewardInputStruct is a structure needed in RL algorithms to compute the next state and reward

%% Disturbance scenario definition
% The disturbance process is an AR(6) process with t-distributed innovation
% process

% coefficients of an AR(6) process
p=zeros(6,1);
p(1) = 0.5*exp(1j*2*pi*(-0.4));
p(2) = 0.6*exp(1j*2*pi*(-0.2));
p(3) = 0.7*exp(1j*2*pi*(0));
p(4) = 0.4*exp(1j*2*pi*(0.1));
p(5) = 0.5*exp(1j*2*pi*(0.3));
p(6) = 0.6*exp(1j*2*pi*(0.35));
rho = poly(p);
sigma2_c = 1;                                           % Variance of the process
lambda = 2;                                             % lambda value of the t distributed innovation process
scale=1/(lambda-1);
Ntrans=100*ceil(3*log(10)/abs( log( max(abs(p)) ) ) );  % number of samples that must be discarded (transitory phase)

%% Target Scenario definition
% In this section the user can choose one of the possible scenarios by
% setting the parameter TargetScenario.
%
% Tmax:                     Total number of discrete time instants
% NuTargetIndex:            Index of ν array defining the Target position
% TargetSNRdB:              Target SNR
% ScenarioChangeInstant:    time instants when the scenario changes

TargetScenario=2;             

switch TargetScenario
    case 1
        
        % Dynamic scenario with fixed SNR
        Tmax = 400;                                                                        
        NuTargetIndex={[5,13],[13],[13,17],[17]};           
        TargetSNRdB = {[-18,-21],[-21],[-21,-20],[-20]};     
        ScenarioChangeInstant=[101,201,301];                
    
    case 2
        
        % Static scenario with fixed SNR
        Tmax = 1*1e2; 
        NuTargetIndex={[7,16]};                                                
        TargetSNRdB = {[-20,-20]}; 
        ScenarioChangeInstant=[];        
                
    case 3
        
        % Scenario with two targets that have fixed angular position, but
        % variable SNR
        Tmax = 2*1e2;  
        NuTargetIndex=num2cell([7 16]'*ones(1,Tmax),1);                                                
        TargetSNRdB = num2cell(ones(2,1)*[linspace(-30,-20,floor(Tmax/2)), linspace(-20,-30,ceil(Tmax/2))],1); 
        ScenarioChangeInstant=2:Tmax;

    otherwise
        error('The selected scenario doesn''t exist');
end

% Number of subscenarios in each scenario
NumberVaryingScenarios=size(NuTargetIndex,2);

% Effective amplitude of the target, computed from the SNR and noise power
TargetEffAmpl=cellfun(@(x) sqrt(sigma2_c*10.^(x'/10)),TargetSNRdB,'UniformOutput',false);

% Start and stop instant of each subscenario
ScenarioStartInstant=[1, ScenarioChangeInstant];
ScenarioStopInstant=[ScenarioChangeInstant-1, Tmax];

%% Monte Carlo simulation 

% Parameters
MC_iter = 1*1e0;                                    % number of Monte Carlo runs
LoopOver='LoopOverNt';                   % Select which variable to loop on
LoopVar=LoopOver(9:end);                            % string added to the saved results' name

% OuterLoop is simply a script to simplify the readability of the code by
% not repeating the same lines for each case 

switch LoopOver
    case 'LoopOverNt'

        loop_array=Nt_array;
        OuterLoop
        
    case 'LoopOverPFA'

        loop_array=P_FA_array;
        OuterLoop
        
    case 'LoopOverEpsilon'

        loop_array=epsilon;
        OuterLoop
                
    case 'LoopOverAlpha'
        
        loop_array=alpha;
        OuterLoop

    case 'LoopOverEpsilonGreedy'

        loop_array=epsilonGreedy_array;
        OuterLoop

    case 'LoopOverQInitialization'

        loop_array=Qinitialization;
        OuterLoop

    otherwise
        error('Selected LoopOver variable doesn''t exist');
end

%% Figures

% Figure parameters
line_width=2;
FontSize=18;

% Select which algorithm's results to plot
OrthogonalPlot=1;
OptimalPlot=1;
AdaptivePlot=1;
SARSAPlot=1;

% Starting figure counters
nFigSingle=101;
Marker_array={'-','--',':','-.'};

% Discrete time array
time_array=1:Tmax;

if OrthogonalPlot
    
    nFigShared=1;
    Algorithm='Orthogonal';
    Colour='b';
    Marker='-';

    DetectionFrequency=DetectionFrequency_Ort;
    BPaverage=BP_Ort;
    PdLegend=sprintf('%s',Algorithm);

    Plot_Pd
end

if OptimalPlot

    nFigShared=1;
    Algorithm='Optimal';
    Colour='k';
    Marker='-';

    DetectionFrequency=DetectionFrequency_Opt;
    BPaverage=BP_Opt;
    PdLegend=sprintf('%s',Algorithm);

    Plot_Pd
    
end

if AdaptivePlot

    nFigShared=1;
    Algorithm='Adaptive';
    Colour='g';
    Marker='-';
    
    DetectionFrequency=DetectionFrequency_Adaptive;    
    BPaverage=BPaverage_Adaptive;
    PdLegend=sprintf('%s',Algorithm);

    Plot_Pd

end

if isnumeric(loop_array)
   loop_array=num2cell(loop_array); 
end

for n=1:Nloop
    
    StartIndex=strfind(FileName_Results,LoopVar)+length(LoopVar);
    StopIndex=strfind(FileName_Results,'_Scenario');
    if isnumeric(loop_array{n})
        FileName_Results=[FileName_Results(1:StartIndex), sprintf('%s',num2str(loop_array{n})) , FileName_Results(StopIndex:end)];
    else
        if ischar(loop_array{n})
            FileName_Results=[FileName_Results(1:StartIndex), sprintf('%s',loop_array{n}) , FileName_Results(StopIndex:end)];
        else
            FileName_Results=[FileName_Results(1:StartIndex), sprintf('%s',func2str(loop_array{n})) , FileName_Results(StopIndex:end)];
        end
    end
    FileName_Results(find(FileName_Results=='.'))='';

    
    % Retrieve Path of the folder where the results will be stored
    PathName_Results=what('Results');
    PathName_Results=PathName_Results.path;

    load([PathName_Results '/' FileName_Results])
        
    if SARSAPlot

        Algorithm='SARSA';
        Colour='r';
    
        Marker=Marker_array{mod(n-1,4)+1};
        nFigShared=1;
        DetectionFrequency=DetectionFrequency_SARSA;
        BPaverage=BPaverage_SARSA;
        QcubeAvg=QcubeAvg_SARSA;
        if isnumeric(loop_array{n})
            PdLegend=sprintf('%s, %s = %.2f',Algorithm,LoopVar,loop_array{n});
            LoopVarValue=num2str(loop_array{n});
        else
            if ischar(loop_array{n})
                LoopVarValue=loop_array{n};
                LoopVarValue=strrep(LoopVarValue,'_','\_');
                PdLegend=sprintf('%s, %s',Algorithm,LoopVarValue);
            else
                LoopVarValue=func2str(loop_array{n});
                LoopVarValue=strrep(LoopVarValue,'_','\_');
                PdLegend=sprintf('%s, %s',Algorithm,LoopVarValue);
            end

        end

        Plot_Pd

        rewardAvg=rewardAvg_SARSA;
        RelErrQ_average=RelErrQ_average_SARSA;

        PlotSetRL

    end
    
end

%% Plot of the Power Spatial Density of the noise process with vertical
% lines in the position where the target are placed
deltaf=1e-3;
f = -0.5:deltaf:0.5;
switch TargetScenario
    case 1
    
        for u=1:NumberVaryingScenarios
            if ~isempty(NuTargetIndex{u})
                plot_normalized_PSD2(rho);
                vline(nu_array(NuTargetIndex{u}),'--r');
            end
        end
    
    case 2
        
        plot_normalized_PSD2(rho);
        vline(nu_array(NuTargetIndex{u}),'--r');
                
    case 3
        
        plot_normalized_PSD2(rho);
        vline(nu_array(NuTargetIndex{1}),'--r'  );

        figure()
        plot(time_array, cell2mat(TargetSNRdB), 'r','LineWidth',line_width);
        grid on;
        xlabel('Discrete time k','FontSize',FontSize);
        ylabel('$SNR_{dB,k}$','fontsize',18,'interpreter','latex');
        title('$SNR_{dB,k} temporal evolution$','FontSize',FontSize,'interpreter','latex');
        set(gca,'FontSize',FontSize)
    
end
