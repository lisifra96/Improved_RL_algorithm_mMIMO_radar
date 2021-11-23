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

% Script to compute the loop over Loop variable 

Nloop=length(loop_array);
 
if ~isequal(LoopOver,'LoopOverNt')
    Nt=1e2;
    if ~isequal(Processing,'online')
        % Loading of the presaved W cube (not required in the online processing)
        FileName_Wcubes=['Wcube_nuBins_',num2str(nuBinsNumber),'_Nt_',num2str(Nt),...
        '_MaxTarget_',num2str(MaxDetectableTargetsOffline),'_Ptot_',num2str(Ptot)];
        load([PathName_Wcubes '/' FileName_Wcubes]);
    end
end

for n=1:Nloop

    if isequal(LoopOver,'LoopOverNt')
        Nt=Nt_array{n};
        
        if ~isequal(Processing,'online')
            % Loading of the presaved W cube (not required in the online processing)
            FileName_Wcubes=['Wcube_nuBins_',num2str(nuBinsNumber),'_Nt_',num2str(Nt),...
            '_MaxTarget_',num2str(MaxDetectableTargetsOffline),'_Ptot_',num2str(Ptot)];
            load([PathName_Wcubes '/' FileName_Wcubes]);
        end
    end
    Nr=Nt;
    M=Nt*Nr;
    TruncationLag=ceil(M^(1/4));                            % Truncation lag for the estimate of the noise covariance matrix

    aT_Mat=exp(1i*2*pi*(0:Nt-1)'*nu_array');                % each column correspond to the steering vector for the associated ν
    aR_Mat=exp(1i*2*pi*(0:Nr-1)'*nu_array');

    W_ort=sqrt(Ptot/Nt)*eye(Nt);                            % uniform and orthogonal weight matrix 
    X=W_ort'*conj(aT_Mat);                                  % Auxiliary variable used to compute the beampattern
    BP_Ort=dot(X,X).'*ones(1,Tmax);                         % Beampattern
    
    if isequal(LoopOver,'LoopOverPFA')
        P_FA=P_FA_array{n};
    else
        P_FA=1e-4;
    end
    Thresh = chi2inv(1-P_FA,2);                 % Threshold (valid in the asymptotic massive MIMO regime)
    GetStateRewardInputStruct.Thresh=Thresh;    % GetStateRewardInputStruct is a structure needed in RL algorithms to compute the next state and reward

    if isequal(LoopOver,'LoopOverEpsilon')
        if isequal(epsilon{n},'VaryingEpsilon')
            VaryingEpsilonFlag=1;                   % This flag is set to 1 if the varying epsilon algorithm must be used
        else
            SARSAparam.epsilon=epsilon{n};
            VaryingEpsilonFlag=0;
        end
    else
         VaryingEpsilonFlag=1;
    end
    
    if isequal(LoopOver,'LoopOverEpsilonGreedy')
        epsilonGreedy=epsilonGreedy_array{n};
    else
        epsilonGreedy=@quasi_epsilonGreedy_wTargetRecovery;                   % epsilon greedy default function
    end
    
    if isequal(LoopOver,'LoopOverAlpha')
        if isequal(alpha{n},'VaryingAlpha')
            VaryingAlphaFlag=1;                     % This flag is set to 1 if the varying alpha algorithm must be used
        else
            SARSAparam.alpha=alpha{n};
            VaryingAlphaFlag=0;
        end
    else
         VaryingAlphaFlag=1;
    end
    
    
    if isequal(LoopOver,'LoopOverQInitialization')
        switch Qinitialization{n}
            case 'zero'
                Q_SARSAinit=zeros(MaxDetectableTargets+1,MaxDetectableTargets+1);       
            case 'identity'
                Q_SARSAinit=eye(MaxDetectableTargets+1,MaxDetectableTargets+1); 
            otherwise
                Q_SARSAinit=zeros(MaxDetectableTargets+1,MaxDetectableTargets+1);       
        end
    else
        Q_SARSAinit=eye(MaxDetectableTargets+1,MaxDetectableTargets+1); 
    end
    
    % Online mode:  online computation of the optimization algorithm to
    %               obtain the W matrix from the detected target position array
    % Offline mode: The W matrix are precomputed for each possible target
    %               position combination.
    switch Processing
        case 'online'
            MonteCarlo_online
        case 'offline'
            MonteCarlo_offline;
        otherwise
            error('Selected Processing parameter doesn''t exist');
    end

    %% Saving of the results
    % This section saves the results in a single matlab object
    
    if isnumeric(loop_array{n})
        FileName_Results=sprintf('Results_%s_%s_%s_%s_%s_%s_%s_Scenario_%s',year,month,day,hour,minute,LoopVar,num2str(loop_array{n}),num2str(TargetScenario));
    else
        if ischar(loop_array{n})
            FileName_Results=sprintf('Results_%s_%s_%s_%s_%s_%s_%s_Scenario_%s',year,month,day,hour,minute,LoopVar,loop_array{n},num2str(TargetScenario));                    
        else
            FileName_Results=sprintf('Results_%s_%s_%s_%s_%s_%s_%s_Scenario_%s',year,month,day,hour,minute,LoopVar,func2str(loop_array{n}),num2str(TargetScenario));        
        end
    end
    FileName_Results(find(FileName_Results=='.'))='';
    save([PathName_Results '/' FileName_Results])
    
end
