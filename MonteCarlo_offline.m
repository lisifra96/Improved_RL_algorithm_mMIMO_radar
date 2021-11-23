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

%% Wcube and beampattern computation for the optimal algorithm 
Wcube_opt=zeros(Nt,Nt,NumberVaryingScenarios);          % cube containing the optimal W for each varying scenario
BP_Opt=zeros(nuBinsNumber,Tmax);                        % Beampattern evaluated with the optimal W for each time instant

% this part computes the Wcube and BP in the optimal case (Supposing the radar knows the position of the targets)
for u=1:NumberVaryingScenarios                              % loop over varying scenarios 
    Wcube_opt(:,:,u)=getWfromTargetIndexes_offline(WcubeStruct,NuTargetIndex{u});
    X=Wcube_opt(:,:,u)'*conj(aT_Mat);
    BP_Opt(:,ScenarioStartInstant(u):ScenarioStopInstant(u))=dot(X,X).'*ones(1,ScenarioStopInstant(u)-ScenarioStartInstant(u)+1);
end

%% Orthogonal algorithm: variable definition
% Definition of some variables
WaldTest_Stat_Ort=zeros(nuBinsNumber,1);                % array containing the Wald test statistic

% Monte Carlo averaged results variable initialization
DetectionFrequency_Ort=zeros(nuBinsNumber,Tmax);        % Detection frequency averaged over MC trials

%% Optimal algorithm: variable definition
% Definition of some variables
WaldTest_Stat_Opt=zeros(nuBinsNumber,1);                % array containing the Wald test statistic

% Monte Carlo averaged results variable initialization
DetectionFrequency_Opt=zeros(nuBinsNumber,Tmax);        % Detection frequency averaged over MC trials

%% Adaptive algorithm: variable definition
% Definition of some variables
WaldTest_Stat_Adaptive=zeros(nuBinsNumber,1);           % array containing the Wald test statistic
DetectionArray_Adaptive=zeros(nuBinsNumber,1);          % double(WaldTest_Stat>Thresh)

% Monte Carlo averaged results variable initialization
BPaverage_Adaptive=zeros(nuBinsNumber,Tmax);            % Beam pattern averaged over MC trials
DetectionFrequency_Adaptive=zeros(nuBinsNumber,Tmax);   % Detection frequency averaged over MC trials

%% SARSA algorithm: variable definition
% Definition of some variables
WaldTest_Stat_SARSA=zeros(nuBinsNumber,1);              % array containing the Wald test statistic
DetectionArray_SARSA=zeros(nuBinsNumber,1);             % double(WaldTest_Stat>Thresh)
state_array=zeros(1,Tmax);                              % array containing each state in a single MC run
action_array=zeros(1,Tmax);                             % array containing each action in a single MC run

% Monte Carlo averaged results variable initialization
WaldTest_StatMat_SARSA=zeros(nuBinsNumber,Tmax);        % Matrix containing all the WaldTest statistics averaged over MC trials
BPaverage_SARSA=zeros(nuBinsNumber,Tmax);               % Beam pattern averaged over MC trials
DetectionFrequency_SARSA=zeros(nuBinsNumber,Tmax);      % Detection frequency averaged over MC trials
QcubeAvg_SARSA=zeros(MaxDetectableTargets+1,MaxDetectableTargets+1,Tmax);   % Cube containing all the Q matrices over time averaged over MC trials
RelErrQ_average_SARSA=zeros(1,Tmax);                    % norm of the difference of the Q matrices at following time instants averaged over MC trials
Qnorm_average_SARSA=zeros(1,Tmax);                      % norm of the Q matrices at following time instants averaged over MC trials
rewardAvg_SARSA=zeros(1,Tmax);                          % reward averaged over MC trials
StateFreq=zeros(MaxDetectableTargets+1,Tmax);           % Matrix containing the probability of visiting each state at each time averaged ovre MC trials
ActionFreq=zeros(MaxDetectableTargets+1,Tmax);          % Matrix containing the probability of selecting each action at each time averaged ovre MC trials
epsilonValuesArrayAvg=zeros(1,Tmax);                    % array containing the ε value at each time instant averaged over all MC trials (adaptive ε algorithm)
alphaValuesArrayAvg=zeros(1,Tmax);                      % array containing the α value at each time instant averaged over all MC trials (adaptive α algorithm)

%% Monte Carlo loop
for m=1:MC_iter
    
    %% Adaptive algorithm: variable initialization
    W_Adaptive=W_ort;                           % The initial value of W is set to W_ort

    %% SARSA algorithm: variable initialization
    Q_SARSA=Q_SARSAinit;                        % Q matrix initialization
    state_prev_SARSA=1;                         % state initialization
    action_prev_SARSA=1;                        % action initialization
    W_SARSA=W_ort;                              % The initial value of W is set to W_ort
    if VaryingEpsilonFlag
        SARSAparam.epsilon=EpsilonInitial;      % ε initialization (adaptive ε algorithm)
    end
    if VaryingAlphaFlag
        SARSAparam.alpha=AlphaInitial;          % α initialization (adaptive α algorithm)
    end
    ExplorationFlagArray=ones(3,1);             % This array contains the last three samples of the ExplorationFlag
    PrevReward=0;                               % reward value at the previous iteration                                                            
    AbsDiffReward=0;                            % absolute value of the difference between the current reward and the previous reward                                                             
    
    %% Loop over time and angular bins
    
    for u=1:NumberVaryingScenarios                              % loop over varying scenarios         
        for t=ScenarioStartInstant(u):ScenarioStopInstant(u)    % loop over time within the given scenario
            
            % Generation of the array containing the coefficients that
            % account for the two way path loss and Radar Cross Section of
            % the targets
            Amplitude_array=GenAmplArray(NuTargetIndex{u},TargetEffAmpl{u},nuBinsNumber);            
            
            for l=1:nuBinsNumber                                % loop over angular bins
                
                % Generation of the noise vector
                noise=AR_gen_t_dist(M,Ntrans,rho,sigma2_c,lambda,scale);
                                
                %% Orthogonal algorithm: Wald test statistic computation
                h=kron(W_ort.'*aT_Mat(:,l),aR_Mat(:,l));
                y=Amplitude_array(l)*h+noise;
                NormSquare_h=sum(abs(h).^2);
                hat_alpha=sum(y.*conj(h))/NormSquare_h;
                WaldTest_Stat_Ort(l)=Wald_test(M,h,y,hat_alpha,TruncationLag,NormSquare_h);

                %% Optimal algorithm: Wald test statistic computation
                h=kron(Wcube_opt(:,:,u).'*aT_Mat(:,l),aR_Mat(:,l));
                y=Amplitude_array(l)*h+noise;
                NormSquare_h=sum(abs(h).^2);
                hat_alpha=sum(y.*conj(h))/NormSquare_h;
                WaldTest_Stat_Opt(l)=Wald_test(M,h,y,hat_alpha,TruncationLag,NormSquare_h);

                %% Adaptive algorithm: Wald test statistic computation
                h=kron(W_Adaptive.'*aT_Mat(:,l),aR_Mat(:,l));
                y=Amplitude_array(l)*h+noise;
                NormSquare_h=sum(abs(h).^2);
                hat_alpha=sum(y.*conj(h))/NormSquare_h;
                WaldTest_Stat_Adaptive(l)=Wald_test(M,h,y,hat_alpha,TruncationLag,NormSquare_h);
                
                %% SARSA algorithm: Wald test statistic computation
                h=kron(W_SARSA.'*aT_Mat(:,l),aR_Mat(:,l));
                y=Amplitude_array(l)*h+noise;
                NormSquare_h=sum(abs(h).^2);
                hat_alpha=sum(y.*conj(h))/NormSquare_h;
                WaldTest_Stat_SARSA(l)=Wald_test(M,h,y,hat_alpha,TruncationLag,NormSquare_h);                
            end

            %% Orthogonal algorithm: MC loop variables update
            DetectionFrequency_Ort(:,t)=DetectionFrequency_Ort(:,t)+double(WaldTest_Stat_Ort>Thresh);    

            %% Optimal algorithm: MC loop variables update
            DetectionFrequency_Opt(:,t)=DetectionFrequency_Opt(:,t)+double(WaldTest_Stat_Opt>Thresh);    

            %% Adaptive algorithm: detection and W selection
            DetectionArray_Adaptive=WaldTest_Stat_Adaptive>Thresh;

            % The following instruction select the optimal W
            % matrix given the current action
            [~,TargetIndexes]=maxk(WaldTest_Stat_Adaptive,min(sum(DetectionArray_Adaptive),MaxDetectableTargets));
            [W_Adaptive] = getWfromTargetIndexes_offline(WcubeStruct,TargetIndexes);

            %% Adaptive algorithm: MC loop variables update
            DetectionFrequency_Adaptive(:,t)=DetectionFrequency_Adaptive(:,t)+DetectionArray_Adaptive;
    
            % Beampattern update
            X=W_Adaptive'*conj(aT_Mat);
            BPaverage_Adaptive(:,t)=BPaverage_Adaptive(:,t)+dot(X,X).';

            %% SARSA algorithm: detection and W selection
            DetectionArray_SARSA=WaldTest_Stat_SARSA>Thresh;

            % Generation of the struct GetStateRewardInputStruct
            GetStateRewardInputStruct.WaldTest_Stat=WaldTest_Stat_SARSA;
            GetStateRewardInputStruct.DetectionArray=DetectionArray_SARSA;

            % Generation of the struct SARSAinStruct
            SARSAinStruct.state_prev=state_prev_SARSA;
            SARSAinStruct.action_prev=action_prev_SARSA;
            SARSAinStruct.Q=Q_SARSA;

            % SARSA algorithm
            [SARSAoutStruct] = SARSA(GetStateReward,epsilonGreedy,GetStateRewardInputStruct,SARSAparam,SARSAinStruct);

            % The following instruction select the optimal W
            % matrix given the current action
            [~,TargetIndexes]=maxk(WaldTest_Stat_SARSA,SARSAoutStruct.action-1);   % array containing the indexes of the angular bins associated to the (action-1) highest decision statistics
            W_SARSA(:,:) = getWfromTargetIndexes_offline(WcubeStruct,TargetIndexes);

            %% SARSA algorithm: Adaptive ε and α algorithms
            % Update of the array containing the last 3 ExplorationFlag
            % values
            ExplorationFlagArray(1:2)=ExplorationFlagArray(2:3);
            ExplorationFlagArray(3)=SARSAoutStruct.ExplorationFlag;

            % Computation of the absolute value of the difference
            % between the current and the previous reward value
            AbsDiffReward=abs(SARSAoutStruct.reward-PrevReward);
            AdaptiveAlgorithmUpdateBit=~(ExplorationFlagArray(1)||ExplorationFlagArray(2));

            % Adaptive epsilon algorithm
            if VaryingEpsilonFlag

                epsilonValuesArrayAvg(t)=epsilonValuesArrayAvg(t)+SARSAparam.epsilon;

                if AdaptiveAlgorithmUpdateBit
                    if AbsDiffReward<VaryingEpsilonThresh1
                        SARSAparam.epsilon=max(SARSAparam.epsilon*EpsilonDecreaseFactor,EpsilonMin);
                    else
                        if AbsDiffReward>VaryingEpsilonThresh2
                            SARSAparam.epsilon=EpsilonMax;
                        else
                            SARSAparam.epsilon=min(SARSAparam.epsilon*EpsilonIncreaseFactor1,EpsilonMax);
                        end
                    end
                end
            end

            % Adaptive alpha algorithm
            if VaryingAlphaFlag

                alphaValuesArrayAvg(t)=alphaValuesArrayAvg(t)+SARSAparam.alpha;

                if AdaptiveAlgorithmUpdateBit
                    if AbsDiffReward<VaryingAlphaThresh1
                        SARSAparam.alpha=max(SARSAparam.alpha*AlphaDecreaseFactor,AlphaMin);
                    else
                        if AbsDiffReward>VaryingAlphaThresh2
                            SARSAparam.alpha=AlphaMax;
                        else
                            SARSAparam.alpha=min(SARSAparam.alpha*AlphaIncreaseFactor1,AlphaMax);
                        end
                    end
                end
            end
            PrevReward=SARSAoutStruct.reward;

            %% SARSA algorithm: update of state, action and Q
            % Q update
            Q_SARSA(:,:)=SARSAoutStruct.Q;

            % Update state and action
            state_prev_SARSA=SARSAoutStruct.state;
            action_prev_SARSA=SARSAoutStruct.action;

            %% SARSA algorithm: MC loop variables update

            rewardAvg_SARSA(t)=rewardAvg_SARSA(t)+SARSAoutStruct.reward;
            DetectionFrequency_SARSA(:,t)=DetectionFrequency_SARSA(:,t)+DetectionArray_SARSA;
            QcubeAvg_SARSA(:,:,t)=QcubeAvg_SARSA(:,:,t)+Q_SARSA;

            state_array(t)=SARSAoutStruct.state;
            action_array(t)=SARSAoutStruct.action;
            WaldTest_StatMat_SARSA(:,t)=WaldTest_Stat_SARSA;

            % Compute Beampattern
            X=W_SARSA'*conj(aT_Mat);
            BPaverage_SARSA(:,t)=BPaverage_SARSA(:,t)+dot(X,X).';
                        
        end
    end
    
    %% SARSA algorithm: computation of StateFreq and ActionFreq matrices
    for x=1:MaxDetectableTargets+1
       
        StateFreq(x,:)=StateFreq(x,:)+double(state_array==x);
        ActionFreq(x,:)=ActionFreq(x,:)+double(action_array==x);
        
    end
    
end
%% Orthogonal algorithm: MC loop variables computation
DetectionFrequency_Ort=DetectionFrequency_Ort/MC_iter;

%% Optimal algorithm: MC loop variables computation
DetectionFrequency_Opt=DetectionFrequency_Opt/MC_iter;

%% Adaptive algorithm: MC loop variables computation
DetectionFrequency_Adaptive=DetectionFrequency_Adaptive/MC_iter;
BPaverage_Adaptive=BPaverage_Adaptive/MC_iter;

%% SARSA algorithm: MC loop variables computation
DetectionFrequency_SARSA=DetectionFrequency_SARSA/MC_iter;
BPaverage_SARSA=BPaverage_SARSA/MC_iter;
QcubeAvg_SARSA=QcubeAvg_SARSA/MC_iter;
rewardAvg_SARSA=rewardAvg_SARSA/MC_iter;
StateFreq=StateFreq/MC_iter;
ActionFreq=ActionFreq/MC_iter;
epsilonValuesArrayAvg=epsilonValuesArrayAvg/MC_iter;
alphaValuesArrayAvg=alphaValuesArrayAvg/MC_iter;
RelErrQ_average_SARSA(1,1)=norm(QcubeAvg_SARSA(:,:,1)-Q_SARSAinit,'fro')/norm(QcubeAvg_SARSA(:,:,1),'fro');
RelErrQ_average_SARSA(1,2:Tmax)=vecnorm(vecnorm(QcubeAvg_SARSA(:,:,2:Tmax)-QcubeAvg_SARSA(:,:,1:(Tmax-1)),2,1),2,2)./vecnorm(vecnorm(QcubeAvg_SARSA(:,:,2:Tmax),2,1),2,2);
