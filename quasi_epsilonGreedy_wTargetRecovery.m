%  Copyright 2021 Francesco Lisi
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

function [action,ExplorationFlag] = quasi_epsilonGreedy_wTargetRecovery(epsilon,stateStruct,Q)
% This function chooses one action according to the quasi ε-greedy policy with
% target recovery. Refer to the article ... for the definition of this
% policy,
% INPUT VARIABLES
% epsilon:          ε real number within the interval (0,1)
% stateStruct:      structure containing previous state and current state
% Q:                current state-action value function
% OUTPUT VARIABLES
% action:           selected action

    PrevState=stateStruct.PrevState;
    CurrState=stateStruct.CurrState;
    
    % When the system loses at least one target it performs the optimal
    % action associated to the previous state to recover the lost targets more rapidly.
    if CurrState<PrevState
        [~,action]=max(Q(PrevState,:));
        ExplorationFlag=1;
    else
        
        % quasi e-greedy policy
        tmp=rand(1);
        [~,MaxIndex]=max(Q(CurrState,:));
        if tmp<epsilon
            Naction=size(Q,2);
            if CurrState~=Naction
                RandomActionSet=setdiff(CurrState:Naction,MaxIndex);
                action=RandomActionSet(randsample(length(RandomActionSet),1));
                ExplorationFlag=1;
            else
                action=MaxIndex;
                ExplorationFlag=0;
            end
        else
            action=MaxIndex;
            ExplorationFlag=0;
        end
        
    end

end
