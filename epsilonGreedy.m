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

function [action,ExplorationFlag] = epsilonGreedy(epsilon,stateStruct,Q)
% This function chooses one action according to the ε-greedy policy. The policy
% consist in selecting the best action with probability 1-ε, and choosing a
% random action between the remaining ones with probability ε.
% INPUT VARIABLES
% epsilon:          ε real number within the interval (0,1)
% stateStruct:      structure containing previous state and current state
% Q:                current state-action value function
% OUTPUT VARIABLES
% action:   selected action

    CurrState=stateStruct.CurrState;
    tmp=rand(1);                        
    [~,MaxIndex]=max(Q(CurrState,:));                                   % optimal action index
    if tmp<epsilon  % Exploration
        Naction=size(Q,2);
        RandomActionSet=setdiff(1:Naction,MaxIndex);                    % set containing the actions to randomly choose from
        action=RandomActionSet(randsample(length(RandomActionSet),1));  % selection of the random action
        ExplorationFlag=1;
    else            % Exploitation
        action=MaxIndex;
        ExplorationFlag=0;
    end

end

