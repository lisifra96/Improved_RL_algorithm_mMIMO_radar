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

function [SARSAoutStruct] = SARSA(GetStateReward,epsilonGreedy,GetStateRewardInputStruct,SARSAparam,SARSAinStruct)
% This function computes the SARSA algorithm using epsilonGreedy policy.
% INPUT VARIABLES
% GetStateReward:                   function handler to the function that computes
%                                   the new state and reward from GetStateRewardInputStruct
% epsilonGreedy:                    function handler to the ε-greedy policy function
% GetStateRewardInputStruct:        Input structure of GetStateReward (refer to specific
%                                   GetStateReward to see what input parameters it needs)
% SARSAparam.epsilon:               ε value
% SARSAparam.gamma:                 γ value
% SARSAparam.alpha:                 α value
% SARSAinStruct.state_prev:         MDP state at previous iteration
% SARSAinStruct.action_prev:        MDP action at previous iteration
% SARSAinStruct.Q:                  Q matrix
% OUTPUT VARIABLES
% SARSAoutStruct.state:             new state
% SARSAoutStruct.action:            new action
% SARSAoutStruct.reward:            reward 
% SARSAoutStruct.Q:                 Updated Q
% SARSAoutStruct.ExplorationFlag:   This flag is equal to 1 if the system
%                                   decides to explore the environment

epsilon=SARSAparam.epsilon;
gamma=SARSAparam.gamma;
alpha=SARSAparam.alpha;

state_prev=SARSAinStruct.state_prev;
action_prev=SARSAinStruct.action_prev;
Q=SARSAinStruct.Q;

% Computation of the state and reward
[state,reward] = GetStateReward(GetStateRewardInputStruct);

% State struct is a structure used by the quasi e-greedy policy with target
% recovery
stateStruct.CurrState=state;
stateStruct.PrevState=SARSAinStruct.state_prev;

% computation of the action following one of the possible e-greedy policies
[action,ExplorationFlag] = epsilonGreedy(epsilon,stateStruct,Q);

% update of the Q matrix
Q(state_prev,action_prev)=Q(state_prev,action_prev)...
    +alpha*(reward+gamma*Q(state,action)-Q(state_prev,action_prev));

% output structure creation
SARSAoutStruct.Q=Q;
SARSAoutStruct.state=state;
SARSAoutStruct.action=action;
SARSAoutStruct.reward=reward;
SARSAoutStruct.ExplorationFlag=ExplorationFlag;

end

