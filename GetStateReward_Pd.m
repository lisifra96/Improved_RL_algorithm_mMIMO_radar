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

function [state,reward] = GetStateReward_Pd(GetStateRewardInputStruct)
% This function computes the state and the reward at a given iteration
% based on the decision statistic of the robust Wald-type test.
% INPUT VARIABLES
% GetStateRewardInputStruct.WaldTest_Stat:      Array containing the decision statistic in each angular bin
% GetStateRewardInputStruct.DetectionArray:     double(WaldTest_Stat>Thresh)  
% GetStateRewardInputStruct.MaxTargetNumber:    Maximum number of detectable targets
% GetStateRewardInputStruct.Thresh:             Wald Test threshold
% OUTPUT VARIABLES
% state:            MDP state
% reward:           MDP reward

WaldTest_Stat=GetStateRewardInputStruct.WaldTest_Stat;
MaxTargetNumber=GetStateRewardInputStruct.MaxTargetNumber;
NoverThresh=sum(GetStateRewardInputStruct.DetectionArray);                              % Number of detected targets

Pd_approx=marcumq(sqrt(WaldTest_Stat),sqrt(GetStateRewardInputStruct.Thresh));          % computation of the estimate of the Pd
state=min(MaxTargetNumber,NoverThresh)+1;                                               % MDP state computation
[~,ActionIndex] = maxk(WaldTest_Stat,MaxTargetNumber);                                  % ActionIndex contains the indexes of the MaxTargetNumber highest values in WaldTest_Stat 
ActionIndex2=ActionIndex(1:(state-1));                                                  % ActionIndex2 contains the indexes of the (state-1) highest values in WaldTest_Stat

reward=sum(Pd_approx(ActionIndex2))-sum(Pd_approx(setdiff(ActionIndex,ActionIndex2)));  % Computation of the reward

end