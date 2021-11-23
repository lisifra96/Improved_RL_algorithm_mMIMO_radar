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

% This script plots various figures relative to RL algorithms from the data evaluated in main

% Plot of the average reward
Plot_Reward

% Plot of the average relative Q difference
Plot_Qdiff

% Plot of the Frobenius norm of Q
Plot_Qnorm

% Plot of the epsilon value in the varying epsilon case
if VaryingEpsilonFlag
    Plot_epsilonTemporalEvolution
end

% Plot of the alphan value in the varying alpha case
if VaryingAlphaFlag
    Plot_alphaTemporalEvolution
end