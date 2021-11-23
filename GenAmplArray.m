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

function [AmplitudeArray] = GenAmplArray(NuTargetIndex,TargetEffAmpl,nuBinsNumber)
% This function computes at each time instant an array of simulated
% swerling 0 target amplitudes according to the scenario.
% INPUT VARIABLES
% NuTargetIndex:    Row vector containing the angular bins of the present target    
% TargetEffAmpl:    Row vector containing the associated effective
%                   amplitude, i.e. the square root of the power
% nuBinsNumber:     Number of bins of the ν array
% OUTPUT VARIABLE
% AmplitudeArray:   Column vector containing nuBinsNumber element that are
%                   the amplitudes of the target in the bins with index
%                   NuTargetIndex and 0 in the bins where there is no
%                   target

TargetNumber=length(NuTargetIndex);
AmplitudeArray=zeros(nuBinsNumber,1);
if TargetNumber~=0
    AmplitudeArray(NuTargetIndex)=reshape(TargetEffAmpl,TargetNumber,1).*exp(1i*2*pi*rand(TargetNumber,1));
end

end

