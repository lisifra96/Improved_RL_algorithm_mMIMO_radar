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

function [W] = getWfromTargetIndexes_offline(WcubeStruct,TargetIndexes)
% This function extracts the W matrix from Wcube corresponding to the
% optimization problem for the target present in TargetIndexes. 
%
% To understand the structure of Wcube an explicative example is provided.
% Suppose that nuBinsNumber=5 and MaxDetectableTargets=2. Then Wcube will
% contain nchoosek(5,0)+nchoosek(5,1)+nchoosek(5,2)=16 elements corresponding
% to all the possible combinations of angular bins up to
% MaxDetectableTargets=2 out of nuBinsNumber=5. The first layer of Wcube
% corresponds to W_ort, i.e. the solution associated to the case n=0, where
% n is the number of angular bins in the set Omega (see article). The
% following nchoosek(5,1)=5 layers correspond to the solutions of the
% optimization problem associated to each possible combination of n=1
% targets' angular position, whose order is specified in CombosCell{n+1=2}.
% The following nchoosek(5,2)=10 layers correspond to the solution of the
% optimization problem associated to each possible combination of n=2
% targets' angular position, whose order is specified in CombosCell{n+1=3}.
% If the user wants to find the solution associated to a couple of angular
% bins, he can compute the corresponding Wcube index as index_1+index_2.
% index_1 is the index associated to the first solution with n=2, in the
% example this corresponds to nchoosek(0,5)+nchoosek(1,5)+1=7. index_2
% corresponds to the index associated to the specific choice of the two
% angular bins and can be retrieved by comparing the vector containing the
% two angular bins with the vector CombosCell{n+1=3}.
%
% INPUT VARIABLES
% WcubeStruct.Wcube:            Cube containing the previously computed solution of the
%                               optimization problem for every possible target angle
%                               combination (see W_opt_generator). Each
%                               layer corresponds to one solution.
% WcubeStruct.BinomialCoeffSum: Vector containing the cumulative sum of the binomial coefficients
%                               up to MaxDetectableTargets. 
% WcubeStruct.CombosCell:       Cell containing all the combinations of n elements out
%                               of nuBinsNumber
% TargetIndexes:                Array containing the indexes of the angular
%                               bins in the set Omega
% OUTPUT VARIABLE
% W:                            Optimal weight matrix 

NumberTarget=length(TargetIndexes);
TargetIndexes=reshape(sort(TargetIndexes,'ascend'),1,NumberTarget);
StartIndex=[1, WcubeStruct.BinomialCoeffSum+1];

% These instructions set the TargetIndexes array to [0] if it is an empty
% set, while leaves it unchanged if is not empty
TargetIndexes=[TargetIndexes, 0];
TargetIndexes=TargetIndexes(1:end-(TargetIndexes(1)~=0));


index_1=StartIndex(NumberTarget+1);
[~,index_2]=max(sum(WcubeStruct.CombosCell{NumberTarget+1}==TargetIndexes,2));
W_Index=index_1+index_2-1;
W=WcubeStruct.Wcube(:,:,W_Index);

end
