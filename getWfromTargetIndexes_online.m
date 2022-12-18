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

function [W] = getWfromTargetIndexes_online(Nt,A,Ptot,W_ort)
% This function computes the W matrix corresponding to the target indexes
% contained in TargetIndexes.
%
% INPUT VARIABLES
% A:            Matrix whose columns correspond to the transmit steering 
%               vector associate to each angular bin in the set omega
% Nt:           Number of transmitters
% Ptot:         total transmitted power
% Wort:         sqrt(Ptot/Nt)*eye(Nt)
% OUTPUT VARIABLE
% W:            Optimal weight matrix 

if isempty(A)
    W=W_ort;
else
    %[~,W]=Alg2v1(Nt,A,Ptot);
    W=Closed_Form_W(A,Ptot);
end

end