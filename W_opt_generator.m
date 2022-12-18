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


% This script computes all the W matrix that solve the optimization problem
% for all different possible steering vectors with Nt=Nr fixed and number of
% angular bins fixed

nuBinsNumber=10;                                    % Number of ν bins 
nu_array=-0.5+(0:nuBinsNumber-1)'/nuBinsNumber;     % Array defining the grid of ν=d/λ*sin(θ) values
Nt_array=[1e2];                                     % Array containing the number of transmit antennas
loopNt=length(Nt_array);
Ptot=1;                                             % Total transmit power
MaxDetectableTarget=3;                              % Max number of target that can be detected
CombosCell=cell(MaxDetectableTarget+1,1);           % cell containing all possible combinations of n elements over a set with nuBinsNumber elements
CombosCell{1}=0;                                    % The Combos Cell in the case of no target present is seto to 0 by default
PathName=what('Wcubes');
PathName=PathName.path;

%% l=1 case

l=1;
tic
Nt=Nt_array(l);
aT_Mat=exp(1i*2*pi*(0:Nt-1)'*nu_array');            % each column correspond to the steering vector for the associated ν

UpperBoundPartialSum=ceil(nchoosek(nuBinsNumber,MaxDetectableTarget)*(nuBinsNumber-MaxDetectableTarget+1)/(nuBinsNumber-2*MaxDetectableTarget+1));
W_CUBE=zeros(Nt,Nt,UpperBoundPartialSum);           % The number of layers is set to a higher value because there isn't a closed form to compute the sum of binomial coefficients
W_CUBE(:,:,1)=sqrt(Ptot/Nt)*eye(Nt);                % In case there are no detected target the W matrix is the identity

BinomialCoeffSum=ones(1,MaxDetectableTarget+1);     % This is a vector containing the sum of the binomial coefficients up to MaxDetectableTarget

% cycle to create the W cube

for n=1:MaxDetectableTarget
    combos = combntns(1:nuBinsNumber,n);                % matrix containing all the combinations of n elements out of nuBinsNumber
    CombosCell{n+1}=combos;
    Ncombos=nchoosek(nuBinsNumber,n);                   % Binomial(nuBinsNumber,n) 
    ThisBinomialCoeffSum=BinomialCoeffSum(n);
    parfor m=1:Ncombos
        A=aT_Mat(:,combos(m,:));                        % Matrix conaining all the steering vectors associated to the given combination
%         [~, W]=Alg2v1(Nt,A,Ptot);      % optimal W computation (CVX)
        W=Closed_Form_W(A,Ptot);                    % optimal W computation (Closed-form)
        W_CUBE(:,:,ThisBinomialCoeffSum+m)=W;
    end
    BinomialCoeffSum(n+1)=BinomialCoeffSum(n)+Ncombos;
end

W_CUBE=W_CUBE(:,:,1:BinomialCoeffSum(end));             % reduction of W_CUBE dimension to the real one

WcubeStruct.Wcube=W_CUBE;
WcubeStruct.BinomialCoeffSum=BinomialCoeffSum;
WcubeStruct.CombosCell=CombosCell;

FileName=['Wcube_nuBins_',num2str(nuBinsNumber),'_Nt_',num2str(Nt),...
    '_MaxTarget_',num2str(MaxDetectableTarget),'_Ptot_',num2str(Ptot)];
save([PathName '/' FileName],...
    'WcubeStruct','-v7.3');
toc

%% l>1 case (to avoid the computation of combos and BinomialCoefficientSum for every Nt)

if loopNt>1
    for l=2:loopNt
        tic
        Nt=Nt_array(l);
        aT_Mat=exp(1i*2*pi*(0:Nt-1)'*nu_array');            % each column correspond to the steering vector for the associated ν

        W_CUBE=zeros(Nt,Nt,UpperBoundPartialSum);           % The number of layers is set to a higher value because there isn't a closed form to compute the sum of binomial coefficients
        W_CUBE(:,:,1)=sqrt(Ptot/Nt)*eye(Nt);                % In case there are no detected target the W matrix is the identity

        % cycle to create the W cube

        for n=1:MaxDetectableTarget
            combos = CombosCell{n+1};                             % matrix containing all the combinations of n elements out of nuBinsNumber
            Ncombos=nchoosek(nuBinsNumber,n);                     % Binomial(nuBinsNumber,n) 
            ThisBinomialCoeffSum=BinomialCoeffSum(n);
            parfor m=1:Ncombos
                A=aT_Mat(:,combos(m,:));                    % Matrix conaining all the steering vectors associated to the given combination
%                 [~, W]=Alg2v1(Nt,A,Ptot);      % optimal W computation (CVX)
                W=Closed_Form_W(A,Ptot);                    % optimal W computation (Closed-form)
                W_CUBE(:,:,ThisBinomialCoeffSum+m)=W;
            end
        end

        W_CUBE=W_CUBE(:,:,1:BinomialCoeffSum(end));             % reduction of W_CUBE dimension to the real one

        WcubeStruct.Wcube=W_CUBE;
        
        FileName=['Wcube_nuBins_',num2str(nuBinsNumber),'_Nt_',num2str(Nt),...
            '_MaxTarget_',num2str(MaxDetectableTarget),'_Ptot_',num2str(Ptot)];
        save([PathName '/' FileName],...
            'WcubeStruct','-v7.3');
        toc
    end
end

clear nuBinsNumber nu_array Nt_array loopNt Ptot MaxDetectableTarget CombosCell...
    l Nt aT_Mat W_CUBE BinomialCoeffSum n combos Ncombos ThisBinomialCoeffSum m...
    A W UpperBoundPartialSum;
