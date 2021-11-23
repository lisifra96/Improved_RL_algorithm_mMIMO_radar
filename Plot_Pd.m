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

% This script plots the Pd of the targets

switch TargetScenario
    
    case 1
        
        for u=1:NumberVaryingScenarios
            for s=1:length(NuTargetIndex{u})
                                
                figure(nFigShared)
                plot(time_array,DetectionFrequency(NuTargetIndex{u}(s),time_array),'Color',Colour,'LineStyle',Marker,'LineWidth',line_width,'DisplayName',PdLegend);
                legend('off')
                legend('show')
                legend({},'FontSize',FontSize)
                grid on;
                hold on;
                axis([1, Tmax, 0, 1]);
                xlabel('Discrete time k','FontSize',FontSize);
                ylabel('P_D','FontSize',FontSize);
                title(sprintf('Detection Frequency of the target corresponding to \\nu=%.2f',nu_array(NuTargetIndex{u}(s))),'FontSize',FontSize);
                set(gca,'FontSize',FontSize)
                nFigShared=nFigShared+1;
                
            end
        end
        
    case 2
        
        NumberOfTargets=length(NuTargetIndex{1});
        for s=1:NumberOfTargets
            
            figure(nFigShared)
            plot(time_array,DetectionFrequency(NuTargetIndex{1}(s),time_array),'Color',Colour,'LineStyle',Marker,'LineWidth',line_width,'DisplayName',PdLegend);
            legend('off')
            legend('show')
            legend({},'FontSize',FontSize)
            grid on;
            hold on;
            axis([1, Tmax, 0, 1]);
            xlabel('Discrete time k','FontSize',FontSize);
            ylabel('P_D','FontSize',FontSize);
            title(sprintf('Detection Frequency of the target corresponding to \\nu=%.2f',nu_array(NuTargetIndex{1}(s))),'FontSize',FontSize);
            set(gca,'FontSize',FontSize)
            nFigShared=nFigShared+1;
            
        end
        
    case 3
        
        NumberOfTargets=length(NuTargetIndex{1});
        for s=1:NumberOfTargets
            
            figure(nFigShared)
            plot(time_array,DetectionFrequency(NuTargetIndex{1}(s),time_array),'Color',Colour,'LineStyle',Marker,'LineWidth',line_width,'DisplayName',PdLegend);
            legend('off')
            legend('show')
            legend({},'FontSize',FontSize)
            grid on;
            hold on;
            xlabel('Discrete time k','FontSize',FontSize);
            ylabel('P_D','FontSize',FontSize);
            axis([1, Tmax, 0, 1]);
            title(sprintf('Detection Frequency of the target corresponding to \\nu=%.2f',nu_array(NuTargetIndex{1}(s))),'FontSize',FontSize);
            set(gca,'FontSize',FontSize)
            nFigShared=nFigShared+1;
            
        end
        
end
