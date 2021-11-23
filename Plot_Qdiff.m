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

% Plot of relative Q difference

figure(nFigSingle)
plot(time_array, RelErrQ_average,Colour,'LineWidth',line_width);
grid on;
xlabel('Discrete time k','FontSize',FontSize);
ylabel('$\frac{\big\vert\big\vert\textbf{Q}_k-\textbf{Q}_{k-1}\big\vert\big\vert}{\big\vert\big\vert\textbf{Q}_k\big\vert\big\vert}$','fontsize',FontSize,'interpreter','latex');
title(sprintf('Q relative difference, %s = %s',LoopVar,LoopVarValue),'FontSize',FontSize);
set(gca,'FontSize',FontSize)
nFigSingle=nFigSingle+1;

