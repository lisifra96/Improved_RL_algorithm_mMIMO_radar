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

% Script to plot the alpha value temporal evolution for the varying
% alpha algorithm

figure(nFigSingle)
plot(time_array, alphaValuesArrayAvg, 'r','LineWidth',line_width);
axis([1, Tmax, AlphaMin, AlphaMax]);
grid on;
xlabel('Discrete time k','FontSize',FontSize);
ylabel('$\alpha_k$','fontsize',18,'interpreter','latex');
title('$\textbf{$\alpha$ temporal evolution}$','FontSize',FontSize,'interpreter','latex');
set(gca,'FontSize',FontSize)
nFigSingle=nFigSingle+1;
