function S_norm = plot_normalized_PSD2(rho)
%   This function plots the normalized PSD of an AR process whose
%   coefficients are stored in rho

deltaf=1e-3;
f = -0.5:deltaf:0.5;
den=0;

for i = 2:length(rho)
    den = den + rho(i)*exp(-1j*2*pi*(i-1)*f);
end
S = 1./abs(1+den).^2;
S_norm = S/max(S);

%% Plot
line_width=2;
FontSize=26;

figure
plot(f,(S_norm),'LineWidth',line_width)
xlabel('$\nu$','interpreter','latex','FontSize',FontSize,'interpreter','latex')
ylabel('$\frac{S_c(\nu)}{\max{S_c(\nu)}}$','interpreter','latex','FontSize',FontSize)
grid on
%title('\textbf{Normalized PSD}','FontSize',FontSize,'interpreter','latex')
set(gca,'FontSize',FontSize)

end

