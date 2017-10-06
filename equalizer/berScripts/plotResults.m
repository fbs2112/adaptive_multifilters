clear;
clc;
close all;




addpath(['.' filesep 'results']);



fileVector = {'resultsBER1','resultsBER2'};

for k = 1:length(fileVector)

    load([fileVector{k} '.mat']);
    figure;

    h = plot(randn(10,1));
    availableMarkers = set( h, 'Marker' );

    for j = 1:size(ber,1)


        for i = 1:size(ber,2)
            h = semilogy(SNR,squeeze(ber(j,5,:)));
            h.Marker = availableMarkers{i};
            hold on

        end

    end

    xlabel('SNR [dB]','interpreter','latex');
    ylabel('BER','interpreter','latex');

    H = legend('$\bar{\gamma}_{\mathrm{NL}} = \bar{\gamma}_{\mathrm{L}}$','$\bar{\gamma}_{\mathrm{NL}} = 1.5\bar{\gamma}_{\mathrm{L}}$','$\bar{\gamma}_{\mathrm{NL}} = 2\bar{\gamma}_{\mathrm{L}}$',...
        '$\bar{\gamma}_{\mathrm{NL}} = 2.5\bar{\gamma}_{\mathrm{L}}$','$\bar{\gamma}_{\mathrm{NL}} = 3\bar{\gamma}_{\mathrm{L}}$','$\bar{\gamma}_{\mathrm{NL}} = 3.5\bar{\gamma}_{\mathrm{L}}$',...
        '$\bar{\gamma}_{\mathrm{NL}} = 4\bar{\gamma}_{\mathrm{L}}$');

    set(H,'interpreter','latex','location','SouthWest')
    xlim([min(SNR) max(SNR)]);
    ylim([1e-8 1]);
    grid on
end







rmpath(['.' filesep 'results']);
