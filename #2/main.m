% Main script

clear;
clc;
close all;

% object wavefront reconstruction
wavefront_10 = reconstruct('object10pixels.png', true); 
wavefront_14 = reconstruct('object14pixels.png', false); %

%%
% Comparison of amplitude and phase values
figure;
subplot(1, 2, 1);
plot(mean(abs(wavefront_10),1), 'g');  
hold on;
plot(mean(abs(wavefront_14),1), 'r');  
hold off;
title('Crosssection Mean of Amplitude values');
xlabel('Spatial Position');
ylabel('Amplitude');
legend('10px', '14px');

subplot(1, 2, 2);
plot(mean(angle(wavefront_10),1), 'g');  
hold on;
plot(mean(angle(wavefront_14),1), 'r');  
hold off;
title('Crosssection Mean of Phase values');
xlabel('Spatial Position');
ylabel('Phase');
legend('10px', '14px');

 