% Performs 2D Fourier filtering and Wavefront Propagation respectively
function [angular_spectrum] = reconstruct(path,flag)
    % in meters
    pixel_size = 3.45E-6;
    dist = 0.11; % ray and camera
    wavelength = 532E-9;

    disp(['path: ', path]);
    F = fourier_filtering(path, flag);

    U0 = ifft2(ifftshift(F)); % object wavefront in spatial domain
    
    subplot(3, 3, 7:9);
    imshow(log(1 + abs(U0)), []);
    title('6. Object Wavefront in Spatial Domain');
    sgtitle(['2D Fourier Filtering for ',path]);

    distances = -0.4:0.02:-0.2; % back propagation
    distances = -0.28; % optimal value is choosen visually

    for i = 1:length(distances)
   
        distance = distances(i);
        angular_spectrum = AngularSpectrum(U0, distance, wavelength, pixel_size);
 
        figure;
        subplot(1, 2, 1);
        imagesc(abs(angular_spectrum));
        colormap(gray);
        colorbar;
        title('Amplitude of Reconstructed Object Wavefront');
        
        subplot(1, 2, 2);
        imagesc(angle(angular_spectrum));
        colormap(gray);
        colorbar;
        title('Phase of Reconstructed Object Wavefront');
        sgtitle(['Reconstruction of ', path, ' when distance is ',num2str(distance)]);

    end
end

