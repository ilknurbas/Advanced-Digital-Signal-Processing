% Performs 2D Fourier filtering by cropping the image
function F = fourier_filtering(path, flag)

    F = fftshift(fft2(imread(path)));
    F_abs = abs(F);  % in freq domain
    [col, row] = size(F);
    figure;
    subplot(3, 3, 1);
    imshow(log(1 + F_abs), []);  
    title('0. Log Magnitude Spectrum');
        
    % zero out all terms expect the second one that contains object's wavefront
        
    % 1. zero out left
    F(:, 1:(row/2)) = 0;
    % figure;
    subplot(3, 3, 2);
    imshow(log(1 + abs(F)), []);  
    title('1. Left Region Zeroed in Freq Domain');
        
    % 2. zero out upper
    F(1:col/2, :) = 0;
    % figure;
    subplot(3, 3, 3);
    imshow(log(1 + abs(F)), []);  
    title('2. Top Region Zeroed in Freq Domain');
        
    % 3. zero out center vertical strip
    [~, y] = meshgrid(1:row, 1:col);
    square = min(row, col) * 0.2; % visually selected
     if ~flag
        square = min(row, col) * 0.05; 
     end
    % mask = (abs(x - row/2) <= square/2) & (abs(y - col/2) <= square/2);
    mask = (abs(y - col/2) <= square/2);
    F(mask) = 0 ; 
    % figure;
    subplot(3, 3, 4);
    imshow(log(1 + abs(F)), []);  
    title('3. Vertical Strip Region Zeroed in Freq Domain');
        
    % 4. move to center
    [~, idx] = max(abs(F(:)));
    [term_row, term_col] = ind2sub(size(F), idx);
    F = circshift(F, [col/2 - term_row, row/2 - term_col]);
    % figure;
    subplot(3, 3, 5);
    imshow(log(1 + abs(F)), []); 
    title('4. Centering Max Frequency Component');
        
    % 5. get center
    [x, y] = meshgrid(1:row, 1:col);
    square = min(row,col) * 0.20; % visually selected
    if ~flag
        square = min(row,col) * 0.15; % visually selected
    end
    mask = (abs(x - row/2) <= square/2) & (abs(y - col/2) <= square/2);
    F(~mask) = 0 ; 
    % figure;
    subplot(3, 3, 6);
    imshow(log(1 + abs(F)), []);  
    title('5. Preserve Center Square Region');                    
    
end

