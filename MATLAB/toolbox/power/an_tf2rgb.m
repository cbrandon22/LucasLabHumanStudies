function an_tf2rgb

%Sample power spectral plot


% Say this is the given matrix:
G = rand(10,10)*300 - 100;
% Use IMAGESC to plot G.
figure('pos',[100 100 1000 800]);
% colormap(copper) % You realize this affects final image (H)?
subplot(1,2,1)
imagesc(G);
title('IMAGESC (MxN)')
% Now make an RGB image that matches display from IMAGESC:
C = colormap;  % Get the figure's colormap.
L = size(C,1);
% Scale the matrix to the range of the map.
Gs = round(interp1(linspace(min(G(:)),max(G(:)),L),1:L,G));
H = reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.
subplot(1,2,2)
image(H)  % Does this image match the other one?
title('IMAGE (MxNx3)')