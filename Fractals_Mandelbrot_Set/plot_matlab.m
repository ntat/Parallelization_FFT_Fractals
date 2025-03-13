clear all
%close all
clc
fid = fopen('test', 'r')
zm = fread(fid, 'int32');
A = reshape(zm,4096,4096);
imagesc(log(A))

%colormap( [jet();flipud( jet() );0 0 0] );
%print(gcf,'foo3.png','-dpng','-r2000');