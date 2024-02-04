clear all
close all

num_colors = 64;  % Number of colors
cmap = interp1([1, num_colors/2, num_colors], [0, 0, 1; 1, 1, 0; 1, 0, 0], 1:num_colors);  % Gradient from red to yellow to blue

DAMAS = 0; % Draw DAMAS or beamforming

f0=[200 250 315 400 500 630 800 1000 1250 1600 2000];
tt=7;
m=51;n=51;l=100;w=100;Z=75;
Y0=72; % Distance from microphone array center to rotation center plane
alpha=0/180*pi; % Azimuth angle, in radians, counterclockwise as positive
ii=1:m;jj=1:n;
x1(ii)=(ii-(m+1)/2)*l/(m-1);
z1(jj)=25.5+(jj-1)*w/(n-1);
[x2,z2]=meshgrid(x1,z1);
x=x2'*cos(alpha);y=Y0+x*tan(alpha);z=z2'; % x, y, z are absolute coordinates, calculating coordinates; x2' and z2' are plotting coordinates
% Grid relative to hub relative coordinates, calculating coordinates
xw=x;yw=y-Y0;zw=z-Z;
R=44; % Impeller rotation radius
% Impeller outline, plotting coordinates
kk=1:360;
xo(kk)=R*cos((kk-1)/180*pi);
zo(kk)=R*sin((kk-1)/180*pi)+Z;
step1=num2str(f0(tt));
if DAMAS == 0
    foldername{1,1}=strcat(step1,'Hz_1.xlsx');
    imag_title = ['\bf{Beamforming (no diagonals) ' num2str(f0(tt)) ' Hz}'];
else
    foldername{1,1}=strcat(step1,'Hz_2_400iterations.xlsx');
    imag_title = ['\bf{DAMAS ' num2str(f0(tt)) ' Hz 400 iterations}'];
end
WM=xlsread(foldername{1,1});
figure(1);title(imag_title);hold on;
[im,jm]=find(WM==max(WM(:)));
imax(1)=(im-(m+1)/2)*l/(m-1);jmax(1)=25.5+(jm-1)*w/(n-1);
Rmax(1)=sqrt(imax(1)*imax(1)+(jmax(1)-Z)*(jmax(1)-Z));
% plot(x(:,j),z(:,j),'k','linewidth',1.5);
% plot(x(:,mod(j+n/3,n)),z(:,mod(j+n/3,n)),'k','linewidth',1.5);
% plot(x(:,mod(j+2*n/3,n)),z(:,mod(j+2*n/3,n)),'k','linewidth',1.5);
summit=abs(max(WM(:)));
% v=[summit-3:0.1:summit];
%v=[summit-12:0.2:summit-6,summit-5.8:0.1:summit];
v=summit-12:0.1:summit;
contour(x,z,WM,v);hold on;axis equal
colormap(cmap);
plot(xo,zo,'k-');
%temp1=caxis;

if DAMAS == 0
    file_name = ['Beamforming_Nodiagnal_' num2str(f0(tt)) '_Hz'];
    print(gcf,file_name,'-dpng','-r300') % 300 dpi, i.e. 300 dots per inch
else
    file_name = ['DAMAS_' num2str(f0(tt)) '_Hz'];
    print(gcf,file_name,'-dpng','-r300') % 300 dpi, i.e. 300 dots per inch
end
