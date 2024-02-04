clc;
clear all
close all
tic;

num_colors = 64;  % Number of colors
cmap = interp1([1, num_colors/2, num_colors], [0, 0, 1; 1, 1, 0; 1, 0, 0], 1:num_colors);  % Gradient from blue to yellow to red

t1=[0.06 0.55 0.06 0.55];t2=[0.55 0.55 0.05 0.05];t3=[0 100 1000 5000];
set(gcf,'position',[500 200 520 520])

% Allocate space for the large matrix
l=1.52;w=1.52;% Length and width of the grid
m=77;n=77;% Dimensions of the grid in x and y directions; must be larger than the dimensions of the airfoil
p=32; % Number of microphones
N0=m*n;% Total number of grid points
A=zeros(N0);
dx=zeros(m,n,p);
dy=zeros(m,n,p);
dz=zeros(m,n,p);
rnl2(m,n,p)=0;% Distance squared between microphone and grid point
tn(m,n,p)=0;% Propagation time of sound between microphone and grid point
rnel(m,n,p)=0;% Distance between microphone and grid point
e=zeros(m,n,p);
et=zeros(N0,p);
X=zeros(1,N0);
Y=zeros(1,N0);
Y1=zeros(1,N0);
X1=zeros(1,N0);
XM=zeros(m,n);
x(m,n)=0; % Initial value of grid x coordinate
y(m,n)=0; % Initial value of grid y coordinate

% Initial parameter settings
NI=5000;% Total number of iterations
c0=344;% Speed of sound
Urr=[0 0 0];% Incoming flow velocity vector
M=Urr/c0;% Mach number
beta=1-M*M';
f=10000;% Frequency
z0=1.5;% Distance between grid center and array center;
% Array radius (to be modified)
R1=0.1;

% Grid coordinates
ii=1:m;jj=1:n;
x1(ii)=(ii-(m+1)/2)*l/(m-1);
y1(jj)=(jj-(n+1)/2)*w/(n-1);
[x2,y2]=meshgrid(x1,y1);
x=x2';y=y2';

% Microphone coordinates (microphones may not yet be aligned with the airfoil, which is tragic)
sita=2*pi/p;% Angle interval between microphones
ii=1:p;
xm(ii)=R1*cos(sita*(ii-1));ym(ii)=R1*sin(sita*(ii-1));% Coordinates of the microphones

fid=fopen('x.dat','wt');% Write file path
matrix=x;               
[m,n]=size(matrix);
 for ii=1:1:m
   for jj=1:1:n
      if jj==n
        fprintf(fid,'%g\n',matrix(ii,jj));
     else
       fprintf(fid,'%g\t',matrix(ii,jj));
      end
   end
end
fclose(fid);

fid=fopen('y.dat','wt');% Write file path
matrix=y;                        
[m,n]=size(matrix);
 for ii=1:1:m
   for jj=1:1:n
      if jj==n
        fprintf(fid,'%g\n',matrix(ii,jj));
     else
       fprintf(fid,'%g\t',matrix(ii,jj));
      end
   end
end
fclose(fid);

% Calculate distance vector, distance, propagation time from grid point to microphone
for ii=1:m
    for jj=1:n
        for kk=1:p
        dx(ii,jj,kk)=xm(kk)-x(ii,jj);
        dy(ii,jj,kk)=ym(kk)-y(ii,jj);
        dz(ii,jj,kk)=z0;
        rn=[dx(ii,jj,kk);dy(ii,jj,kk);dz(ii,jj,kk)];      
        rnl2(ii,jj,kk)=rn'*rn;% Calculate the square of the distance from grid to microphone
        tn(ii,jj,kk)=1/beta*(-M*rn+((M*rn)^2+beta*rnl2(ii,jj,kk))^(1/2))/c0;
        rne=[dx(ii,jj,kk);dy(ii,jj,kk);dz(ii,jj,kk)]-Urr'*tn(ii,jj,kk);
        rnel(ii,jj,kk)=c0*tn(ii,jj,kk)+M*rne;
        end
    end
end

% Calculate control function
for kk=1:p
    for ii=1:m       
        e(ii,:,kk)=exp(2*i*pi*f*tn(ii,:,kk)).*rnel(ii,:,kk)/z0;
    end
end

% Calculate control vector
oo=1;
for jj=1:n
    for ii=1:m
        et(oo,1:p)=e(ii,jj,:);
        oo=oo+1;
    end
end

% Calculate the A matrix, the most important step in the algorithm
for ii=1:N0
    for jj=1:N0
        A(ii,jj)=(conj(et(ii,:))*(1./et(jj,:))'*((1./et(jj,:))*et(ii,:).'))/(p*p);        
    end
end

% Simulate the position of the sound source
data=xlsread('xy data.xls')/100;
data(:,1)=data(:,1)-0.74;
i=1:length(data(:,1));
X((data(i,2)/w*(n-1)+(n-1)/2)*m+data(i,1)/l*(m-1)+(m+1)/2)=4;% Initial value means the square of effective sound pressure

% Calculate the Y of each grid point
for oo=1:N0
    Y(oo)=A(oo,1:N0)*X(1:N0)';
end

% Convert the one-dimensional array to a two-dimensional array
oo=1;
for jj=1:n
    for ii=1:m
        YM(ii,jj)=Y(oo);
        oo=oo+1;
    end
end

% Convert the energy unit of the grid points to decibels
for ii=1:m
    W(ii,:)=10*log10(YM(ii,:)/4*10^10);
end

% Find the maximum position
subplot(2,2,1);
set(gca,'ytick', [-0.5 0 0.5])
title(['\bf{Beamforming}']);hold on;
[im,jm]=find(W==max(W(:)));
imax(1)=x(im,jm);jmax(1)=y(im,jm);

% Draw the contour plot
summit=abs(max(W(:)));
v=summit-40:0.08:summit;
contour(x,y,W,v);hold on;axis equal;axis([-0.8,0.8,-0.8,0.8]);
colormap(cmap);

fid=fopen('W1.dat','wt');% Write file path
matrix=W;                        
[m,n]=size(matrix);
 for ii=1:1:m
   for jj=1:1:n
      if jj==n
        fprintf(fid,'%g\n',matrix(ii,jj));
     else
       fprintf(fid,'%g\t',matrix(ii,jj));
      end
   end
end
fclose(fid);

X=Y;% Initial conditions
ff=2;% Initialize plot number
for ite=1:1:NI
    X1=X;
    X(1)=Y(1)-A(1,2:N0)*X1(2:N0)';
    if X(1)<0
        X(1)=0;
    end
    for ii=2:N0-1
        X(ii)=Y(ii)-A(ii,1:ii-1)*X(1:ii-1)'-A(ii,ii+1:N0)*X1(ii+1:N0)';
        if X(ii)<0
            X(ii)=0;
        end
    end
    X(N0)=Y(N0)-A(N0,1:N0-1)*X(1:N0-1)';
    if X(N0)<0
        X(N0)=0;
    end    
    X1=X;
    for ii=N0-1:2
        X(ii)=Y(ii)-A(ii,1:ii-1)*X(1:ii-1)'-A(ii,ii+1:N0)*X1(ii+1:N0)';
        if X(ii)<0
            X(ii)=0;
        end        
    end
    X(1)=Y(1)-A(1,2:N0)*X1(2:N0)';
    if X(1)<0
        X(1)=0;
    end
    
    % Display every Ni steps
    if ite==100||ite==1000||ite==5000
        ite
        
        % Convert the one-dimensional array to a two-dimensional array
        oo=1;
        for jj=1:n
            for ii=1:m
                XM(ii,jj)=X(oo);
                oo=oo+1;
            end
        end
        
        % Convert the energy unit of the grid points to decibels
        for ii=1:m
            for jj=1:n
                if XM(ii,jj)==0
                    W(ii,jj)=0; 
                else
                    W(ii,jj)=10*log10(XM(ii,jj)/4*10^10);
                end
            end
        end
        
        % Plot
        subplot(2,2,ff);
        set(gca,'ytick', [-0.5 0 0.5])
        title(['\bf{DAMAS  ' num2str(ite) ' Iterations}']);
        hold on;
        [im,jm]=find(W==max(W(:)));
        imax(ff)=x(im,jm);jmax(ff)=y(im,jm);
        % Draw the contour plot
        summit=abs(max(W(:)));
        v=summit-40:0.08:summit;
        contour(x,y,W,v);hold on;axis equal;axis([-0.8,0.8,-0.8,0.8]);
        colormap(cmap);
        
        tmp=strcat('W',num2str(ff),'.dat');
        fid=fopen(tmp,'wt');% Write file path
        matrix=W;                        
        [m,n]=size(matrix);
         for ii=1:1:m
           for jj=1:1:n
              if jj==n
                fprintf(fid,'%g\n',matrix(ii,jj));
             else
               fprintf(fid,'%g\t',matrix(ii,jj));
              end
           end
        end
        fclose(fid);        
        ff=ff+1;
    end
end

toc;

h = colorbar('Location', 'southoutside');
set(h, 'Position', [0.1, 0.045, 0.8, 0.015]);

print(gcf,'UCAS_10000','-dpng','-r300')
