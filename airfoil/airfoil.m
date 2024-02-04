% Final program for wind tunnel simulation, considering uniform flow.
clc;
clear all
close all
tic;

num_colors = 64;  % Number of colors
cmap = interp1([1, num_colors/2, num_colors], [0, 0, 1; 1, 1, 0; 1, 0, 0], 1:num_colors);  % Interpolate colors
low = 70; high = 100;

t1=[0.05 0.55 0.05 0.55];t2=[0.55 0.55 0.05 0.05];t3=[0 100 1000 5000];

% Allocate space for large matrices
l=0.6;w=0.26;% Length and width of the grid;
m=61;n=27;% Dimensions of the grid in x and y directions; must be larger than the dimensions of the airfoil
p=32; % Number of microphones
N0=m*n;% Total number of grid points
A=zeros(N0);
dx=zeros(m,n,p);
dy=zeros(m,n,p);
dz=zeros(m,n,p);
rnl2(m,n,p)=0;% Distance squared between microphone and grid points
tn(m,n,p)=0;% Propagation time of sound between microphone and grid points
rnel(m,n,p)=0;% Distance between microphone and grid points
e=zeros(m,n,p);
et=zeros(N0,p);
X=zeros(1,N0);
Y=zeros(1,N0);
Y1=zeros(1,N0);
X1=zeros(1,N0);
XM=zeros(m,n);
x(m,n)=0; % Initial value of grid x coordinates
y(m,n)=0; % Initial value of grid y coordinates

% Initial parameter settings
NI=5000;% Total number of iterations
xss=0.15*ones(1,15);yss=[-0.07:0.01:0.07];% Simulated sound source
c0=344;% Speed of sound 
Urr=[0 0 0];% Inflow velocity vector
M=Urr/c0;% Mach number at the wind tunnel outlet
beta=1-M*M';
f=2000;% Frequency
l1=0.074;w1=0.15;% Chord length and wind tunnel diameter of the airfoil;
s=8;t=16;% Dimensions of the airfoil in x and y directions;
l2=0.012;w2=0.0048;% Length and width of the sawtooth;
d=13;b=31;% Dimensions in the x and counterclockwise directions of the sawtooth (to be calculated);
x0=0.069;% Distance from the airfoil leading edge to the wind tunnel outlet;
xc1=0.16;% Distance from the airfoil trailing edge to the wind tunnel outlet;
z0=1.5;% Distance between the airfoil trailing edge and the array center;
alpha1=12/180*pi;% The angle between the line connecting midpoint of the airfoil trailing edge 
                 % and the array center, and the direction perpendicular to the airfoil plane.
alpha=20/180*pi;% Angle between the airfoil plane and the array plane;
% Array position (to be modified)
R1=0.05;
R2=0.16;
R3=0.30;
R4=0.50;

% Airfoil coordinates
ii=1:s;jj=1:t;
x1(ii)=x0+(ii-1)*l1/(s-1);
y1(jj)=(jj-(t+1)/2)*w1/(t-1);
[x2,y2]=meshgrid(x1,y1);
xa=x2';ya=y2';

% Sawtooth coordinates (don't forget to add constants in the coordinates)
ii=1:d;
xs(ii,1)=x1(s)+(ii-1)*l2/(d-1);
ys(ii,1)=0.003+(1-(t+1)/2)*w1/(t-1);
for jj=3:2:(2*b-1)
    for ii=1:d
        xs(ii,jj)=x1(s)+(ii-1)*l2/(d-1);
        ys(ii,jj)=w2/2/l2*(xs(ii,jj)-x1(s))+(jj-3)/2*w2+0.0054+(1-(t+1)/2)*w1/(t-1);
    end
end
for kk=2:2:2*b
    for ii=1:d
        xs(ii,kk)=x1(s)+(ii-1)*l2/(d-1);
        ys(ii,kk)=-w2/2/l2*(xs(ii,jj)-x1(s)-0.012)+(kk-2)/2*w2+0.003+(1-(t+1)/2)*w1/(t-1);
    end
end

% Grid coordinates
ii=1:m;jj=1:n;
x1(ii)=(ii-1)*l/(m-1);
y1(jj)=(jj-(n+1)/2)*w/(n-1);
[x2,y2]=meshgrid(x1,y1);
x=x2';y=y2';

% Microphone coordinates (microphones may not be placed directly opposite the airfoil yet, what a tragedy)
sita=2*pi/8;% Angle interval between microphones
ii=1:8;
xm(ii)=z0*sin(alpha1)+xc1+R1*cos(sita*(ii-1))*cos(alpha);ym(ii)=R1*sin(sita*(ii-1));zm(ii)=z0*cos(alpha1)-R1*cos(sita*(ii-1))*sin(alpha);% Coordinates of 1st ring microphones
ii=9:16;
xm(ii)=z0*sin(alpha1)+xc1+R2*cos(sita*(ii-9))*cos(alpha);ym(ii)=R2*sin(sita*(ii-9));zm(ii)=z0*cos(alpha1)-R2*cos(sita*(ii-9))*sin(alpha);% Coordinates of 2nd ring microphones
ii=17:24;
xm(ii)=z0*sin(alpha1)+xc1+R3*cos(sita*(ii-17))*cos(alpha);ym(ii)=R3*sin(sita*(ii-17));zm(ii)=z0*cos(alpha1)-R3*cos(sita*(ii-17))*sin(alpha);% Coordinates of 3rd ring microphones
ii=25:32;
xm(ii)=z0*sin(alpha1)+xc1+R4*cos(sita*(ii-25))*cos(alpha);ym(ii)=R4*sin(sita*(ii-25));zm(ii)=z0*cos(alpha1)-R4*cos(sita*(ii-25))*sin(alpha);% Coordinates of 4th ring microphones

fid=fopen('x.dat','wt');%Write file path
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

fid=fopen('y.dat','wt');%Write file path
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

% Compute distance vectors, distances, propagation times from grid points to microphones
for ii=1:m
    for jj=1:n
        for kk=1:p
        dx(ii,jj,kk)=xm(kk)-x(ii,jj);
        dy(ii,jj,kk)=ym(kk)-y(ii,jj);
        dz(ii,jj,kk)=zm(kk);
        rn=[dx(ii,jj,kk);dy(ii,jj,kk);dz(ii,jj,kk)];      
        rnl2(ii,jj,kk)=rn'*rn;% Compute the squared distance from the grid to the microphone
        tn(ii,jj,kk)=1/beta*(-M*rn+((M*rn)^2+beta*rnl2(ii,jj,kk))^(1/2))/c0;
        rne=[dx(ii,jj,kk);dy(ii,jj,kk);dz(ii,jj,kk)]-Urr'*tn(ii,jj,kk);
        rnel(ii,jj,kk)=c0*tn(ii,jj,kk)+M*rne;
        end
    end
end

% Compute control functions
for kk=1:p
    for ii=1:m       
        e(ii,:,kk)=exp(2*1j*pi*f*tn(ii,:,kk)).*rnel(ii,:,kk)/z0;
    end
end

% Compute control vectors
oo=1;
for jj=1:n
    for ii=1:m
        et(oo,1:p)=e(ii,jj,:);
        oo=oo+1;
    end
end

% Compute matrix A, the most important step in the algorithm
for ii=1:N0
    for jj=1:N0
        A(ii,jj)=(conj(et(ii,:))*(1./et(jj,:))'*((1./et(jj,:))*et(ii,:).'))/(p*p); 
        %A(ii,jj)=(conj(et(ii,:))*(1./et(jj,:))'*((1./et(jj,:))*et(ii,:).')-(norm((1./et(jj,:)).*et(ii,:)))^2)/(p*p-p);
    end
end

% Position of simulated sound source
X(round(((yss/w*(n-1)+(n+1)/2)-1)*m+xss/l*(m-1)+1))=0.04;% The initial value represents the square of the effective sound pressure
% X((floor(n*3/8)*m+0.5*m))=4;
% X((floor(n*6/8)*m+0.5*m))=4;

% Compute Y for each grid point
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

% Convert the energy unit of grid points to decibels
for ii=1:m
    W(ii,:)=10*log10(YM(ii,:)/4*10^10);
end

% Find the maximum position
ff = 1;
subplot(2,2,1);set(gca,'position',[t1(ff),t2(ff),0.42,0.38])
title(['\bf{Beamforming}']);hold on;
[im,jm]=find(W==max(W(:)));
imax(1)=x(im,jm);jmax(1)=y(im,jm);

% Draw the contour plot
summit=abs(max(W(:)));
v=summit-3:0.05:summit;
contour(x,y,W,v);hold on;axis equal
temp1=caxis;
caxis([low,high])
colormap(cmap);

% Draw the airfoil
plot(xa(1,:),ya(1,:));plot(xa(end,:),ya(end,:));
plot(xa(:,1),ya(:,1));plot(xa(:,end),ya(:,end));

fid=fopen('W1.dat','wt');%Write file path
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
        
        % Convert the energy unit of grid points to decibels
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
        set(gca,'position',[t1(ff),t2(ff),0.42,0.38])
        title(['\bf{DAMAS  ' num2str(ite) ' Iterations}']);hold on
        [im,jm]=find(W==max(W(:)));
        imax(ff)=x(im,jm);jmax(ff)=y(im,jm);

        % Draw contour plot
        summit = abs(max(W(:)));
        v = summit - 3:0.05:summit;
        contour(x, y, W, v); hold on; axis equal
        %caxis(temp1)
        caxis([low,high])
        colormap(cmap);

        % Draw airfoil
        plot(xa(1,:), ya(1,:)); plot(xa(end,:), ya(end,:));
        plot(xa(:,1), ya(:,1)); plot(xa(:,end), ya(:,end));
        % Output files
        tmp = strcat('W', num2str(ff), '.dat');
        fid = fopen(tmp, 'wt'); % Write file path
        matrix = W;                        
        [m, n] = size(matrix);
         for ii = 1:1:m
           for jj = 1:1:n
              if jj == n
                fprintf(fid, '%g\n', matrix(ii,jj));
              else
                fprintf(fid, '%g\t', matrix(ii,jj));
              end
           end
        end
        fclose(fid); 
        ff = ff + 1;
    end
end

toc;

h = colorbar('Location', 'southoutside');
set(h, 'Position', [0.1, 0.045, 0.8, 0.015]);

print(gcf, 'Airfoil', '-dpng', '-r300')
