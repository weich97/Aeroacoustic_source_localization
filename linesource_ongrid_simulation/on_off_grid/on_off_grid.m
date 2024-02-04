clc;
clear all
close all
tic;

num_colors = 64;  % Number of colors
cmap = interp1([1, num_colors/2, num_colors], [0, 0, 1; 1, 1, 0; 1, 0, 0], 1:num_colors);  % Interpolated color map from red to yellow to blue

t1=[0.06 0.55 0.06 0.55];t2=[0.55 0.55 0.05 0.05];t3=[0 100 1000 5000];
set(gcf,'position',[500 200 520 520])

% Allocate space for large matrix
% Simulate the position of sound sources
data=xlsread('xy1 data.xls')/100;
xs=data(:,1)-0.74;ys=data(:,2);
s=length(data(:,1));
ii=1:s;
Xs(ii)=4;% The initial value indicates the square of effective sound pressure
l=1.52;w=1.52;% Length and width of the grid
m=77;n=77;% Dimensions of the grid in x and y directions; must be larger than the dimensions of the airfoil
p=32; % Number of microphones
N0=m*n;% Total number of grid points
A=zeros(N0);
As=zeros(N0,s);
dx=zeros(m,n,p);
dy=zeros(m,n,p);
dz=zeros(m,n,p);
rnl2(m,n,p)=0;% Square of the distance between microphones and grid points
tn(m,n,p)=0;% Propagation time of sound waves between microphones and grid points
rnel(m,n,p)=0;% Distance between microphones and grid points
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

% Microphone coordinates (microphones may not be placed directly facing the airfoil yet, tragedy!)
sita=2*pi/p;% Angle interval between microphones
ii=1:p;
xm(ii)=R1*cos(sita*(ii-1));ym(ii)=R1*sin(sita*(ii-1));% Microphone coordinates

fid=fopen('x.dat','wt');% File path for writing
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

fid=fopen('y.dat','wt');% File path for writing
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

% Calculate the distance vector, distance, and propagation time from grid points to microphones
for ii=1:m
    for jj=1:n
        for kk=1:p
        dx(ii,jj,kk)=xm(kk)-x(ii,jj);
        dy(ii,jj,kk)=ym(kk)-y(ii,jj);
        dz(ii,jj,kk)=z0;
        rn=[dx(ii,jj,kk);dy(ii,jj,kk);dz(ii,jj,kk)];      
        rnl2(ii,jj,kk)=rn'*rn;% Calculate the square of the distance from the grid to the microphone
        tn(ii,jj,kk)=1/beta*(-M*rn+((M*rn)^2+beta*rnl2(ii,jj,kk))^(1/2))/c0;
        rne=[dx(ii,jj,kk);dy(ii,jj,kk);dz(ii,jj,kk)]-Urr'*tn(ii,jj,kk);
        rnel(ii,jj,kk)=c0*tn(ii,jj,kk)+M*rne;
        end
    end
end

% Calculate the distance vector, distance, and propagation time from sound sources to microphones
for ii=1:s
    for kk=1:p
    dxs(ii,kk)=xm(kk)-xs(ii);
    dys(ii,kk)=ym(kk)-ys(ii);
    dzs(ii,kk)=z0;
    rn=[dxs(ii,kk);dys(ii,kk);dzs(ii,kk)];      
    rnl2s(ii,kk)=rn'*rn;% Calculate the square of the distance from the grid to the microphone
    tns(ii,kk)=1/beta*(-M*rn+((M*rn)^2+beta*rnl2s(ii,kk))^(1/2))/c0;
    rne=[dxs(ii,kk);dys(ii,kk);dzs(ii,kk)]-Urr'*tns(ii,kk);
    rnels(ii,kk)=c0*tns(ii,kk)+M*rne;
    end
end

% Calculate the control function
for kk=1:p
    for ii=1:m       
        e(ii,:,kk)=exp(2*i*pi*f*tn(ii,:,kk)).*rnel(ii,:,kk)/z0;
    end
end

% Calculate the sound source control function (originally the control vector)
for kk=1:p     
    es(:,kk)=exp(2*i*pi*f*tns(:,kk)).*rnels(:,kk)/z0;
end

% Calculate the control vector
oo=1;
for jj=1:n
    for ii=1:m
        et(oo,1:p)=e(ii,jj,:);
        oo=oo+1;
    end
end

% Calculate the grid A matrix, the most important step in the algorithm
for ii=1:N0
    for jj=1:N0
        A(ii,jj)=(conj(et(ii,:))*(1./et(jj,:))'*((1./et(jj,:))*et(ii,:).'))/(p*p);        
    %A(ii,jj)=(conj(et(ii,:))*(1./et(jj,:))'*((1./et(jj,:))*et(ii,:).')-(norm((1./et(jj,:)).*et(ii,:)))^2)/(p*p-p);
    % For A, it needs to be explained that the two items in the middle are multiplied to obtain the primed matrix.
    % If A is viewed as a combination of row vectors, then each row of A should traverse the scanning region of grid points.
    % The set of each row forms the sound intensity of all grid points in the scanning region.
    end
end

% Calculate the sound source As matrix
for ii=1:N0
    for jj=1:s
        As(ii,jj)=(conj(et(ii,:))*(1./es(jj,:))'*((1./es(jj,:))*et(ii,:).'))/(p*p);        
    end
end

% Calculate the Y of each grid point
for oo=1:N0
    Y(oo)=As(oo,1:s)*Xs(1:s)';
end

% Convert a one-dimensional array to a two-dimensional array
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
subplot(2,2,1);
set(gca,'ytick', [-0.5 0 0.5])
title(['\bf{Beamforming}']);
hold on;
[im,jm]=find(W==max(W(:)));
imax(1)=x(im,jm);jmax(1)=y(im,jm);

% Draw contour plot
summit=abs(max(W(:)));
v=summit-40:0.08:summit;
contour(x,y,W,v);hold on;axis equal;axis([-0.8,0.8,-0.8,0.8]);
colormap(cmap);

fid=fopen('W1.dat','wt');% File path for writing
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

X=Y;% Initial condition
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
        
        % Convert a one-dimensional array to a two-dimensional array
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
        set(gca,'ytick', [-0.5 0 0.5])
        title(['\bf{DAMAS  ' num2str(ite) ' Iterations}']);
        hold on;
        [im,jm]=find(W==max(W(:)));
        imax(ff)=x(im,jm);jmax(ff)=y(im,jm);
        % Draw contour plot
        summit=abs(max(W(:)));
        v=summit-40:0.08:summit;
        contour(x,y,W,v);hold on;axis equal;axis([-0.8,0.8,-0.8,0.8]);
        colormap(cmap);
        
        tmp=strcat('W',num2str(ff),'.dat');
        fid=fopen(tmp,'wt');% File path for writing
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

print(gcf,'UCAS_10000_nongrid','-dpng','-r300')
