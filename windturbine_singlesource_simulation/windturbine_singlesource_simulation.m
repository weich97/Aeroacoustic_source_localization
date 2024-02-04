%Simulation program considering mean flow.
clc;
close all
clear all
tic;

num_colors = 64;  % Number of colors
cmap = interp1([1, num_colors/2, num_colors], [0, 0, 1; 1, 1, 0; 1, 0, 0], 1:num_colors);  % Gradient from red to yellow to blue

%Allocate space for the large matrix.
m=40;%Radial grid points
n=120;%Circumferential grid points
p=32;%Number of microphones
N0=m*n;%Total number of grid points
A=zeros(N0);
dx=zeros(m,n,p);
dy=zeros(m,n,p);
dz=zeros(m,n,p);
rnl2=zeros(m,n,p);
tn=zeros(m,n,p);
rnel=zeros(m,n,p);
e=zeros(m,n,p);
et=zeros(N0,p);
X=zeros(1,N0);
Y=zeros(1,N0);
X1=zeros(1,N0);
XM=zeros(m,n);

R=38.68;%Impeller radius
R1=1; %Microphone array radius
Y0=65; %Distance from microphone array center to rotation center
c0=344;%Speed of sound 
N=0;%Speed (clockwise is positive)
U0=0;%Inflow velocity
Urr=[0,U0,0];%Velocity vector of inflow
f=1000;%Frequency
NI=5000;%Total number of iterations

%Initialize grid point coordinates
for ii=1:m
    jj=1:n;
    x(ii,:)=ii/m*R*cos((jj-1)/n*2*pi);
    z(ii,:)=ii/m*R*sin((jj-1)/n*2*pi);
end

%Initialize coordinates of the microphone circular array
ii=1:p;  
sita_2=(ii-1)/p*2*pi;
x1(ii)=R1*cos(sita_2);
y1(ii)=0;
z1(ii)=R1*sin(sita_2);

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

fid=fopen('z.dat','wt');%Write file path
matrix=z;                        
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

%Calculate distance vector, distance, propagation time from grid points to microphones
for ii=1:m
    for jj=1:n
        Ur=[N*pi/30*z(ii,jj),0,-N*pi/30*x(ii,jj)];
        M=Urr/c0;
        M1=(Urr-Ur)/c0;
        bet2(ii,jj)=1-M*M';
        for kk=1:p
        dx(ii,jj,kk)=x1(kk)-x(ii,jj);
        dy(ii,jj,kk)=-Y0;
        dz(ii,jj,kk)=z1(kk)-z(ii,jj);
        rn=[dx(ii,jj,kk);dy(ii,jj,kk);dz(ii,jj,kk)];
        rnl2(ii,jj,kk)=rn'*rn;%Calculate the square of the distance from the grid to the microphone
        tn(ii,jj,kk)=1/bet2(ii,jj)*(-M*rn+((M*rn)^2+bet2(ii,jj)*rnl2(ii,jj,kk))^(1/2))/c0;
        rne=[dx(ii,jj,kk);dy(ii,jj,kk);dz(ii,jj,kk)]-Urr'*tn(ii,jj,kk);
        rnel(ii,jj,kk)=c0*tn(ii,jj,kk)+M1*rne;
        end
    end
end
    
%Calculate control functions
rc=Y0;%Dimensionless control function
for kk=1:p
    for ii=1:m
        e(ii,:,kk)=exp(2*i*pi*f*tn(ii,:,kk)).*rnel(ii,:,kk)/rc;
    end
end

%Calculate control vectors
oo=1;
for jj=1:n
    for ii=1:m
        et(oo,1:p)=e(ii,jj,:);
        oo=oo+1;
    end
end

%Calculate matrix A, the most important step in the algorithm
for ii=1:N0
    for jj=1:N0
    A(ii,jj)=(conj(et(ii,:))*(1./et(jj,:))'*((1./et(jj,:))*et(ii,:).'))/(p*p);
    end
end

%Simulation of the source position
X((floor(n*1/8)*m+0.5*m))=4;%Initial value denotes the square of effective sound pressure
% X((floor(n*3/8)*m+0.5*m))=4;
% X((floor(n*6/8)*m+0.5*m))=4;

%Calculate Y for each grid point
for oo=1:N0
    Y(oo)=A(oo,1:N0)*X(1:N0)';
end

%Convert one-dimensional array to two-dimensional array
oo=1;
for jj=1:n
    for ii=1:m
        YM(ii,jj)=Y(oo);
        oo=oo+1;
    end
end

%Convert energy unit of grid points to decibels
for ii=1:m
    W(ii,:)=10*log10(YM(ii,:)/4*10^10);
end

%Plot
subplot(2,2,1);title(['Beamforming  ' num2str(f) 'Hz']);hold on;
[i,j]=find(W==max(W(:)));
imax(1)=round(i/m*R);jmax(1)=round(j/n*360);
plot(x(:,j),z(:,j),'k','linewidth',1.5);
plot(x(:,mod(j+n/3,n)),z(:,mod(j+n/3,n)),'k','linewidth',1.5);
plot(x(:,mod(j+2*n/3,n)),z(:,mod(j+2*n/3,n)),'k','linewidth',1.5);
summit=max(W(:));
v=summit-40:0.08:summit;
contour(x,z,W,v);
colormap(cmap);
axis equal
plot(x(m,:),z(m,:),'-');
temp1=caxis;
set(gca,'FontSize',9)
xrange_lower=-40; % sets the lower bound for the x-axis when plotting
xrange_upper=40; % sets the upper bound for the x-axis when plotting
axis([xrange_lower xrange_upper xrange_lower xrange_upper]);

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

X=Y;%Initial condition
ff=2;%Plot number initialization
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
    
    %Display every Ni steps
    if ite==100||ite==1000||ite==5000
        ite
        
        %Convert one-dimensional array to two-dimensional array
        oo=1;
        for jj=1:n
            for ii=1:m
                XM(ii,jj)=X(oo);
                oo=oo+1;
            end
        end
        
        %Convert energy unit of grid points to decibels
        for ii=1:m
            for jj=1:n
                if XM(ii,jj)==0
                    W(ii,jj)=0; 
                else
                    W(ii,jj)=10*log10(XM(ii,jj)/4*10^10);
                end
            end
        end
        
        %Plot
        subplot(2,2,ff);title(['DAMAS  ' num2str(f) 'Hz  ' num2str(ite) ' Iterations']);hold on;
        [i,j]=find(W==max(W(:)));
        imax(ff)=round(i/m*R);jmax(ff)=round(j/n*360);
        plot(x(:,j),z(:,j),'k','linewidth',1.5);
        plot(x(:,mod(j+n/3,n)),z(:,mod(j+n/3,n)),'k','linewidth',1.5);
        plot(x(:,mod(j+2*n/3,n)),z(:,mod(j+2*n/3,n)),'k','linewidth',1.5);
        summit=max(W(:));
        v=summit-40:0.08:summit;
        contour(x,z,W,v);axis equal
        colormap(cmap);
        plot(x(m,:),z(m,:),'-', 'linewidth',1.5);
        caxis(temp1)
        set(gca,'FontSize',9)
        xrange_lower=-40; % sets the lower bound for the x-axis when plotting
        xrange_upper=40; % sets the upper bound for the x-axis when plotting
        axis([xrange_lower xrange_upper xrange_lower xrange_upper]);        
        tmp=strcat('W',num2str(ff),'.dat');
        fid=fopen(tmp,'wt');%Write file path
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

% Add a shared color bar
h = colorbar('Location', 'southoutside');
set(h, 'Position', [0.1, 0.045, 0.8, 0.015]);
%h.Label.String = 'Colorbar Label';

toc;

print(gcf,'windturbine_singlesource_simulation','-dpng','-r300') % 300 dpi, i.e. 300 dots per inch
