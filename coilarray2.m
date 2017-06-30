%-------------------------------------------------------------------------%
%-------Cicular Coil Magnetic Field Simulation--------------------------%
%-------------------------------------------------------------------------%

clc
close all
clear all

%-------------------------------------------------------------------------%
%---Coil is in the X-Y plane and Magnetic Field is Evaluated -------------%
%-------------at every point in the Y-Z plane(X=0)------------------------%
%-------------------------------------------------------------------------%
coils = 1;
Nz=51;  % No. of grids in Z-axis
Ny=51;  % No. of grids in Y-axis
N=25;   % No of grids in the coil ( X-Y plane)
Ra=3;    % Radius of the coil in the X-Y plane
I=6;    % current in the coil
u0=1;   % for simplicity, u0 is taken as 1 (permeability)
phi=-pi/2:2*pi/(N-1):3*pi/2; % For describing a circle (coil)
xoffset = 0;
yoffset = -3;%-3

zoffset = 0;

yp(1:51)=-25:1:25; % Y-coordinates of the plane 
zp(1:51)=-25:1:25;% Z-coordinates of the plane 

% Y(1:Ny,1:Nz)=0; % This array is for 1-d to 2-d conversion of coordinates
% Z(1:Ny,1:Nz)=0; 
% 
% for i=1:Ny
%     Y(i,:)=yp(i); % y-coordinates 2-d form (entire row(i) = yp(i))
% end
% for i=1:Nz 
%     Z(:,i)=zp(i); % z-coordinates 2-d form (enitre col(i) = zp(i))
% end

BXtot = zeros(Ny, Nz);
BYtot = zeros(Ny, Nz);
BZtot = zeros(Ny, Nz);
BRtot = zeros(Ny, Nz);
%---------------------------------------------------- --------------------
%-variable "a" in for loop is along Y and variable "b" is along Z-axis----
%-------------------------------------------------------------------------
for w=1:coils
Xc=Ra*cos(phi) + xoffset; % X-coordinates of the coil
Yc=Ra*sin(phi) + yoffset; % Y-coordinates of the coil
null=zeros(N);
Zc=Ra*cos(phi) + xoffset; % Z-coordinates of the coil

figure(1)
plot3(Xc,Yc,null,'linewidth',3)
hold on

yoffset = yoffset + 6;

Xc=Ra*cos(phi) + xoffset; % X-coordinates of the coil
Yc=Ra*sin(phi) + yoffset; % Y-coordinates of the coil
null=zeros(N);
Zc=Ra*cos(phi) + xoffset; % Z-coordinates of the coil

plot3(Xc,Yc,null,'linewidth',3)
hold on

yoffset = 0;
zoffset = -3;

Xc=Ra*cos(phi) + xoffset; % X-coordinates of the coil
Yc=Ra*sin(phi) + zoffset; % Y-coordinates of the coil
null=zeros(N);
Zc=Ra*cos(phi) + zoffset; % Z-coordinates of the coil


plot3(Xc,null,Yc,'linewidth',3)

zoffset = zoffset + 6;


Xc=Ra*cos(phi) + xoffset; % X-coordinates of the coil
Yc=Ra*sin(phi) + zoffset; % Y-coordinates of the coil
null=zeros(N);
Zc=Ra*cos(phi) + zoffset; % Z-coordinates of the coil

plot3(Xc,null,Yc,'linewidth',3)

%I = -I;

%axis([-10 10 -10 10])
xlabel('X-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
zlabel('Z-axis','fontsize',14)
title('Loop Co-ordinates','fontsize',14)
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
h = get(gca, 'ylabel');
fh = figure(1); 
set(fh, 'color', 'white'); 
grid on
hold on

for a=1:Ny  
for b=1:Nz 


    
%-------------------------------------------------------------------------    
% calculate R-vector from the coil(X-Y plane)to Y-Z plane where we are 
% interested to find the magnetic field and also the dl-vector along the
% coil where current is flowing
% (note that R is the position vector pointing from coil (X-Y plane) to
% the point where we need the magnetic field (in this case Y-Z plane).)
% dl is the current element vector which will make up the coil------------
%-------------------------------------------------------------------------
%up is +ve vector down is -ve vector 
for i=1:N-1
Rx(i)=-0.5*(Xc(i)+Xc(i+1));
Ry(i)=(yp(a)-(0.5*(Yc(i)+Yc(i+1))));
Rz(i)=zp(b);
dlx(i)=Xc(i+1)-Xc(i);
dly(i)=Yc(i+1)-Yc(i);
end
Rx(N)=-0.5*(Xc(N)+Xc(1));
Ry(N)=(yp(a)-(0.5*(Yc(N)+Yc(1))));
Rz(N)=zp(b);
dlx(N)=-Xc(N)+Xc(1);
dly(N)=-Yc(N)+Yc(1);

%--------------------------------------------------------------------------
% for all the elements along coil, calculate dl cross R -------------------
% dl cross R is the curl of vector dl and R--------------------------------
% XCross is X-component of the curl of dl and R (dl x R), similarly I get Y and Z- 
%--------------------------------------------------------------------------
%(yz - zy)i - (xz - zx)j + (xy - yx)k
%(dly*Rz)i - (dlx*Rz)j + (dlx*Ry - dly*Rx)k 
for i=1:N
Xcross(i)=dly(i).*Rz(i);
Ycross(i)=-dlx(i).*Rz(i);
Zcross(i)=(dlx(i).*Ry(i))-(dly(i).*Rx(i));
R(i)=sqrt(Rx(i).^2+Ry(i).^2+Rz(i).^2);
end

%-------------------------------------------------------------------------
% this will be the biot savarts law equation------------------------------
%--------------------------------------------------------------------------

Bx1=(I*u0./(4*pi*(R.^3))).*Xcross;
By1=(I*u0./(4*pi*(R.^3))).*Ycross;
Bz1=(I*u0./(4*pi*(R.^3))).*Zcross;
%--------------------------------------------------------------------------
% now we have  magnetic field from all current elements in the form of an
% array named Bx1,By1,Bz1, now its time to add them up to get total
% magnetic field 
%-------------------------------------------------------------------------
BX(a,b)=0;       % Initialize sum magnetic field to be zero first
BY(a,b)=0;
BZ(a,b)=0;


%--------------------------------------------------------------------------
% here we add all magnetic field from different current elements which 
% make up the coil
%--------------------------------------------------------------------------

for i=1:N   % loop over all current elements along coil    
    BX(a,b)=BX(a,b)+Bx1(i);
    BY(a,b)=BY(a,b)+By1(i);
    BZ(a,b)=BZ(a,b)+Bz1(i);
end

%-------------------------------------------------------------------------

    BXtot(a,b) = BXtot(a,b) + BX(a,b);
    BYtot(a,b) = BYtot(a,b) + BY(a,b);
    BZtot(a,b) = BZtot(a,b) + BZ(a,b);

    BRtot(a,b) = sqrt(BYtot(a,b).^2 + BZtot(a,b).^2);
    %sqauring loses negative sign 
end
end
end

BYtot = BYtot + BZtot;
BXtot = BXtot + BXtot;
BZtot = BZtot + BYtot;



%------------------------------------------------------------------------
%---BX is a null field and only BY and BZ are color mapped...............
%--------------------------------------------------------------------------



figure(2)
lim1=min(min(BRtot));
lim2=max(max(BRtot));
steps=(lim2-lim1)/100;
%contour(zp,yp,BRtot,lim1:steps:lim2) %plots contours with strength 
q = quiver(zp,yp,BZtot,BYtot,3)
axis([-10 10 -10 10])
xlabel('Z-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
title('Single loop','fontsize',14)
colorbar('location','eastoutside','fontsize',14);
h=gca; %get current axes 
get(h,'FontSize') 
set(h,'FontSize',14)
h = get(gca, 'ylabel');
fh = figure(2); 
set(fh, 'color', 'white'); 
% 
% 
% figure(3)
% lim1=min(min(BYtot));
% lim2=max(max(BYtot));
% steps=(lim2-lim1)/100;
% contour(zp,yp,BYtot,lim1:steps:lim2)
% axis([-25 25 -25 25])
% xlabel('Z-axis','fontsize',14)
% ylabel('Y-axis','fontsize',14)
% title('BY component','fontsize',14)
% colorbar('location','eastoutside','fontsize',14);
% h=gca; 
% get(h,'FontSize') 
% set(h,'FontSize',14)
% h = get(gca, 'ylabel');
% fh = figure(3); 
% set(fh, 'color', 'white');
% 
% figure(4)
% lim1=min(min(BZtot));
% lim2=max(max(BZtot));
% steps=(lim2-lim1)/100;
% contour(zp,yp,BZtot,lim1:steps:lim2) %plots contours with strength 
% axis([-25 25 -25 25])
% xlabel('Z-axis','fontsize',14)
% ylabel('Y-axis','fontsize',14)
% title('BZ component','fontsize',14)
% colorbar('location','eastoutside','fontsize',14);
% h=gca; %get current axes 
% get(h,'FontSize') 
% set(h,'FontSize',14)
% h = get(gca, 'ylabel');
% fh = figure(2); 
% set(fh, 'color', 'white'); 
hold on
%figure(5);

c = q.Color;
q.Color = 'red';
%axis([-25 25 -25 25]);
%xlabel('Z-axis','fontsize',14);
% ylabel('Y-axis','fontsize',14);
% title('B-field Vector flow','fontsize',14);
% h=gca; 
% get(h,'FontSize');
% set(h,'FontSize',14);
% h = get(gca, 'ylabel');
% fh = figure(4); 
% set(fh, 'color', 'white'); 

%--------------------------------------------------------------------------
