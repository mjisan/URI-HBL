%%% This came from Kun
%%% Catrina modifed it a bit and added the nc stuff at the end
% To run it needs a output file from boundary_parametric and it will
% write over the variables and can then be used to run an ideal 
% boundary_model case

clc;clear;close all
%% Set up The entier domain and constant values
nx0=81; % Total grid points in x dir
ny0=81; % Total grid points in y dir
cenx=41; % TC center x index
ceny=41; % TC center y index

%hurrican touch boundary (2) or in middle of a lot of space (4)
grid = 2;

%Initial parameters
ug0 = 0; % TC translation speed
vmax = 40; % TC max wind       
rmax = 40e3; % RMW
B= 1.3; % Parameter B in Holland's Model 
%fcor=sw_f(20);
fcor=1e-4;

%IDK what this is exactly make sure its at 8 though 
rc=8*rmax; %5*rmax 10*rmax

stepx=10*1e3; % grid spacing in x dir
stepy=10*1e3; % grid spacing in y dir

%% environmental wind
%This is what is happening in the backgroud of the hurricane, it should
%just all be zero becasue rn we are only looking at a hurricane that does
%not move

%Azimuth: every value is set to the TC translation speed
ug1(1:nx0,1:ny0)=ug0;
%Radial: every value is set to zero
vg1(1:nx0,1:ny0)=0;

%Account for coriolis effect
pgfx1=-fcor*vg1;
pgfy1=fcor*ug1;

%% TC wind
%This looks at the hurricane that will be placed in the environmental wind,
%it is in a different cordinate system, built off of the radius

e = exp(1); 
f = fcor;               %coriolis effect
dr = 1e3;               %step size 
r = rmax-dr:dr:rmax*10; %Length is 362
r_in = 0:dr:rmax;       %Length is 41 %This is never used...

% Vg(r) outside of RMW 
A=rmax.^B;
c=(B/e).^0.5;
dp=(vmax./c).^2;
vg=(A*B*dp*exp(-A./r.^B)./r.^B+r.^2*f^2/4).^0.5-r*f/2; %Wind speed at each point along the radius

% Vg(r) inside of RMW - this is all at one single point, it is used in
% calculating 
a = vg(2);
b = (vg(3)+vg(1)-2*vg(2))/dr^2;
c3 = (a+b/2*rmax^2)/rmax^3;
c2 = b/2-3*rmax*c3;
c1 = -b*rmax+3*c3*rmax^2;

%%%%%%%%%%%
% For model
% So I think this is where the r of the TC is set up with in the
% environment and the way I found this code it is set up to make the r the
% entire length of the grid - 1 so the boundary around the hurrican is only
% one extra grid point

% Original code:
% dr = stepx; %Environment grid step
% r_in = dr:dr:rmax;
% r = (rmax + dr):dr:( ( (nx0-1)/2) *dr ); %r=[rmax+dr:dr:(nx0-1)/2*dr];
% ri=find(r==rc);
% nr=length(r)

% I am trying to change it so the background environmental wind extends out
% the size of the r something like below

% ---------------------------------
%                                  
% 
% 
%                 _ 
%             /        \ 
%           /            \
%         /               \
% - - - -|- - - - + - - - -|- - - -
%        \                /
%          \             /
%            \         /
%                 - 
% 
% 
% 
% ---------------------------------

%Altered Code:
dr = stepx; %Environment grid step
r_in = dr:dr:rmax;
r = (rmax + dr):dr:( ( (nx0-1)/grid) *dr );
ri=find(r==rc);
nr=length(r);
% Note - I only changed this little chunk right here and it creates what I
% wanted, I don't think it makes sense to alter the r anywhere else so it
% should be correct 

vg_in = c1*r_in+c2*r_in.^2+c3*r_in.^3;
vg = (A*B*dp*exp(-A./r.^B)./r.^B+r.^2*f^2/4).^0.5-r*f/2;
vg(ri:end) = linspace(vg(ri),0,nr-ri+1); % Here, we let vg decreases to 0 linearly 

r=[r_in,r];
va1=[vg_in,vg];
pgfa1=-(f+va1./r).*va1;

%%
fig=figure;
set(fig,'units','inches','position',[0 0 8 6]);

plot([0,r/1e3],[0,va1],'k','linewi',2);
%set(gca,'xlim',[0. 200],'xtick',[0:50:200],'ylim',[0 60],'ytick',[0:10:60])
set(gca,'fontsize',20)
ylabel('Vg (m/s)') 
xlabel('Radius (km)') 

%%
x1 = (1:nx0)*stepx;
x1 = x1 - mean(x1);
y1 = x1;

[y,x]=meshgrid(y1,x1);
dist=( (x.^2) + (y.^2) ).^0.5;
%dist is a matrix where the distance from the center is calculated at every
%point
sina=y./dist;
cosa=x./dist;

va(1:nx0,1:ny0)=0;
pgfa(1:nx0,1:ny0)=0;

% This loop seems to check each grid point and make sure that things are
% going to the right location and only things with in the radious are
% altered - so every thing in the bg env is zero
for i=1:nx0
    for j=1:ny0
        for ri=1:length(r)-1
            if dist(i,j)>=r(ri) & dist(i,j)<r(ri+1)
                va(i,j)   = va1(ri) + (va1(ri+1) - va1(ri)) / dr * (dist(i,j)-r(ri));
                pgfa(i,j) = pgfa1(ri) + (pgfa1(ri+1) - pgfa1(ri)) / dr * (dist(i,j)-r(ri));
            end
        end
    end
end

ug2=-sina.*va;
vg2=cosa.*va;
pgfx2=cosa.*pgfa;
pgfy2=sina.*pgfa;

ug2(cenx,ceny)=0;
vg2(cenx,ceny)=0;
pgfx2(cenx,ceny)=0;
pgfy2(cenx,ceny)=0;

%% nonlinear factor
dug2dx(1:nx0,1:ny0)=0;
dvg2dx(1:nx0,1:ny0)=0;

for i=2:nx0-1
    for j=2:ny0-1
        vel=ug1(i,j)+ug2(i,j);
        if vel>=0
           dug2dx(i,j)=(ug2(i,j)-ug2(i-1,j))/stepx;
           dvg2dx(i,j)=(vg2(i,j)-vg2(i-1,j))/stepx;
        else
           dug2dx(i,j)=(ug2(i+1,j)-ug2(i,j))/stepx;
           dvg2dx(i,j)=(vg2(i+1,j)-vg2(i,j))/stepx;            
        end
    end
end

%dug2dx(2:nx0-1,:)=(ug2(3:end,:)-ug2(1:end-2,:))/2/stepx;
%dvg2dx(2:nx0-1,:)=(vg2(3:end,:)-vg2(1:end-2,:))/2/stepx;

%pgfx3=ug1.*dug2dx;
%pgfy3=ug1.*dvg2dx;

%% combine the three components
ug=ug1+ug2;
vg=vg1+vg2;

%pgfx=pgfx1+pgfx2+pgfx3;
%pgfy=pgfy1+pgfy2+pgfy3;

%% derive pressure field
um=ug;
vm=vg;

pgfx(1:nx0,1:ny0)=0;
pgfy(1:nx0,1:ny0)=0;

for i=2:nx0-1
   for j=2:ny0-1

      if um(i,j)>=0 
         difu=vm(i,j)-vm(i-1,j);
      else
         difu=vm(i+1,j)-vm(i,j);
      end

      if vm(i,j)>=0 
         difv=vm(i,j)-vm(i,j-1);
      else
         difv=vm(i,j+1)-vm(i,j);
      end

      %difu=0.5*(vm(i+1,j)-vm(i-1,j));
      %difv=0.5*(vm(i,j+1)-vm(i,j-1));
      
      pgfy(i,j)=um(i,j)*difu/stepx+vm(i,j)*difv/stepy+fcor*um(i,j);
      
      if um(i,j)>=0 
         difu=um(i,j)-um(i-1,j);
      else
         difu=um(i+1,j)-um(i,j);
      end

      if vm(i,j)>=0 
         difv=um(i,j)-um(i,j-1);
      else
         difv=um(i,j+1)-um(i,j);
      end

      %difu=0.5*(um(i+1,j)-um(i-1,j));
      %difv=0.5*(um(i,j+1)-um(i,j-1));

      pgfx(i,j)=um(i,j)*difu/stepx+vm(i,j)*difv/stepy-fcor*vm(i,j);

   end
end

%% writent out
f1=fopen(['ug','.dat'],'w');
f2=fopen(['vg','.dat'],'w');
f3=fopen(['pgfx','.dat'],'w');
f4=fopen(['pgfy','.dat'],'w');

for i=1:nx0
    for j=1:ny0
        fprintf(f1,'%11.8f',ug(i,j)) ;
        fprintf(f1,'\n');

        fprintf(f2,'%11.8f',vg(i,j)) ;
        fprintf(f2,'\n');

        fprintf(f3,'%11.8f',pgfx(i,j)) ;
        fprintf(f3,'\n');

        fprintf(f4,'%11.8f',pgfy(i,j)) ;
        fprintf(f4,'\n');
    end
end
fclose(f1);
fclose(f2);
fclose(f3);
fclose(f4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check balance
tau = 5.
um(1:nx0,1:ny0)=0;
vm(1:nx0,1:ny0)=0;
wm(1:nx0,1:ny0)=0;
um1(1:nx0,1:ny0)=0;
vm1(1:nx0,1:ny0)=0;

um=ug;
vm=vg;
um1=um;
vm1=vm;

%%
for t=1:1

    t/3600

%% v eqn
for i=2:nx0-1
   for j=2:ny0-1

      if um(i,j)>=0 
         difu=vm(i,j)-vm(i-1,j);
      else
         difu=vm(i+1,j)-vm(i,j);
      end

      if vm(i,j)>=0 
         difv=vm(i,j)-vm(i,j-1);
      else
         difv=vm(i,j+1)-vm(i,j);
      end

%      difu=0.5*(vm(i+1,j)-vm(i-1,j));
%      difv=0.5*(vm(i,j+1)-vm(i,j-1));
 
      vm1(i,j)=vm(i,j)+...
      tau*(...
          -um(i,j)*difu/stepx...
          -vm(i,j)*difv/stepy...
          +pgfy(i,j)...
          -fcor*um(i,j)...
          );
   end
end

%% u eqn
for i=2:nx0-1
   for j=2:ny0-1
       
      if um(i,j)>=0 
         difu=um(i,j)-um(i-1,j);
      else
         difu=um(i+1,j)-um(i,j);
      end

      if vm(i,j)>=0 
         difv=um(i,j)-um(i,j-1);
      else
         difv=um(i,j+1)-um(i,j);
      end

%      difu=0.5*(um(i+1,j)-um(i-1,j));
%      difv=0.5*(um(i,j+1)-um(i,j-1));

      um1(i,j)=um(i,j)+...
      tau*(...
          -um(i,j)*difu/stepx...
          -vm(i,j)*difv/stepy...
          +pgfx(i,j)...
          +fcor*vm(i,j)...
          );
   end
end

um=um1;
vm=vm1;

%%
difu=um-ug;
difv=vm-vg;
cmaxu=max(max(abs(difu)));
cmaxv=max(max(abs(difv)));

figure

subplot(121)
contourf(x/1e3,y/1e3,difu);
caxis([-cmaxu cmaxu])
colorbar
axis equal
axis ([-100 100 -100 100])

subplot(122)
contourf(x/1e3,y/1e3,difv);
caxis([-cmaxv cmaxv])
colorbar
axis equal
axis ([-100 100 -100 100])

end


%%
nx2 = nx0+1;
ny2 = ny0+1;
ug_test = zeros(nx2,ny2,108);
ug_format = ug;
ug_format(nx2,:) = 0;
ug_format(:,ny2) = 0;

for t = 1:108
ug_test(:,:,t) = ug_format(:,:);
end

ug_test_sampel = ug_test(:,:,8);

vg_test = zeros(nx2,ny2,108);
vg_format = vg;
vg_format(nx2,:) = 0;
vg_format(:,ny2) = 0;

for t = 1:108
vg_test(:,:,t) = vg_format(:,:);
end

vg_test_sampel = vg_test(:,:,8);

%Add Pressure!!!!
pgfx_test = zeros(nx2,ny2,108);
pgfx_format = pgfx;
pgfx_format(nx2,:) = 0;
pgfx_format(:,ny2) = 0;

for t = 1:108
pgfx_test(:,:,t) = pgfx_format(:,:);
end
pgfx_test_sampel = pgfx_test(:,:,8);

pgfy_test = zeros(nx2,ny2,108);
pgfy_format = pgfy;
pgfy_format(nx2,:) = 0;
pgfy_format(:,ny2) = 0;

for t = 1:108
pgfy_test(:,:,t) = pgfy_format(:,:);
end

%%
cd('D:\Eyewall_Local\matlab_script\used_files');

% %File One
file_name_1 = 'boundary_parametric.nc';
% 
% %File Two
ncwrite(file_name_1, 'utop', ug_test)
% 
% %File Three
ncwrite(file_name_1, 'vtop', vg_test)
% 
% %File Four
ncwrite(file_name_1, 'pgfx', pgfx_test)
% 
% %File Five
ncwrite(file_name_1, 'pgfy', pgfy_test)

ncdisp(file_name_1)
% 
% test = ncread(file_name_1,'vtop', [1,1,5], [Inf, Inf, 1]);