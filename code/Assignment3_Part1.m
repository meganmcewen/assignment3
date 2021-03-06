%Part 1
%This code uses Monte Carlo to map moving electrons in a static electric
%field.

%constants
xrange = 2e-7; %size of area in x
yrange = 1e-7; %size of area in y
n = 1000; %number of particles
m0 = 9.10938356e-31; %electron mass
m = 0.26*m0;
T = 300; %temperature (K)
k = 1.380648e-23; %Boltzmann constant
q = 1.6e-19; %electron charge
tau = 0.2e-12;
iter = 100; %number of iterations to run the simulation
eConc = 10e-15; %electron concentration
maxt = 0;

%generate empty arrays
Px = zeros(n,1);
Py = zeros(n,1);
Vx = zeros(n,1);
Vy = zeros(n,1);
V_XScatter = zeros(n,1);
V_YScatter = zeros(n,1);

%initialize particles
Px(:,1) = xrange*rand(n,1);
Py(:,1) = yrange*rand(n,1);

%calculate electric field
voltX = 1; %voltage in x dimension
voltY = 0; %voltage in y dimension

%set up EField
EFieldX = voltX/xrange;
EFieldY = voltY/yrange;

%calculate force on electrons due to field
Fx = EFieldX*q;
Fy = EFieldY*q;

%calculate acceleration on electrons
AccX = Fx/m0;
AccY = Fy/m0;

%calculate vTH
vTH = sqrt(2*k*T/m);

%generate random velocity
randAngle = 2*pi*rand(n,1);
Vx(:,1) = vTH * cos(randAngle);
Vy(:,1) = vTH  *sin(randAngle);

%get Gaussian distribution for velocity components (scattering)
V_XScatter(:,1) = vTH.*randn(n,1);
V_YScatter(:,1) = vTH.*randn(n,1);

%check average is close to vTH
Avg = sqrt(V_XScatter.^2 + V_YScatter.^2);
VAvg = mean(Avg);

%time loop
timeStep = 1e-14;
time = 1:iter;

%set up temperature calculation
temp = zeros(iter,1);


%begin particle updating loop
for i = 1:iter
    
    %update dt
    dt = timeStep;
    
    %save old position
    Px_old = Px;
    Py_old = Py;
    
    %update velocity due to electric field
    Vx = Vx + AccX*dt;
    Vy = Vy + AccY*dt;
   
    %scattering
    Pscat = 1-exp(-dt/tau);
    ind = Pscat > rand(n,1);

    Vx(ind) = sqrt((k*T)/m).*randn(sum(ind),1);
    Vy(ind) = sqrt((k*T)/m).*randn(sum(ind),1);
   
   %update position
   Px = Px + Vx*timeStep; %updates position in x
   Py = Py + Vy*timeStep; %updates position in y
    
   %boundaries
   
   %x hitting right side
   id = Px >= xrange;
   Px(id) = Px(id) - xrange;
   Px_old(id) = Px_old(id) - xrange;

   %x hitting left side
   id = Px <= 0;
   Px(id) = Px(id) + xrange;
   Px_old(id) = Px_old(id) + xrange;
 
   %bouncing y off top/bottom
   Vy(Py >= yrange) = Vy(Py >= yrange) * -1;
   Vy(Py <= 0) = Vy(Py <= 0) * -1;
   Py(Py>yrange) = yrange-(Py(Py>yrange)-yrange);
   
   %temperature plot
   VAvg = mean(Vx.^2 + Vy.^2);
   temp(i) = (1/2)*(m*(VAvg))*(1/k);
   
   %calculate current
   mu = VAvg/EFieldX;
   I(i) = q*mu*EFieldX*eConc/(xrange*yrange);
   maxt = maxt + dt;
   
  %plot
  figure(1)
   plot([Px_old'; Px';] ,[Py_old'; Py';],"b");
   title('Particles with Scattering and Drift');
   axis ([0 200e-9 0 100e-9]);
   hold on
   drawnow
end

%plot current vs time
figure(2)
Irange = linspace(0,dt*maxt,length(I));
plot(Irange,I)
title('Current vs Time in X');
hold on


%temperature and electron density maps
x_region = linspace(0,xrange,10);
y_region = linspace(0,yrange,10);

x_bin = discretize(Px,x_region);
y_bin = discretize(Py,y_region);

temp_bin = zeros(10,10);
e_bin = zeros(10,10);

for i = 1:10
    for j  = 1:10
        inI = x_bin == i;
        inJ = y_bin == j;
        inBin = inI & inJ;
        
        sum_i = sum(Vx(inBin))/(dt/vTH);
        sum_j = sum(Vy(inBin))/(dt/vTH);
     
        sum_e = sum(inBin);
        
        avg = sqrt((sum_i)^2 + (sum_j)^2);
        
        temp_bin(i,j) = (m*0.5*avg.^2)/(k*2);
        e_bin(i,j) = sum_e;
    end
end

figure(3)
surf(e_bin);
title(['Electron Density Map after ', num2str(iter),' iterations for ' ,num2str(n),' particles']);
colorbar

figure(4)
surf(temp_bin);
title(['Temperature Map for ',num2str(n),' particles']);
colorbar

