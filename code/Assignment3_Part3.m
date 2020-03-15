%this code combines both parts, using the finite difference method to
%generate an electric field which is then used as an input to the Monte
%Carlo setup. This code is used for Part 2c and all of Part 3.

%set up Finite Difference
%constants
L = 200; %size of matrix in x
W = 100; %size of matrix in y
Lb = 40; %length of contact
Wb = 40; %width of contact
sigmaInside = 1e-2; %sigma value for contacts
sigmaOutside = 1;
V0 = 0.8;

%create initial matrices
G = sparse(L*W,L*W);
F = zeros(1,L*W);
sigma = ones(W,L);

%create contacts
contactLeft = 80;
contactRight = 120;
contactBottom = 40;
contactTop = 60;

%set up sigma for contacts
sigma(:,:) = sigmaOutside;
sigma(contactTop:W, contactLeft:contactRight) = sigmaInside;
sigma(1:contactBottom, contactLeft:contactRight) = sigmaInside;

figure(2)
surf(sigma)

%G matrix
for x=1:L
    for y=1:W
        
        %mapping equation
        n = y +(x-1)*W;
        
        if(x==1) 
         
            G(n,n) = 1;
            F(n) = 0.8;
            
        elseif(x==L)
            
            G(n,n) = 1;
            F(n) = 0;
            
        elseif(y == 1)
            
            %local mapping
            nyp = y+1+(x-1)*W;
            nxp = y+(x)*W;
            nxm = y+(x-2)*W;

            sig_yp =(sigma(y,x)+sigma(y+1,x))/2;
            sig_xp=(sigma(y,x)+ sigma(y,x+1))/2;
            sig_xm =(sigma(y,x)+ sigma(y,x-1))/2;
            
            G(n,n)= -(sig_yp+sig_xp+sig_xm);
            G(n,nyp)= sig_yp;
            G(n,nxp)=sig_xp;
            G(n,nxm)= sig_xm;
            
        elseif(y==W)
            
            %local mapping
            nxp = y+(x)*W;
            nxm = y+(x-2)*W;
            nym = y-1+(x-1)*W;

            sig_xp=(sigma(y,x)+ sigma(y,x+1))/2;
            sig_xm =(sigma(y,x)+ sigma(y,x-1))/2;
            sig_ym =(sigma(y,x)+ sigma(y-1,x))/2;
           
            G(n,n)=-(sig_ym+sig_xp+sig_xm);
            G(n,nym)=sig_ym;
            G(n,nxp)=sig_xp;
            G(n,nxm)=sig_xm;
            
        else
            
            %local mapping
            nyp = y+1+(x-1)*W;
            nxp = y+(x)*W;
            nxm = y+(x-2)*W;
            nym = y-1+(x-1)*W;

            sig_yp =(sigma(y,x)+sigma(y+1,x))/2;
            sig_xp=(sigma(y,x)+ sigma(y,x+1))/2;
            sig_xm =(sigma(y,x)+ sigma(y,x-1))/2;
            sig_ym =(sigma(y,x)+ sigma(y-1,x))/2;
        
            G(n,n)=-(sig_yp+sig_ym+sig_xp+sig_xm);
            G(n,nyp)= sig_yp;
            G(n,nym)= sig_ym;
            G(n,nxp)= sig_xp;
            G(n,nxm)= sig_xm;
            
        end
    end
end

V = G\F';

for x=1:L
    for y=1:W
        n = y+(x-1)*W;
        VMatrix(y,x) = V(n);
    end
end

[Ex_Map, Ey_Map] = gradient(-VMatrix*1e9);

figure(1)
quiver(Ex_Map, Ey_Map);
title('Electric Field Map');

figure(2)
mesh(VMatrix)
colormap
title('Voltage Map');

%now do Monte Carlo
%constants
xrange = 2e-7; %size of area in x
yrange = 1e-7; %size of area in y
m0 = 9.10938356e-31; %electron mass
m = 0.26*m0;
T = 300; %temperature (K)
k = 1.380648e-23; %Boltzmann constant
q = 1.6e-19; %electron charge
tau = 0.2e-12;
iter = 100; %number of iterations to run the simulation
eConc = 10e-15; %electron concentration
maxt = 0;
n = 50;

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

%calculate force on electrons due to field
Fx = Ex_Map*q;
Fy = Ey_Map*q;

%calculate acceleration on electrons
AccX = Fx/m;
AccY = Fy/m;

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

%define boxes for bottleneck
boxLeft = 0.8e-7;
boxRight = 1.2e-7;
boxBottom = 0.6e-7;
boxTop = 0.4e-7;

%define inside of boxes
inX = (boxLeft < Px) & (Px < boxRight);
inY = (boxBottom < Py) | (boxTop > Py);
inBox = inX&inY;

%check particles are not initially in boxes and move if needed
numInBox = 1;
while(numInBox > 0)
   
    Px(inBox) = xrange*rand(1,numInBox); 
  
    inX = (boxLeft < Px) & (Px < boxRight);
    inY = boxBottom < Py | boxTop > Py;
    inBox = inX&inY;
    
    numInBox = sum(inBox);
end

%begin particle updating loop
for i = 1:iter
    
    %update dt
    dt = timeStep;
    
    %save old position
    Px_old = Px;
    Py_old = Py;
    
    %update velocity due to electric field
    px_r = Px.*1e9;
    py_r = Py.*1e9;
    
    %get rounded position numbers to index the acceleration matrix
    for i = 1:n
    if (round(py_r(i)) <= 0)
        py_r(i) = 1;
    end
    if(round(px_r(i)) <= 0)
        px_r(i) = 1;
    end
    if(round(py_r(i)) > W)
        py_r(i) = W - 1;
    end
    if(round(px_r(i)) > L)
        px_r(i) = L - 1;
    end
    end
    
    Acc_Px = zeros(n,1);
    Acc_Py = zeros(n,1);
    
    %get new acceleration vector sized to the number of particles by grabbing the acceleration for each
    %particle in the x and y direction
    for i =  1:n
        Acc_Px(i) = AccX(round(py_r(i)),round(px_r(i)));
        Acc_Py(i) = AccY(round(py_r(i)),round(px_r(i)));
    end
    
    %calulate velocity due to electric field
    Vx = Vx + Acc_Px.*dt;
    Vy = Vy + Acc_Py.*dt;
   
    %scattering
    Pscat = 1-exp(-dt/tau);
    ind = Pscat > rand(n,1);

    Vx(ind) = sqrt((k*T)/m).*randn(sum(ind),1);
    Vy(ind) = sqrt((k*T)/m).*randn(sum(ind),1);
   
   %update position
   Px = Px + Vx*timeStep; %updates position in x
   Py = Py + Vy*timeStep; %updates position in y
    
   %boundaries
   %bouncing off boxes if box boundary encountered
   inX = (boxLeft < Px) & (Px < boxRight);
   inY = boxBottom < Py | boxTop > Py;
   inBox = inX&inY;
    
   betweenBoxes = (Py_old > boxTop)&(Py_old < boxBottom);
   Vy(inBox&betweenBoxes) = -1*Vy(inBox&betweenBoxes);
   Vx(inBox&~betweenBoxes) = -1*Vx(inBox&~betweenBoxes);
   
   for (i = 1:n)
       if (boxLeft < Px(i)) & (Px(i) < boxRight) & (boxBottom < Py(i) | boxTop > Py(i))
           Px(i) = Px_old(i);
           Py(i) = Py_old(i);
       end
   end
   
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
   E_avg = sqrt(Ex_Map.^2 + Ey_Map.^2);
   mu = VAvg/E_avg(i);
   I(i) = q*mu*E_avg(i)*eConc/(xrange*yrange);
   maxt = maxt + dt;
   
  %plot
  figure(3)
   plot([Px_old'; Px';] ,[Py_old'; Py';],"b");
   title('Particles with Scattering and Drift');
   axis ([0 200e-9 0 100e-9]);
   hold on
   drawnow
   %add boxes to particle plot
   rectangle('Position',[0.8e-7 0 0.4e-7 0.4e-7])
   rectangle('Position',[0.8e-7 0.6e-7 0.4e-7 0.4e-7])

end

%plot current vs time
figure(4)
Irange = linspace(0,dt*maxt,length(I));
plot(Irange,I)
title('Current vs Time in X');
hold on

%temperature map
x_region = linspace(0,xrange,100);
y_region = linspace(0,yrange,100);

x_bin = discretize(Px,x_region);
y_bin = discretize(Py,y_region);

e_bin = zeros(100,100);

for i = 1:100
    for j  = 1:100
        inI = x_bin == i;
        inJ = y_bin == j;
        inBin = inI & inJ;
        
        sum_i = sum(Vx(inBin))/(dt/vTH);
        sum_j = sum(Vy(inBin))/(dt/vTH);
     
        sum_e = sum(inBin);
        
        e_bin(i,j) = sum_e;
    end
end

figure(5)
surf(e_bin);
title(['Electron Density Map after ', num2str(iter),' iterations for ' ,num2str(n),' particles']);
colorbar

figure(6)



