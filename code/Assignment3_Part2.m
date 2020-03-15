clear all

%this code uses the finite difference method to generate an electric field
%solution. It is used for Parts 2a and 2b only.

%constants
L = 200; %size of matrix in x
W = 100; %size of matrix in y
Lb = 40; %length of contact
Wb = 40; %width of contact
sigmaInside = 1e-2; %sigma value for contacts
sigmaOutside = 1;
q = 1.6e-19; %electron charge
m0 = 9.10938356e-31; %electron mass
T = 300; %temperature (K)
k = 1.380648e-23; %Boltzmann constant
m = 0.26*m0;
eConc = 10e-15; %electron concentration
iter = 50; %number of iterations to run the simulation
tau = 0.2e-12;

%create initial matrices
G = sparse(L*W,L*W);
F = zeros(1,L*W);
sigma = ones(W,L);

%create contacts
contactLeft = (L/2) - (Lb/2);
contactRight = (L/2) + (Lb/2);
contactBottom = Wb;
contactTop = W - Wb;

%set up sigma for contacts
sigma(:,:) = sigmaOutside;
sigma(contactTop:W, contactLeft:contactRight) = sigmaInside;
sigma(1:contactBottom, contactLeft:contactRight) = sigmaInside;

%G matrix
for x=1:L
    for y=1:W
        
        %mapping equation
        n = y +(x-1)*W;
        
        if(x==1) 
         
            G(n,n) = 1;
            F(n)= 0.8;
            
        elseif(x==L)
            
            G(n,n) = 1;
            F(n)=0;
            
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

VMatrix = zeros(W,L);
%populate V matrix solution
for (x = 1:L)
    for (y = 1:W)
        n = y + (x-1)*W;
        VMatrix(y,x) = V(n);
    end
end

%solve for E 
[Ex_FDM1, Ey_FDM1] = gradient(-VMatrix*1e9);


%plots
figure(1)
surf(VMatrix);
title('Voltage Map');
colormap

figure(2)
surf(sigma);
title('Contacts inside Box');
view(2)

figure(3)
quiver(Ex_FDM1,Ey_FDM1)
title('Electric Field Plot');
view(0,90);


