% MD of fluids by S Tseng

% A molecular dynamics (MD) program to simulate the motion of a collection
% of identical particles in two dimensions. 

% In regards to conservation laws, without any outside forces, we would 
% expect for this system to have its total energy (K+P) to be conserved 
% at each time step. With the use of periodic boundary conditions, 
% the number of particles, N, is conserved as well. 

% Using a soft potential function like in (1) does not seem to keep 
% all particles within the boundary though. 

clear all
close all

% Parameter Definitions
N=100; % Number of particles
sigma=9.5; % Radius where phi(r_ij)=0 
m=7; % mass of each particle
eps=150; % Interaction strength
epsw=10.0*eps; 
L=20; % Length of a wall of square box,particles are in [-L/2,L/2]x[-L/2,L/2]

% Time, Steps
steps=10^3; % Number of steps 
tmax=10; 
t0=0; 
dt=(tmax-t0)/steps;

% Initial positions of particles, r_i=(x_i,y_i)
r=L*(rand(N,2)-0.5); 
lim=[-1.5*(L/2),1.5*(L/2)]; % Plot range
figure;
plot(r(:,1),r(:,2),'o','markersize',40)
hold on
line([-1.5*L,1.5*L],[L/2,L/2])
line([-1.5*L,1.5*L],[-L/2,-L/2])
line([L/2,L/2],[-1.5*L,1.5*L])
line([-L/2,-L/2],[-1.5*L,1.5*L])
xlim(lim)
ylim(lim)
grid on
hold off

[i,j]=meshgrid(1:N,1:N); % Initialize i,j values 

% Total force on particle i = f_i + f_wall
Ftotal=force_p(r,N,i,j,sigma,eps)+force_w(r,N,L,epsw);    
   
% v_0 initialized randomly, v_init=v_{-1/2}
v_0=randn(N,2);
v_init=v_0-(dt/2)*Ftotal*(1/m); 
v=v_init;

for step=1:steps    
    
    % Leap-frog Verlet Algorithm
    v=v+Ftotal*(1/m)*dt;
    r=r+dt*v;
    
    % Plot
    plot(r(:,1),r(:,2),'.','markersize',65)
    xlim(lim)
    ylim(lim)
    hold on
    line([-1.5*L,1.5*L],[L/2,L/2])
    line([-1.5*L,1.5*L],[-L/2,-L/2])
    line([L/2,L/2],[-1.5*L,1.5*L])
    line([-L/2,-L/2],[-1.5*L,1.5*L])
    xlim(lim)
    ylim(lim)
    grid on
    hold off
    gf=getframe;
    
    % Total force on particle i = f_i + f_wall
    Ftotal=force_p(r,N,i,j,sigma,eps)+force_w(r,N,L,epsw);   
    
    if (mod(step,200) ~= 0)
        continue 
    else
        disp('steps = '); disp(step);
    end
    
end

% figure;
% scatter(r(:,1),r(:,2),100,'r','filled')
% xlim(lim)
% ylim(lim)
% grid on

function f_i=force_p(r,N,i,j,sigma,eps)

% Differences between all particle coordinates r_ij
dr_ij(:,:,1)=reshape(r(i,1)-r(j,1),N,N); % Array for (x_i-x_j) for all i,j
dr_ij(:,:,2)=reshape(r(i,2)-r(j,2),N,N); % Array for (y_i-y_j) for all i,j 

% Distance between all particles, but ignoring diagonal
norm_r=sqrt(dr_ij(:,:,1).^2+dr_ij(:,:,2).^2)+eye(N);

% Partial Derivatives
Dx=dr_ij(:,:,1)./norm_r(:,:); % Partial derivative of r_ij wrt x
Dy=dr_ij(:,:,2)./norm_r(:,:); % Partial derivative of r_ij wrt y
Dphi=-2*eps*(sigma-norm_r); % Partial phi before truncation
Dphi(logical(eye(size(Dphi))))=0; % Zeroing diagonal of partial phi 
%I=eye(N); Dphi=Dphi.*(~I); % Can also do this to zero out phi diagonal

% Based on analytical solution for partial of phi, any element in Dphi>0
% before truncation should be set to zero to account for the truncation
% distance
Dphi(Dphi>0)=0;

% f is the force between particles, f_ij
f(:,:,1)=-1*Dphi.*Dx;
f(:,:,2)=-1*Dphi.*Dy;

% f_i = Particle forces acting on particle i
f_i(:,1)=sum(f(:,:,1),2);
f_i(:,2)=sum(f(:,:,2),2);
end


function fw=force_w(r,N,L,epsw)
% Boundary condition physical walls
% Difference between each wall and particle
dr_w(:,1)=r(:,1)-(L/2)*ones(N,1);
dr_w(:,2)=r(:,2)-(L/2)*ones(N,1);
dr_w(:,3)=r(:,1)-(-1)*(L/2)*ones(N,1);
dr_w(:,4)=r(:,2)-(-1)*(L/2)*ones(N,1);
normr_w=abs(dr_w); % Distance between each particle and wall
Dxyw=dr_w./normr_w; % Partial derivatives wrt x,y
Dphiw=-2*epsw*((L/10)-normr_w);
Dphiw(Dphiw>0)=0;
% Getting wall forces 
fw1=-1*Dphiw.*Dxyw;
fw(:,1)=fw1(:,1)+fw1(:,3); % sum of wall forces in x direction
fw(:,2)=fw1(:,2)+fw1(:,4); % sum of wall forces in y direction
end



