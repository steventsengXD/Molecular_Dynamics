% MD of fluids by S Tseng

% A molecular dynamics (MD) program to simulate the motion of a collection
% of identical particles in two dimensions. 

% In regards to conservation laws, without any outside forces, we would 
% expect for this system to have its total energy (K+P) to be conserved 
% at each time step. With the use of periodic boundary conditions, 
% the number of particles, N, is conserved as well. 

clear all
close all

% Parameter Definitions
N=100; % Number of particles
sigma=1.5; % Radius where phi(r_ij)=0 
m=1; % mass of each particle
eps=500; % Interaction strength
L=30; % Length of a wall of square box,particles are in [-L/2,L/2]x[-L/2,L/2]

% Time, Steps
steps=10^3; % Number of steps 
tmax=10; 
t0=0; 
dt=(tmax-t0)/steps;

% Initial positions of particles, r_i=(x_i,y_i)
r=L*(rand(N,2)-0.5); 
lim=[-1.5*(L/2),1.5*(L/2)]; % Plot range
figure;
plot(r(:,1),r(:,2),'.','markersize',40)
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
   
% v_0 initialized randomly, v_init=v_{-1/2}
v_0=randn(N,2);
v_init=v_0-(dt/2)*force(r,N,i,j,sigma,eps)*(1/m); 
v=v_init;

for step=1:steps    
    
    % Leap-frog Verlet Algorithm
    v=v+force(r,N,i,j,sigma,eps)*(1/m)*dt;
    r=r+dt*v;

    % Periodic Boundary Conditions 
    % When particle leaves LxL box, particle's image enters on opposite side
    r=mod(r+(L/2),L)-(L/2); 
    
    % Plot
    plot(r(:,1),r(:,2),'.','markersize',40)
    hold on
    line([-1.5*L,1.5*L],[L/2,L/2])
    line([-1.5*L,1.5*L],[-L/2,-L/2])
    line([L/2,L/2],[-1.5*L,1.5*L])
    line([-L/2,-L/2],[-1.5*L,1.5*L])
    xlim(lim)
    ylim(lim)
    grid on
    gf=getframe;
    hold off
      
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


function f_i=force(r,N,i,j,sigma,eps)

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





