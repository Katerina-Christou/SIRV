close all;
clear all;

%% --- Parameter Initialization ---
% Diffusion parameters for each sub-population
delta_u = 0.01; %Infected diffusion parameter
delta_z = 0.01; %Susceptible diffusion parameter
delta_r = 0.01; %Recovered diffusion parameter
delta_v = 0.01; %Vaccinated diffusion parameter

% Epidemiological parameters
beta = 0.75; %Risk of infection
lambda = 0.20; %Recovery rate of infectious individuals
alpha = 0.035; %Virus-induced average fatality rate
epsilon = 0.8; %Vaccine inefficacy
rho = 0.06; %Vaccination rate

% Advection/Pressure related parameters
mu = 3; %Interface condition parameter
Ksi = 0.02; %Cross-diffusion pressure due to self-isolation
Ksd = 0.002; %Self-diffusion pressure due to social distancing

%% --- Spatial and Temporal Discretisation ---
N = 100; %Number of Spatial nodes
dt = 0.0001; %Time-step value
T = 37.5; %Final time
Ntimesteps = T / dt; %Number of timesteps
dx = 1 / N; %Space-step value
x = linspace(0, 1, N + 1);
h = 31; %Interface position

%% --- Initial Conditions ---
%Infected Population density
u(1:h)=0.53*((x(1:h)).*(0.3-x(1:h)));
u(1:16)=0.012;
u(h:101)=0;
%Susceptible Population density
z=ones(1,N+1);
%Recovered Population density
r=zeros(1,N+1);
%Vaccinated Population density
v=zeros(1,N+1);
%Combined Population density
total=u+z+r+v;

 
%% --- Initial Plotting ---
figure('Name','Initial Conditions');
plot(x, total,'b','DisplayName','Total Population');
hold on;
plot(x,u,'g','DisplayName','Infected (u)');
plot(x,z,'r','DisplayName','Susceptible (z)');
plot(x,r,'m','DisplayName','Recovered (r)');
plot(x,v,'c','DisplayName','Vaccinated (v)');
xlabel('x');
ylabel('Sub-population Densities');
title('Initial Population Densities');
legend('Location','best');
xlim([0 1]);
grid on;
hold off;

%% --- Pre-computation for Mass and Theta (based on initial conditions) ---
%Theta Calculation LHS of Interface
theta1= trapz(x(1:h),total(1:h));
%Mass Calculation LHS of Interface
mass1(1:h-1)=(1/theta1)*(total(1:h-1).*(x(2:h)-x(1:h-1)));
mass1(h-1)=(1/theta1)*(total(h-1).*(x(h-1)-x(h-2))); 
mass1(h)=(1/theta1)*(total(h).*(x(h)-x(h-1)));
%Theta Calculation RHS of Interface
theta2= trapz(x(h:N+1),total(h:N+1));
%Mass Calculation RHS of Interface
mass2(h)=(1/theta2)*(total(h).*(x(h+1)-x(h)));
mass2(h+1:N)=(1/theta2)*(total(h+1:N).*(x(h+2:N+1)-x(h+1:N)));
mass2(N)=(1/theta2)*(total(N).*(x(N)-x(N-1)));        
mass2(N+1)=(1/theta2)*(total(N+1).*(x(N+1)-x(N)));
%% --- Time Stepping Loop ---
for t=1:Ntimesteps
%Interface Velocity Calculation    
w(h)=-(delta_u*mu)*((u(h)-u(h-1))/(x(h)-x(h-1)));

%dtheta Calculation LHS of Interface
dtheta1 = (delta_u*(u(h)-u(h-1))+delta_z*(z(h)-z(h-1))+delta_r*(r(h)-r(h-1))+delta_v*(v(h)-v(h-1)))/(x(h)-x(h-1))+...
           (Ksi *(z(h)+r(h)+v(h)+u(h)))*(u(h)-u(h-1))/(x(h)-x(h-1))+(Ksd*(z(h)+r(h)+v(h))+Ksi*u(h))*((z(h)-z(h-1))+(r(h)-r(h-1))+(v(h)-v(h-1)))/(x(h)-x(h-1))+...
           trapz(x(1:h),(-alpha*u(1:h)))+total(h)*w(h);
%dtheta Calculation RHS of Interface
dtheta2 = -(delta_z*(z(h+1)-z(h))+delta_r*(r(h+1)-r(h))+delta_v*(v(h+1)-v(h)))/(x(h+1)-x(h))-...
            Ksi*z(h)*(u(h+1)-u(h))/(x(h+1)-x(h))-Ksd*z(h)*((z(h+1)-z(h))+(r(h+1)-r(h))+(v(h+1)-v(h)))/(x(h+1)-x(h))-total(h)*w(h);

%Velocity is zero at the LHS Boundary
w(1)=0;
%Velocity Calculation LHS of Interface (LHS boundary as the anchor point)
for k=2:h-1 
w(k)=(1/total(k))*((-(delta_u*(u(k)-u(k-1))+delta_z*(z(k)-z(k-1))+delta_r*(r(k)-r(k-1))+delta_v*(v(k)-v(k-1))))/(x(k)-x(k-1))-...
    Ksi*(z(k)+r(k)+v(k)+u(k))*(u(k)-u(k-1))/(x(k)-x(k-1))-Ksd*(z(k)+v(k)+r(k))*((z(k)-z(k-1))+(r(k)-r(k-1))+(v(k)-v(k-1)))/(x(k)-x(k-1))-Ksi*u(k)*((z(k)-z(k-1))+(r(k)-r(k-1))+(v(k)-v(k-1)))/(x(k)-x(k-1))-...
    trapz(x(1:k),(-alpha*u(1:k)))+sum(mass1(1:k))*dtheta1);
end
%Velocity Calculation RHS of Interface (RHS boundary as the anchor point)
for k=h+1:N 
w(k)=(1/total(k))* -((delta_z*(z(k)-z(k-1))+delta_r*(r(k)-r(k-1))+delta_v*(v(k)-v(k-1)))/(x(k)-x(k-1))+...
    Ksd*(z(k)+r(k)+v(k))*((z(k)-z(k-1))+(r(k)-r(k-1))+(v(k)-v(k-1)))/(x(k)-x(k-1))+sum(mass2(k:N+1))*dtheta2);
end
%Velocity is zero at the RHS Boundary
w(N+1)=0;

%Updating Node Position using the Exponential Time-stepping scheme
for i=1:N
    deltax(i)=x(i+1)-x(i);
    deltaw(i)=w(i+1)-w(i);
end
for i=1:N
    deltax_new(i)=deltax(i)*exp(dt*(deltaw(i)/deltax(i)));
end
for i=2:N
x_new(i)=x(1)+sum(deltax_new(1:i-1));    
end
%Boundary conditions on LHS & RHS
x_new(N+1)=x(N+1);
x_new(1)=x(1);
%Updating theta on LHS & RHS
theta1=theta1+dt*dtheta1; 
theta2=theta2+dt*dtheta2; 

%Updating the population densities
for i=2:h-2
 u_new(i)=(1/(x_new(i+1)-x_new(i)))*((x(i+1)-x(i))*u(i)+(dt)*(delta_u*(((u(i+1)-u(i))/(x(i+1)-x(i)))-((u(i)-u(i-1))/(x(i)-x(i-1))))+...
     Ksi*u(i+1)*((u(i+1)-u(i))+(z(i+1)-z(i))+(r(i+1)-r(i))+(v(i+1)-v(i)))/(x(i+1)-x(i))-Ksi*u(i-1)*((u(i)-u(i-1))+(z(i)-z(i-1))+(r(i)-r(i-1))+(v(i)-v(i-1)))/(x(i)-x(i-1))+...
     trapz(x(i:i+1),beta*u(i:i+1).*z(i:i+1)+epsilon*beta*u(i:i+1).*v(i:i+1)-alpha*u(i:i+1)-lambda*u(i:i+1))+w(i+1)*u(i+1)-u(i)*w(i)));
end
for i=h-1
 u_new(i)=(1/(x_new(i)-x_new(i-1)))*((x(i)-x(i-1))*u(i)+(dt)*(delta_u*(((u(i+1)-u(i))/(x(i+1)-x(i)))-((u(i)-u(i-1))/(x(i)-x(i-1))))+...
     Ksi*u(i+1)*((u(i+1)-u(i))+(z(i+1)-z(i))+(r(i+1)-r(i))+(v(i+1)-v(i)))/(x(i+1)-x(i))-Ksi*u(i-1)*((u(i)-u(i-1))+(z(i)-z(i-1))+(r(i)-r(i-1))+(v(i)-v(i-1)))/(x(i)-x(i-1))+...
     trapz(x(i-1:i),beta*u(i-1:i).*z(i-1:i)+epsilon*beta*u(i-1:i).*v(i-1:i)-alpha*u(i-1:i)-lambda*u(i-1:i))+w(i)*u(i)-u(i-1)*w(i-1)));
end
u_new(h:N+1)=0;

for i=2:h+1
 z_new(i)=(1/(x_new(i+1)-x_new(i)))*((x(i+1)-x(i))*z(i)+(dt)*(delta_z*(((z(i+1)-z(i))/(x(i+1)-x(i)))-((z(i)-z(i-1))/(x(i)-x(i-1))))+...
     Ksi*z(i+1)*(u(i+1)-u(i))/(x(i+1)-x(i))+Ksd*z(i+1)*((z(i+1)-z(i))+(r(i+1)-r(i))+(v(i+1)-v(i)))/(x(i+1)-x(i))-Ksi*z(i-1)*(u(i)-u(i-1))/(x(i)-x(i-1))-Ksd*z(i-1)*((z(i)-z(i-1))+(r(i)-r(i-1))+(v(i)-v(i-1)))/(x(i)-x(i-1))+...
     trapz(x(i:i+1),(-beta*u(i:i+1).*z(i:i+1)-rho*z(i:i+1)))+w(i+1)*z(i+1)-z(i)*w(i)));
end
for i=h+2:N 
 z_new(i)=(1/(x_new(i)-x_new(i-1)))*((x(i)-x(i-1))*z(i)+(dt)*(delta_z*(((z(i+1)-z(i))/(x(i+1)-x(i)))-((z(i)-z(i-1))/(x(i)-x(i-1))))+...
     Ksi*z(i+1)*(u(i+1)-u(i))/(x(i+1)-x(i))+Ksd*z(i+1)*((z(i+1)-z(i))+(r(i+1)-r(i))+(v(i+1)-v(i)))/(x(i+1)-x(i))-Ksi*z(i-1)*(u(i)-u(i-1))/(x(i)-x(i-1))-Ksd*z(i-1)*((z(i)-z(i-1))+(r(i)-r(i-1))+(v(i)-v(i-1)))/(x(i)-x(i-1))+...
     trapz(x(i-1:i),(-beta*u(i-1:i).*z(i-1:i)-rho*z(i-1:i)))+w(i)*z(i)-z(i-1)*w(i-1)));
end

for i=2:h+1
 r_new(i)=(1/(x_new(i+1)-x_new(i)))*((x(i+1)-x(i))*r(i)+(dt)*(delta_r*(((r(i+1)-r(i))/(x(i+1)-x(i)))-((r(i)-r(i-1))/(x(i)-x(i-1))))+...
     Ksi*r(i+1)*(u(i+1)-u(i))/(x(i+1)-x(i))+Ksd*r(i+1)*((z(i+1)-z(i))+(r(i+1)-r(i))+(v(i+1)-v(i)))/(x(i+1)-x(i))-Ksi*r(i-1)*(u(i)-u(i-1))/(x(i)-x(i-1))-Ksd*r(i-1)*((z(i)-z(i-1))+(r(i)-r(i-1))+(v(i)-v(i-1)))/(x(i)-x(i-1))+...
     +trapz(x(i:i+1),(lambda*u(i:i+1)))+w(i+1)*r(i+1)-r(i)*w(i)));
end
for i=h+2:N 
 r_new(i)=(1/(x_new(i)-x_new(i-1)))*((x(i)-x(i-1))*r(i)+(dt)*(delta_r*(((r(i+1)-r(i))/(x(i+1)-x(i)))-((r(i)-r(i-1))/(x(i)-x(i-1))))+...
     Ksi*r(i+1)*(u(i+1)-u(i))/(x(i+1)-x(i))+Ksd*r(i+1)*((z(i+1)-z(i))+(r(i+1)-r(i))+(v(i+1)-v(i)))/(x(i+1)-x(i))-Ksi*r(i-1)*(u(i)-u(i-1))/(x(i)-x(i-1))-Ksd*r(i-1)*((z(i)-z(i-1))+(r(i)-r(i-1))+(v(i)-v(i-1)))/(x(i)-x(i-1))+...
     trapz(x(i-1:i),(lambda*u(i-1:i)))+w(i)*r(i)-r(i-1)*w(i-1)));
end

for i=2:h+1
 v_new(i)=(1/(x_new(i+1)-x_new(i)))*((x(i+1)-x(i))*v(i)+(dt)*(delta_v*(((v(i+1)-v(i))/(x(i+1)-x(i)))-((v(i)-v(i-1))/(x(i)-x(i-1))))+...
     Ksi*v(i+1)*(u(i+1)-u(i))/(x(i+1)-x(i))+Ksd*v(i+1)*((z(i+1)-z(i))+(r(i+1)-r(i))+(v(i+1)-v(i)))/(x(i+1)-x(i))-Ksi*v(i-1)*(u(i)-u(i-1))/(x(i)-x(i-1))-Ksd*v(i-1)*((z(i)-z(i-1))+(r(i)-r(i-1))+(v(i)-v(i-1)))/(x(i)-x(i-1))+...
     trapz(x(i:i+1),(-epsilon*beta*v(i:i+1).*u(i:i+1)+rho*z(i:i+1)))+w(i+1)*v(i+1)-v(i)*w(i)));
end
for i=h+2:N 
 v_new(i)=(1/(x_new(i)-x_new(i-1)))*((x(i)-x(i-1))*v(i)+(dt)*(delta_v*(((v(i+1)-v(i))/(x(i+1)-x(i)))-((v(i)-v(i-1))/(x(i)-x(i-1))))+...
     Ksi*v(i+1)*(u(i+1)-u(i))/(x(i+1)-x(i))+Ksd*v(i+1)*((z(i+1)-z(i))+(r(i+1)-r(i))+(v(i+1)-v(i)))/(x(i+1)-x(i))-Ksi*v(i-1)*(u(i)-u(i-1))/(x(i)-x(i-1))-Ksd*v(i-1)*((z(i)-z(i-1))+(r(i)-r(i-1))+(v(i)-v(i-1)))/(x(i)-x(i-1))+...
     trapz(x(i-1:i),(-epsilon*beta*v(i-1:i).*u(i-1:i)+rho*z(i-1:i)))+w(i)*v(i)-v(i-1)*w(i-1)));
end

%Boundary conditions for the Updated population densities
u_new(1)=u_new(2);
z_new(1)=z_new(2);
r_new(1)=r_new(2);
v_new(1)=v_new(2);
u_new(h)=0;
z_new(N+1)=z_new(N);
r_new(N+1)=r_new(N);
v_new(N+1)=v_new(N);
u_new(h:N+1)=0;
%Update the population densities and node positions
u = u_new;
z = z_new;
r = r_new;
v = v_new;
x=x_new;
total(1:N+1)=u(1:N+1)+z(1:N+1)+r(1:N+1)+v(1:N+1);
     
        
if rem(t,25000)==0
figure(2)
plot(x,total,'b','LineWidth',1)
hold on
xlabel('x') 
ylabel('Total density') 
xlim([0 1])
ylim([0 1.1])

figure(3)
plot(x,u,'g','LineWidth',1)
xlabel('x') 
hold on
ylabel('Infected density') 
xlim([0 1])
ylim([0 1.1])

figure(4)
plot(x,z,'r','LineWidth',1)
hold on
xlabel('x') 
ylabel('Susceptible density')
xlim([0 1])
ylim([0 1.1])

figure(5)
plot(x,r,'m','LineWidth',1)
hold on
xlabel('x') 
ylabel('Recovered density') 
xlim([0 1])
ylim([0 1.1])

figure(6)
plot(x,v,'c')
hold on
xlabel('x') 
ylabel('Vaccinated density') 
xlim([0 1])
ylim([0 1.1])
 end
 end
