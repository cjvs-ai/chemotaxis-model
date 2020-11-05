%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Program to (1) simulate chemotaxis under presence of a %%%
%%% chemoattractant by solving a Keller-Segel type model   %%%
%%% and (2) use this model to fit to experimental data and %%%
%%% extract parameters for best fit. Here we use noise     %%%
%%% superimposed on model to substitute for experimental   %%%
%%% data.                                                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except noisesol Bdata % Generate noise solution to fit only on the first run to ensure reproducability
% clf

tic

m = 0; % Cartesian geometry
x = linspace(0,0.003,75); % Spatial length 1mm
t = linspace(0,800,75); % 10 minute simulation time
sol = pdepe(m,@pde,@pdeic,@pdebc,x,t); % Compute solution to K-S model using pdepe
%%%
figure(1) % Plot surface plot + initial and final bacterial density profile
subplot(2,1,1), surf(x,t,sol);
subplot(2,1,2), plot(x,sol([1,75],:));
%%%

% % Add Gaussian noise to signal to test the fitting
% noisesol = awgn(sol, 15, 'measured');
%%%
% figure(2) % Plot noise
% plot(x,noisesol(100,:));
%%%
% % Late time noise data to fit to
% Bdata = noisesol(100,:);
% % Late time data without noise to test fit on perfect data
% Bdatanonoise = sol(100,:);

% Initial guess for parameters to fit (mu, V, chi_0)
param = [1E-10, 10E-6, 1E-8];

% % Generate random numbers for inital guess between an interval
% mu = 5.9E-10; % Random motility coefficient
% V = 30E-6; % Bacteria swimming speed ms^-1
% chi_0 = 5E-8; % Chemotactic sensitivity m^2 s^-1
% random_mu = ((1E-10+mu)-mu).*rand()+mu;
% random_V = ((1E-6+V)-V).*rand()+V;
% random_chi_0 = ((1E-8)-chi_0).*rand()+chi_0;
% param = [random_mu, random_V, random_chi_0];

% Call sse function with initial param to test if it works
% solparam = pdeparam(x,t,param);
%%%
% figure(3)
% plot(x,solparam(100,:)); hold on; plot(x,Bdata);
%%%

% Optimise parameters using sse to find best fit using fminsearch
% and ignoring tail of data
fun = @(param)sseval_notail(param,x,t,Bdata);
bestparam = fminsearch(fun,param); % Parameters which give best fit

% Sum of squared errors with final params
% sse = sseval_notail(bestparam,x,t,Bdata);

% Plot model with best fit parameters
testsol = pdeparam(x,t,bestparam);
%%%
figure(1)
plot(x,testsol(100,:)); hold on; plot(x,Bdata);
xlabel('x (m)')
ylabel('Disribution B(x)')
legend('Best-fit Solution', 'Simulated Data')
%%%

toc

%%% Sum of squared errors function %%%
% Calls pdeparam to solve model with variable parameters
function sse = sseval(param,x,t,Bdata)
solparam = pdeparam(x,t,param);
sse = sum((Bdata - solparam(100,:)).^2);
end

%%% Sum of squared errors which ignores the tail %%%
function sse = sseval_notail(param,x,t,Bdata)
solparam = pdeparam(x,t,param);
late_solparam = solparam(100,:);
sse = sum((Bdata(200:300) - late_solparam(200:300)).^2);
end

%%% PDE solver for K-S chemotaxis model with adjustable input parameters (vector: param) %%%
function solparam = pdeparam(x,t,param)
solparam = pdepe(0,@pde,@pdeic,@pdebc,x,t);

function [c,f,s] = pde(x,t,B,dBdx)
% Constants
% mu = 5.9E-10; % Random motility coefficient
% V = 30E-6; % Bacteria swimming speed ms^-1
% chi_0 = 5E-8; % Chemotactic sensitivity m^2 s^-1
K_D = 0.125; % Receptor / ligand dissociation constant mM
dCdx = 83; % Chemoattractant concentration gradient = 1000*0.0833333
C = 0.1; % Chemoattractant mM
W = 0.001; % Width in m
v_c = (8*param(2))/(3*pi)*tanh((param(3)*pi*K_D)*dCdx/(8*param(2)*(K_D + (C*x/W))*(K_D + (C*x/W)))); % Net drift velocity
% Bacterial transport equation in c f s components
c = 1;
f = param(1) * dBdx - v_c * B;
s = 0;
end

% Initial condition
function B0 = pdeic(x)
B0 = 1;
end

% Boundary conditions
function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
pl = 0;
ql = 1; %1/mu + v_cl*ul;
pr = 0;
qr = 1; %1/mu + v_cr*ur;
end

end

%%% PDE Solver for Keller-Segel model of chemotaxis %%%
% Define pde to solve numerically
function [c,f,s] = pde(x,t,B,dBdx)
% Constants
mu = 5.9E-10; % Random motility coefficient
V = 15E-6; % Bacteria swimming speed ms^-1
chi_0 = 5E-8; % Chemotactic sensitivity m^2 s^-1
K_D = 0.125; % Receptor / ligand dissociation constant mM
dCdx = 33.3; % Chemoattractant concentration gradient = 1000*0.0833333
C = 0.1; % Chemoattractant mM
W = 0.003; % Width in m
v_c = (8*V)/(3*pi)*tanh((chi_0*pi*K_D)*dCdx/(8*V*(K_D + (C*x/W))*(K_D + (C*x/W)))); % Net drift velocity
% Bacterial transport equation in c f s components
c = 1;
f = mu * dBdx - v_c * B;
s = 0;
end
% Define initial condition
function B0 = pdeic(x)
B0 = 1;
end
% Define boundary conditions
function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
mu = 5.9E-10; % Random motility coefficient
V = 15E-6; % Bacteria swimming speed ms^-1
chi_0 = 5E-8; % Chemotactic sensitivity m^2 s^-1
K_D = 0.125; % Receptor / ligand dissociation constant mM
dCdx = 33.3; % Chemoattractant concentration gradient = 1000*0.0833333
C = 0.1; % Chemoattractant mM
W = 0.003; % Width in m

v_cl = (8*V)/(3*pi)*tanh((chi_0*pi*K_D)*dCdx/(8*V*(K_D+(C*xl/W))*(K_D+(C*xl/W)))); % Net drift velocity left
v_cr = (8*V)/(3*pi)*tanh((chi_0*pi*K_D)*dCdx/(8*V*(K_D+(C*xr/W))*(K_D+(C*xr/W)))); % Net drift velocity right

pl = 0;
ql = 1; %1/mu + v_cl*ul;
pr = 0;
qr = 1; %1/mu + v_cr*ur;
end

