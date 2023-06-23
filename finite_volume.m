% dat
clear all
clf


MM = 64 ;  tend=4 ;  dtout=0.2 ;  factor=0.9 ;  D = 0.1 ;  a=0.0 ; b=1; t0=0;
tqfinal = tend*0.05;
tau1 = 1e6;
dx = 1/MM;
M  = (b-a) * MM;

x = zeros(M+2,1);

x = MESH(M,a,b,dx,x);

%% set the timestep
dtEXPL = dx*dx/(2*D);               % max timestep for stability
dt     = factor*dtEXPL;             % factor<1 increases stability
Nend   = round( (tend-t0)/dt ) +1;	% #of timesteps
tau1 = dt*1.1;
%% initialize
%%%%%%%%%%%%%%%%%%%         Semi-infinite bar         %%%%%%%%%%%%%%%%%%%
% [U,F]   = INIT_erf(M);
%%%%%%%%%%%%%%%%%%%         2 sides Dirichlet         %%%%%%%%%%%%%%%%%%%
% [U,F]   = INIT_const(M,dx,D);
%%%%%%%%%%%%%%%%%%%          2 sides Neumann          %%%%%%%%%%%%%%%%%%%
% [U,F]   = INIT_linear(M,dx,D,x);
%%%%%%%%%%%%%%%%%%%          2 sides MCV          %%%%%%%%%%%%%%%%%%%
[U,F,Q]   = INIT_zero(M,dx,D);
%%%%%%%%%%%%%%%%%%%          2 sides MCV          %%%%%%%%%%%%%%%%%%%
[U1,F1,Q1]   = INIT_zero(M,dx,D);


nsteps  = 1;
time    = 0.0;
tout    = max(dtout,dt);
uEXACT  = zeros(M+2,1);
ERR     = -1.0;
OUTPUT(x,U,uEXACT,time,nsteps,ERR,M)


%% begin time-stepping
for nsteps = 2:Nend
%%%%%%%%%%%%%%%%%%%         Semi-infinite bar         %%%%%%%%%%%%%%%%%%%
%     F = FLUX_erf(F,M,D,dx,U,time,b);
%     U =  PDE_erf(U,F,M,dt,dx,D,time,b);
%%%%%%%%%%%%%%%%%%%         2 sides Dirichlet         %%%%%%%%%%%%%%%%%%%
%     F =  FLUX_GammaD(F,M,D,dx,U);
%     U =  PDE_GammaD(U,F,M,dt,dx,D,time,b);
%%%%%%%%%%%%%%%%%%%          2 sides Neumann          %%%%%%%%%%%%%%%%%%%    
%     F =  FLUX_GammaN(F,M,D,dx,U);
%     U =  PDE_GammaN(U,F,M,dt,dx,D,time,b);
%%%%%%%%%%%%%%%%%%%          MCV         %%%%%%%%%%%%%%%%%%%    
    [F,Q]   =  FLUX_MCV(F,M,D,dt,dx,U,dt,Q,tqfinal,time);
    U =  PDE_MCV(U,F,M,dt,dx,D,time,b,Q,dt);
%%%%%%%%%%%%%%%%%%%          MCV         %%%%%%%%%%%%%%%%%%%    
    [F1,Q1] =  FLUX_MCV(F1,M,D,dt,dx,U1,tau1,Q1,tqfinal,time);
    U1 =  PDE_MCV(U1,F1,M,dt,dx,D,time,b,Q1,tau1);    
    
    if time >= tout
%%%%%%%%%%%%%%%%%%%         Semi-infinite bar         %%%%%%%%%%%%%%%%%%%
%         [ERR,uEXACT] = COMPARE_erf(M,D,x,time,U,uEXACT,ERR);        
%%%%%%%%%%%%%%%%%%%         2 sides Dirichlet         %%%%%%%%%%%%%%%%%%%
%         [ERR,uEXACT] = COMPARE_GammaD(M,D,x,time,U,uEXACT,ERR,dx);
%%%%%%%%%%%%%%%%%%%          2 sides Neumann          %%%%%%%%%%%%%%%%%%% 
%         [ERR,uEXACT] = COMPARE_GammaN(M,D,x,time,U,uEXACT,ERR,dx);
%%%%%%%%%%%%%%%%%%%          MCV          %%%%%%%%%%%%%%%%%%% 
%         [ERR,uEXACT] = COMPARE_GammaN(M,D,x,time,U,uEXACT,ERR,dx);
        
%         OUTPUT(x,U,uEXACT,time,nsteps,ERR,M);
        
        %% Plot
% %         figure
        hold on
        plot(x,U      ,'Color','black')
        plot(x,U1     ,'Color','red')
%         plot(x,uEXACT ,'Color','red'  )
%         legend('U','uEXACT','Location',"best")
%         xlabel 'Length'
%         ylim([0 0.02])
%         ylabel 'U(x,t)'  
        drawnow
%         box on
%         grid on
%         hold off
    tout = tout + dtout;
    end
    time = nsteps*dt;
end % end of time-stepping

% exportgraphics(gcf,'2Neumann.pdf','ContentType','vector')

if time >= tend
    fprintf(1, 'DONE, at time= %f, nsteps = %d \n', time,nsteps);
    fprintf(1, '    , max error = %e \n',           ERR);
else
    fprintf(1, '... out of timesteps: need bigger Nend \n')
end
%%
function x = MESH(M,a,b,dx,x)
x(1) = a;
x(2) = x(1) + dx/2;
for i = 3:(M+1)
    x(i) = x(2) + (i-2)*dx;
end
x(M+2) = b;
end
%%
function OUTPUT(x,U,uEXACT,time,nsteps,ERR,M)
filename = [num2str(time),'.o.prof'];
fileID = fopen(filename,'w');   % output file
fprintf(fileID, '# Profile at time: %9.4f, nsteps=%6d \n', time,nsteps);
fprintf(fileID, '# Error up to this time: %15.6e \n',ERR);
for i = 1:(M+2)
    fprintf(fileID, '%12.4f %16.8g %16.8g \n', x(i), U(i), uEXACT(i));
end
fclose(fileID); % end of the output
end

%% zero bc
function [U,F,Q] = INIT_zero(M,dx,D)
U = zeros(M+2,1);
F = zeros(M+1,1);
Q = zeros(M+1,1);
U(1) = 0;
U(M+2)  = 0;
end
%% FLUX for MCV
function [F,Q] = FLUX_MCV(F,M,D,dt,dx,U,tau,Q,tqfinal,time)
if time < tqfinal
F(1)    = 1;                                % Neumann BC for b
else
F(1)    = 0;                                % Neumann BC for b    
end
for i = 2:M
%     F(i)= -D*(U(i+1) - U(i)  )/(dx);        % Internal face values
%     Q(i)= Q(i) + dt/tau*(F(i)-Q(i));
    F(i)= -D*(F(i)-dt/tau*F(i)(U(i+1) - U(i))/(dx));        % Internal face values
end
Q(M+1)  = -D*(U(M+2) - U(M+1))/(dx/2);
end

%% PDE function for MCV
% 
function U = PDE_MCV(U,F,M,dt,dx,D,time,b,Q,tau)
% U(1)    = U(2)+dx/(2*D)*Q(1);              % Neumann BC for b
for i = 2:(M+1)
    U(i) = U(i) + dt/dx/tau*(Q(i-1)-Q(i));      % Internal U values
end
% U(M+2)  = U(M+1)+dx/(2*D)*Q(M+1);%0;                                % Dirichlet BC for b
end

%% zero bc
function [U,F] = INIT_ND(M,dx,D)
U = zeros(M+2,1);
F = zeros(M+1,1);
U(1) = 0;
U(M+2)  = 0;
end
%% FLUX for ND
function F = FLUX_ND(F,M,D,dt,dx,U,time,b,tau)
U(M+2)  = 0;                                % Dirichlet BC for b
U(1)    = U(2);                             % Neumann BC for b
F(1)    = 1;                                % Neumann BC for b
for i = 2:(M+1)
    F(i)= -D*(U(i+1) - U(i)  )/(dx);        % Internal face values
end
F(M+1)  = -D*(U(M+2) - U(M+1))/(dx/2);      % Dirichlet BC for b
end

%% PDE function for ND
% 
function U = PDE_ND(U,F,M,dt,dx,D,time,b)
U(1)    = U(2)+dx/(2*D)*F(1);              % Neumann BC for b
for i = 2:(M+1)
    U(i) = U(i) + dt/dx*(F(i-1)-F(i));      % Internal U values
end
U(M+2)  = 0;                                % Dirichlet BC for b
end
%% Initialize function for linear temperature distribution
% with one Dirichlet BC
function [U,F] = INIT_erf(M)
U = zeros(M+2,1);   % Initialize U vector
F = zeros(M+1,1);   % Initialize Flux vector
U(1) = 1;           % Dirichlet BC for a
end

%% Flux function for semi-infinite bar
% 
function F = FLUX_erf(F,M,D,dx,U,time,b)
U(M+2)  = erfc(b/2/sqrt(D*time));           % Dirichlet BC for b
F(1)    = -D*(  U(2) -  U(1) )/(dx/2);      % Face value at a
for i = 2:(M+1)
    F(i)= -D*(U(i+1) - U(i)  )/(dx);        % Internal face values
end
F(M+1)  = -D*(U(M+2) - U(M+1))/(dx/2);      % Face value at b
end

%% PDE function for semi-infinite bar
% 
function U = PDE_erf(U,F,M,dt,dx,D,time,b)
for i = 2:(M+1)
    U(i) = U(i) + dt/dx*(F(i-1)-F(i));      % Internal U values
end
U(M+2)  = erfc(b/2/sqrt(D*time));           % Dirichlet BC for b
end

%% Compare with exact solution: u(x,t) = ERFC(0.5*x/SQRT(D*t))
% of:    u_t = D u_xx , 0<x<Inf , t>0 ; u(x,0)=0 ; u(0,t)=1.
function [ERR,uEXACT] = COMPARE_erf(M,D,x,time,U,uEXACT,ERR)
% ERR = 0.0; % Comment to see  max error up to this time
             % Uncomment to see max error at only at tout 
for i = 1:(M+2)
    arg      = 0.5*x(i)/sqrt(D*time);
    uEXACT(i)= erfc(arg);
    ERRi= abs(U(i)-uEXACT(i));
    ERR = max(ERRi,ERR);
end
end

%% Initialize function for linear temperature distribution
% 
function [U,F] = INIT_linear(M,dx,D,x)
U = zeros(M+2,1);
F = zeros(M+1,1);
for i = 1:M+2
U(i) = x(i);
end
end

%% Flux function for two Neumann BC
% 
function F = FLUX_GammaN(F,M,D,dx,U)
F(1)    = 0;
for i = 2:(M+1)
    F(i) = -D*(U(i+1) - U(i)  )/(dx);
end
F(M+1)  = 0;
end

%% PDE function for two Neumann BC
% 
function U = PDE_GammaN(U,F,M,dt,dx,D,time,b)
U(1) = U(2);
for i = 2:(M+1)
    U(i) = U(i) + dt/dx*(F(i-1)-F(i));
end
U(M+2) = U(M+1);
end

%% Compare function for two Neumann BC
% for exact solution of 
% u_t = - D u_xx , 0<x<L , t>0 ; q(0,t)=0 ; q(L,t)=0 ; u(x,0)=x.
function [ERR,uEXACT] = COMPARE_GammaN(M,D,x,time,U,uEXACT,ERR,dx)
% ERR = 0.0; % Comment to see  max error up to this time
             % Uncomment to see max error at only at tout 
for i = 1:(M+2)
arg = 0;
for j = 0:3
    arg = arg + 1/((2*j+1)^2)*exp(-D*((2*j+1)^2)*pi^2*time)*cos((2*j+1)*pi*x(i));
end
    uEXACT(i)= 0.5 - 4/pi^2*arg;
    ERRi= abs(U(i)-uEXACT(i));
    ERR = max(ERRi,ERR);
end
end

%% Initialize function for constant temperature distribution
% Set Dirichlet BC for both boundaries
function [U,F] = INIT_const(M,dx,D)
U = zeros(M+2,1);
F = zeros(M+1,1);
U(1) = 1;
U(M+2)  = 0;
end

%% Flux function for two Dirichlet BC
% 
function F = FLUX_GammaD(F,M,D,dx,U)
F(1)    = -D*(  U(2) -  U(1) )/(dx/2);
for i = 2:(M+1)
    F(i) = -D*(U(i+1) - U(i)  )/(dx);
end
F(M+1)  = -D*(U(M+2) - U(M+1))/(dx/2);
end

%% PDE function for two Dirichlet BC
% 
function U = PDE_GammaD(U,F,M,dt,dx,D,time,b)
for i = 2:(M+1)
    U(i) = U(i) + dt/dx*(F(i-1)-F(i));
end
end

%% Compare function for two Dirichlet BC
% for exact solution of 
% u_t = - D u_xx , 0<x<L , t>0 ; u(0,t)=0 ; u(L,t)=1 ; u(x,0)=0.
function [ERR,uEXACT] = COMPARE_GammaD(M,D,x,time,U,uEXACT,ERR,dx)
% ERR = 0.0; % Comment to see  max error up to this time
             % Uncomment to see max error at only at tout 
for i = 1:(M+2)
arg = 0;
for j = 1:10
    arg = arg + ((-1)^j)/j*exp(-D*j^2*pi^2*time)*sin(j*pi*x(i));
end
    uEXACT(i)= 1*x(i) + 2/pi*arg;
    ERRi= abs(U(i)-uEXACT(i));
    ERR = max(ERRi,ERR);
end
end
