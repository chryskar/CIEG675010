function [p1,p2,p3] = covid_sir(S_0,I_0,R_0,k1,b1,max_days,figure_num)
% function covid_sir(S_0,I_0,R_0,k1,b1,max_days,figure_num)

% Description: This function applies a SIR model to simulate the behavior
% of COVID-19

% INPUTS: 
% S_0: Initial susceptible rate
% I_0: Initial infected rate (very small)
% R_0: Initial recovered rate 
% k1: Average time to recover
% b1: How often a susceptible-infected contact occurs
% max_days: Number of days for which the system should be simulated for
% figure_num: Number of figure to be utilized for plotting

% OUTPUTS:
% p1: variable linked to the first curve plotted on the figure (S)
% p2: variable linked to the first curve plotted on the figure (I)
% p3: variable linked to the first curve plotted on the figure (R)

k = 1/k1; % Define the average period of infectiousness
b = 1/b1; % Frequency of susceptible-infected contacts

syms x(t) [3 1] real; % Define three symbolic variable, each for one the population groups
% Let S = x1(t), I = x2(t) and R = x3(t)
%Eq.(3.1)
eqn1 = diff(x1(t),t) == -b*x1*x2;
%Eq.(3.2)
eqn2 = diff(x2(t),t) == b*x1*x2 - k*x2;
%Eq.(3.3)
eqn3 = diff(x3(t),t) == k*x2;

% Group all equations together
eqns = [eqn1, eqn2, eqn3];
% Group all variables together
vars = [x1(t) x2(t) x3(t)];

% Determine the mass matrix form for the ode45 solver
[M,F] = massMatrixForm(eqns,vars);
f = M\F;

% Set initial condition for the ode45 solver
init_conds = [S_0 I_0 R_0];

% Group together all information about the ODE
odefun = odeFunction(f,vars);
t_total=[]; %In this vector we will store all time variables that the ode45 solver returns
sol=[];     %In this vector we will store all state variables that the ode45 solver returns

i=1;    %Iteration counter for the loop
i_th = max_days; %Maximum Number of Iterations (100 days)
opts = odeset('RelTol',1e-6,'AbsTol',1e-6); %Absolute and Relative error tolerances for ode45
max_step_size = 1; % Time period for which the ode is solved at each step
while(1) % Iterate forever until a certain condition is met
    [t_i,sol_i] = ode45(odefun,[0 max_step_size],init_conds,opts); %Solve the system's equations using the specified ICs for 1 day with the error tolerances "opts" specified earlier
    t_total=[t_total;t_i+((i-1)*max_step_size)];    %Adjust the local time vector, so that it corresponds to the total time elapsed and append it to the total time vector
    sol=[sol;sol_i];    %Append the local solutions for the state to the total state solutions vector
    init_conds=[sol(end,1); sol(end,2); sol(end,3)];    %Set the last state of the system for this iteration, as ICs for the next one
    i = i + 1; % Increase the iteration counter
    if i > i_th % Stop simulation when the desired number of days has been reached
        break
    end
end

fig7 = figure(figure_num); % Open a new figure
hold on; % Keep all plots at the same figure
% Plot the time evolution of the susceptible rate
p1 = plot(t_total(1:50:end),sol(1:50:end,1),'--','color',[1, 1-b1/10, 0]);
% Plot the time evolution of the infected rate
p2 = plot(t_total(1:50:end),sol(1:50:end,2),'-','color',[1, 1-b1/10, 0]);
% Plot the time evolution of the recovered rate
p3 = plot(t_total(1:50:end),sol(1:50:end,3),'-.','color',[1, 1-b1/10, 0]);
end