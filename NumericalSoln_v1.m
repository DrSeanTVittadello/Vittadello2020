% Script for numerical simulations of the functional differential equation
% model in the article Vittadello et al. (2020), A novel mathematical model
% of heterogeneous cell proliferation.
%
% Author: Sean T. Vittadello
%         sean.vittadello@qut.edu.au
%         School of Mathematical Sciences
%         Queensland University of Technology
%         Brisbane, Australia
%
% Last update: 26 February 2020

%% Set the default text interpreter
set(0,'defaultTextInterpreter','latex');

%% Set the parameters for the mathematical model
US = 200; % Maximum cell cycle duration for the slow-proliferating cells (h)
UF = 200; % Maximum cell cycle duration for the fast-proliferating cells (h)
rS = 0.1; % Intrinsic growth rate for the slow-proliferating cells (h^{-1})
rF = 1; % Intrinsic growth rate for the fast-proliferating cells (h^{-1})
alphaS = 1; % Proportion of cell divisions that are symmetric for slow-proliferating cells
alphaF = 1; % Proportion of cell divisions that are symmetric for fast-proliferating cells
bS = 0.1; % Induced switching rate from slow- to fast-proliferating (h^{-1})
bF = 0; % Induced switching rate from fast- to slow-proliferating (h^{-1})
K = 500; % Carrying-capacity density

%% Truncate the Erlang distribution
% Set the parameters for the Erlang distributions (not truncated)
lambdaS = 0.08; % Rate of Erlang distribution for slow-proliferating cells (h^{-1})
kS = 8; % Shape of Erlang distribution for slow-proliferating cells
lambdaF = 0.8; % Rate of Erlang distribution for fast-proliferating cells (h^{-1})
kF = 8; % Shape of Erlang distribution for fast-proliferating cells (h^{-1})

% Calculate the cumulative distribution function (CDF) at the truncation points of the Erlang distributions 
syms n
CDF_S = double((1 - symsum((1/factorial(n))*exp(-lambdaS*US)*(lambdaS*US)^n,n,0,kS-1)));
CDF_F = double((1 - symsum((1/factorial(n))*exp(-lambdaF*UF)*(lambdaF*UF)^n,n,0,kF-1)));

% Probability density function (PDF) for the Erlang distribution
Erlang = @(x,L,k) (L^k * x.^(k-1) .* exp(-L*x))/(factorial(k-1));

% Truncated Erlang distributions
Erlang_trunc_S = @(x) Erlang(x,lambdaS,kS)/CDF_S;
Erlang_trunc_F = @(x) Erlang(x,lambdaF,kF)/CDF_F;

%% Set the parameters for the numerical procedures
% Forward Euler method
h = 0.1; % Time step size (h)
T = 200; % Total time (h)
steps = T/h; % Number of steps
S_hist_fun = @(x) 100*exp(rS*x); % History function for slow-proliferating cells for -US \le x \le 0
F_hist_fun = @(x) 0.0001*exp(rF*x); % History function for fast-proliferating cells for -UF \le x \le 0

% Trapezoidal rule (uniform grid)
N_S = 1000; % Number of subintervals corresponding to slow-proliferating cells
length_S = US/N_S; % Length of subintervals corresponding to slow-proliferating cells
N_F = 1000; % Number of subintervals corresponding to fast-proliferating cells
length_F = UF/N_F; % Length of subintervals corresponding to fast-proliferating cells

%% Numerical solution
% Define the required arrays
S_TR = zeros(N_S + 1,1); % Interpolated solution values for trapezoidal rule - slow-proliferating cells
S_TR_t = zeros(N_S + 1,1); % Time points for trapezoidal rule - slow-proliferating cells
ErlangS_t = zeros(N_S + 1,1); % Time points for truncated Erlang distribution (trapezoidal rule) - slow-proliferating cells
F_TR = zeros(N_F + 1,1); % Interpolated solution values for trapezoidal rule - fast-proliferating cells
F_TR_t = zeros(N_F + 1,1); % Time points for trapezoidal rule - fast-proliferating cells
ErlangF_t = zeros(N_F + 1,1); % Time points for truncated Erlang distribution (trapezoidal rule) - fast-proliferating cells
    
S = zeros(steps+1,1); % Cell density of slow-proliferating cells at time t
S(1,1) = S_hist_fun(0); % Cell density of slow-proliferating cells at time t = 0 h
S_hist = zeros(steps+1,1); % Cell density of slow-proliferating cells at previous times
F = zeros(steps+1,1); % Cell density of fast-proliferating cells at time t
F(1,1) = F_hist_fun(0); % Cell density of fast-proliferating cells at time t = 0 h
F_hist = zeros(steps+1,1); % Cell density of fast-proliferating cells at previous times
time = (0:h:T)'; % time points for forward Euler method

t = 0; % Initialise the time
for i=1:steps
    % Set the time points for the trapezoidal rule
    S_TR_t = (t-US:length_S:t)';
    F_TR_t = (t-UF:length_F:t)';
    
    % Set the time points for the truncated Erlang distribution
    ErlangS_t = t - S_TR_t;
    ErlangF_t = t - F_TR_t;
    
    for j=1:size(S_TR_t,1)
        if S_TR_t(j,1) <= 0
            S_TR(j,1) = S_hist_fun(S_TR_t(j,1));
        else
            index = find(time >= S_TR_t(j,1),1);
            if S_TR_t(j,1) == time(index,1)
                S_TR(j,1) = S(index,1);
            elseif S_TR_t(j,1) == time(index-1,1)
                S_TR(j,1) = S(index-1,1);
            else
                S_TR(j,1) = pchip([time(index-1,1),time(index,1)],[S(index-1,1),S(index,1)],S_TR_t(j,1));
            end
        end
    end
    
    for j=1:size(F_TR_t,1)
        if F_TR_t(j,1) <= 0
            F_TR(j,1) = F_hist_fun(F_TR_t(j,1));
        else
            index = find(time >= F_TR_t(j,1),1);
            if F_TR_t(j,1) == time(index,1)
                F_TR(j,1) = F(index,1);
            elseif F_TR_t(j,1) == time(index-1,1)
                F_TR(j,1) = F(index-1,1);
            else
                F_TR(j,1) = pchip([time(index-1,1),time(index,1)],[F(index-1,1),F(index,1)],F_TR_t(j,1));
            end
        end
    end
    
    % Estimates of the distributed delays at the current time
    delayS = length_S*(0.5*(S_TR(1,1)*Erlang_trunc_S(US) + S_TR(end,1)*Erlang_trunc_S(0)) + sum(S_TR(2:end-1,1).*Erlang_trunc_S(ErlangS_t(2:end-1,1)),1));
    delayF = length_F*(0.5*(F_TR(1,1)*Erlang_trunc_F(UF) + F_TR(end,1)*Erlang_trunc_F(0)) + sum(F_TR(2:end-1,1).*Erlang_trunc_F(ErlangF_t(2:end-1,1)),1));
    
    % Estimates of the solutions at the current time
    S(i+1,1) = S(i,1) + h*(((2*alphaS-1)*rS*delayS + 2*(1-alphaF)*rF*delayF)*(1-(S(i,1)+F(i,1))/K) - ((bS-bF)/K)*S(i,1)*F(i,1));
    F(i+1,1) = F(i,1) + h*((2*(1-alphaS)*rS*delayS + (2*alphaF-1)*rF*delayF)*(1-(S(i,1)+F(i,1))/K) + ((bS-bF)/K)*S(i,1)*F(i,1));
    
    t = i*h; % Next time point
end
%% Plot solutions
figure('name','Numerical simulations')
plot(time,S,'LineStyle','-','LineWidth',1.5,'Marker','none','Color','r')
title(sprintf('S(t) and F(t)'));
xlabel('$t$ (h)');
ylabel('$S(t)$, $F(t)$');
xlim([0 T]);
ylim([0 K]);
xticks([0 T])
yticks([0 K])
title(sprintf('$\\alpha_S = %g$, $\\alpha_F = %g$, $\\beta_S = %g$, $\\beta_F = %g$',alphaS,alphaF,bS,bF),'Interpreter','latex')
set(gca,'FontSize',20)
hold on
plot(time,F,'LineStyle','-','LineWidth',1.5,'Marker','none','Color','b')

% maximum value of the Erlang PDF is at t = (k-1)/L
tmaxS = (kS-1)/lambdaS; 
tmaxF = (kF-1)/lambdaF;

figure('name','Truncated Erlang PDF')
plot(0:0.1:US,Erlang_trunc_S(0:0.1:US),'LineStyle','-','LineWidth',1.5,'Marker','none','Color',[1,0.55,0])
xlabel('Duration, $t$ (h)');
ylabel('Probability density');
xlim([0 US]);
if UF < US
    xticks([0 UF US])
elseif UF == US
    xticks([0 UF])
end
ylim([0 ceil(100*max(Erlang_trunc_S(tmaxS),Erlang_trunc_F(tmaxF)))/100]);
yticks([0 ceil(100*max(Erlang_trunc_S(tmaxS),Erlang_trunc_F(tmaxF)))/100])
title(sprintf('PDF of truncated Erlang distribution for slow and fast proliferators'))
set(gca,'FontSize',20)
hold on
plot(0:0.1:US,Erlang_trunc_F(0:0.1:US),'LineStyle','-','LineWidth',1.5,'Marker','none','Color',[0,0.81,0.82])

figure('name','Erlang PDF')
plot(0:0.1:US,Erlang(0:0.1:US,lambdaS,kS),'LineStyle','-','LineWidth',1.5,'Marker','none','Color',[1,0.55,0])
xlabel('Duration, $t$ (h)');
ylabel('Probability density');
xlim([0 US]);
if UF < US
    xticks([0 UF US])
elseif UF == US
    xticks([0 UF])
end
ylim([0 ceil(100*max(Erlang(tmaxS,lambdaS,kS),Erlang(tmaxF,lambdaF,kF)))/100]);
yticks([0 ceil(100*max(Erlang(tmaxS,lambdaS,kS),Erlang(tmaxF,lambdaF,kF)))/100])
title(sprintf('PDF of Erlang distribution for slow and fast proliferators'))
set(gca,'FontSize',20)
hold on
plot(0:0.1:US,Erlang(0:0.1:US,lambdaF,kF),'LineStyle','-','LineWidth',1.5,'Marker','none','Color',[0,0.81,0.82])
