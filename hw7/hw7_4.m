
close all
clc


% INPUT PARAMETERS   =====================================================
Ne = 800; Ni = 200;         % number of excitatory & inhibitory neurons
Nt = 8000;                  % number of time steps


% Model parameters   =====================================================
re = rand(Ne,1); ri = rand(Ni,1);     % random numbers
a = [0.02 * ones(Ne,1); 0.02 + 0.08*ri];
b = [0.20   *ones(Ne,1); 0.25 - 0.05*ri];
c = [-65+15*re.^2; -65*ones(Ni,1)];
d = [8-6*re.^2   ; 2*ones(Ni,1)];
g1=rand(Ne+Ni,Ne); 
g2=rand(Ne+Ni,Ni);
S_exc=0.625;S_inh=-0.75;
S = [S_exc*g1, S_inh*g2];   % coupling strengths
Cc=0;Ce=0;Sm=0;Gm=00;K_Glio=1;Glio=[];
f=ones(1000,1);freq=zeros(1000,1);

v = -65*ones(Ne+Ni,1); % Initial values of v
u =  b.*v;             % Initial values of u
firings = [];          % spike timings


% Time Evolution of Systems  ==============================================
%%
k_e=[];k_i=[];
S_e=[];S_i=[];

for t = 1:Nt 
    if t==2000
        Cc=0;Ce=0;Sm=0;Gm=00;
    end
    if t<=2000
       S_exc=0.5;S_inh=-1;
       K_exc=5;K_inh=2;
    end
    if t==3000
        Cc=0;Ce=0;Sm=0;Gm=00;
    end
    if t>=2000 && t<=3000
       K_exc=K_exc+2.5/1000;
       K_inh=K_inh-1/1000;
    end
    if t==6000
        Cc=0;Ce=0;Sm=0;Gm=00;
    end
    if t>=3000 && t<=6000
       K_exc=K_exc-(7.5-6.25)/3000;
       K_inh=K_inh+0.5/3000;
    end
    if t>=3000 && t<=3300
       S_exc=S_exc+(0.625-0.5)/300;
       S_inh=S_inh+(1-0.75)/300;
    end

    if t>=6000 && t<=6500
       K_exc=K_exc-(6.25-3.75)/500;
       K_inh=K_inh+(2.5-1.5)/500;
       S_exc=S_exc-(0.625-0.56)/500;
       S_inh=S_inh-(0.875-0.75)/500;
    end
    
    S = [S_exc*g1, S_inh*g2];   % coupling strengths

    delta=0.8;
    f1=randn(Ne,1);f2=randn(Ni,1);
    I = [delta*K_exc*f1;delta*K_inh*f2]; % thalamic input
    fired = find(v>=30);               % indices of spikes
    firings = [firings; t+0*fired,fired];  % time steps  / fired neurons

    v(fired) = c(fired);               % membrane potential
    u(fired) = u(fired)+d(fired);      % recovery potential

    freq(fired)=freq(fired)+f(fired);

    I = I + sum(S(:,fired),2);           % thalamic + synaptic input

    Cc=Cc+20*(-Cc-(2/0.04)*(0.13*(Cc^2/(1+Cc^2))-(Ce^2/(1+Ce^2))*((Cc^4/(0.9^4+Cc^4)))-0.004*Ce)+ (0.31 + 0.006 * Sm))*0.001;
    Ce=Ce+20*25*(0.13*(Cc^2/(1+Cc^2))-(Ce^2/(1+Ce^2))*((Cc^4/(0.9^4+Cc^4)))-0.004*Ce)*0.001;
    Sm=Sm+0.2*((1+tanh(5*(K_exc*norm(f1)-0.45)))*(1-Sm)-Sm/3)*0.001;
    Gm=Gm+40*((1+tanh(10*(Cc-0.5)))*(1-Gm)-Gm/3)*0.001;
    Glio=[Glio,Gm];
    I_Astro=K_Glio*delta*Gm;
    
    v = v+0.5*(0.04*v.^2+5*v+140-u+I+I_Astro); % HALF-STEP: step 0.5 ms
    v = v+0.5*(0.04*v.^2+5*v+140-u+I+I_Astro); %  for numerical stability
    u = u+a.*(b.*v-u); 

    %vE(t) = v(numE);                    % excitatory neuron 
    %vI(t) = v(numI);                    % inhibitory neuron
    k_e=[k_e,K_exc];
    k_i=[k_i,K_inh];
    S_e=[S_e,S_exc];
    S_i=[S_i,S_inh];

end



vE(vE > 30) = 30;
vI(vI > 30) = 30;
minf=min(freq);
tS = 1:Nt; nF = zeros(Nt,1);        % fired neurons at each time step
for t = 1 : Nt
nF(t) = sum((firings(:,1) == t));
end


% GRAPHICS ===============================================================
%%
figure(1)    % raster plot of firong neurons / % firing neurons
set(gcf,'units','normalized','Position',[0.1 0.4 0.32,0.45]);
set(gca,'fontsize',12);
col = [0 0 1];
hP = subplot(1,1,1);   % raster plot
plot(firings(:,1)/1000,firings(:,2),'.','color',col);
xlabel('Time(s)','fontsize',12);
set(hP, 'Position',[0.1300 0.52 0.7750 0.4]);
ylabel('# of neurons','fontsize',12)

subplot(2,1,2);       % percentage of firong neurons each time step
t0=0.001:0.001:8;
pspectrum(nF,t0,'spectrogram','MinThreshold',0,'TimeResolution',0.15,'FrequencyLimits',[0,80]);

colormap jet
colorbar
title('')





