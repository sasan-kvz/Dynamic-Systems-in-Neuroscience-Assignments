%%hw1
%%hw1_1
clc; close all ; clear

vcr = linspace(-150,50,200);
resa = zeros(size(vcr));
resa2 = zeros(size(vcr));

for i = 1:length(vcr)
    
    tspan = [0, 60];
    y0 = zeros(3,1);
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
    
    Vc = -90;
    Vs = vcr(i);
    [t,y] = ode45(@(t,y) ode(t,y,Vc,Vs), tspan, y0,opts);
    
    n = y(:,1);
    m = y(:,2);
    h = y(:,3);
    
    idx = find(t<20,1,'last');
    V = zeros(size(t));
    V(1:idx) = Vc;
    V(idx+1:end) = Vs;
    
    I = 120 .* (m.^3).* h .* (V-120) + 367 .* (n.^4) .* (V+12) + (0.3) .* (V-10.6);
    
    resa(i) = I(end);
    resa2(i) = min(I(idx+1:end));
end
hold on
plot(vcr,resa,'r','LineWidth',3)
plot(vcr,resa2,'b','LineWidth',3)
grid on; hold off;set(gca,'FontSize',10);
xlabel('V [mV]'); ylabel('I [pA]'); legend('I-inf','I0','k','r');


%% hw1_2
clc; close all ; clear


tspan1 = [0, 15];
y0 = zeros(2,1);
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);

a=25; %% we will change it to see when neuron will fire
[t1,y1] = ode45(@(t,y) odehh(t,y,a), tspan1, y0,opts);

subplot 211
plot(t1,y1(:,1),'b','LineWidth',2)
grid on; xlabel('t'); ylabel('V');
set(gca,'FontSize',15);

subplot 212
plot(t1,a*(t1>1),'g','LineWidth',2)
grid on; xlabel('t'); ylabel('I');
set(gca,'FontSize',15); 

%% odehh definition 
function dydt = odehh(t, y ,I)

dydt = zeros(2,1);

v = y(1);
n = y(2);

m_inf = @(V) (1./(1+exp((-20 - V)/15)));
n_inf = @(V) (1/(1+exp((-25 - V)/10)));

if t<3
    Ix = 0;
else
    Ix = I;
end


dydt(1) = Ix - 20 * m_inf(v) * (v-61) - 10 * n * (v+90) - 8 * (v+78);
dydt(2) = (n_inf(v) - n) / 0.15;

disp(t)
return
end

%% ode definition

function dydt = ode(t, y ,Vc, Vs)

dydt = zeros(3,1);

n = y(1);
m = y(2);
h = y(3);

ninf = @(V) (1/(1+exp(((-53) - V)/15)));
minf = @(V) (1/(1+exp(((-40) - V)/9)));
hinf = @(V) (1/(1+exp(((-62) - V)/(-7))));

tawm = @(V) (0.04)+ (0.46).*exp(-((-38)-V).^2./30*30);
tawn = @(V) (1.1) + (4.7).*exp(-((-79)-V).^2./50*50);
tawh = @(V) (1.2) + (7.4).*exp(-((-67)-V).^2./20*20);

if t<20
    Vx = Vc;
else
    Vx = Vs;
end

dydt(1) = (ninf(Vx) - n)./tawn(Vx);

dydt(2) = (minf(Vx) - m)./tawm(Vx);

dydt(3) = (hinf(Vx) - h)./tawh(Vx);

disp(t)
end

 