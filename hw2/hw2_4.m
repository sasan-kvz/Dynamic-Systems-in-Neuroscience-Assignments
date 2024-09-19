%%hw2_4
clc; close all ; clear

y1 = zeros(2,1) ; 
tspan = [0, 20];  
opt = odeset('RelTol',1e-4,'AbsTol',1e-4); 

[~,y] = ode45(@(t,y) nat_m_inf(t,y,10), tspan, y1,opt); 
y1 = y(end, :);
[t, y] = ode45(@(t,y) nat_m_inf(t,y,10), tspan, y1, opt);


plot(t,y(:,1)); 
title('v memberane'); xlabel('t'); ylabel('V')
%% nat_m_inf definition
function dydt = nat_m_inf(t, y ,I)

dydt = zeros(2,1);
v = y(1); 
h= y(2); 

if t<5
    Ix = 0;
else
    Ix = I;
end

hinf= 1/(1+exp((-70-v)/-5.8));
tawh=3/exp((-40-v)/33);
dydt(1) =   Ix - (v+70) - 10*h*(v-60); 
dydt(2) = (hinf- h)/tawh; 

disp(t)
return

