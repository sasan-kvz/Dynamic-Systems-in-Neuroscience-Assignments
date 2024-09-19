%% hw2 
%% hw2_1
clc; close all ; clear;

t_s = [0 200];
opt = odeset('RelTol',1e-6,'AbsTol',1e-6);
%%hw2_1_a
[t,y] = ode45(@(t,y) fhn(t,y,0.1,0.01,0.02,5), t_s, [0;0],opt);

V = y(:,1);
W = y(:,2);

subplot 211
plot(t,V,'b','LineWidth',2)
xlabel('t'); ylabel('V'); hold on ; title('a=0.1')

subplot 212
plot(t,W,'r','LineWidth',2) 
grid on; xlabel('t'); ylabel('W'); 
%%hw2_1_b

[t,y] = ode45(@(t,y) fhn(t,y,-0.1,0.01,0.02,5), t_s, [0;0],opt);

V = y(:,1);
W = y(:,2);

figure
subplot 211
plot(t,V,'b','LineWidth',2)
xlabel('t'); ylabel('V'); hold on ; title('a=-0.1')
subplot 212
plot(t,W,'r','LineWidth',2)
grid on; xlabel('t'); ylabel('W'); 

%% FitzHugh_Nagumo definition
function dydt = fhn(t, y ,a, b, c, I)

dydt = zeros(2,1);

V = y(1);
W = y(2);

dydt (1) = V.*(a-V).*(V-1)-W+I;
dydt (2) = b*V-c*W;

return





