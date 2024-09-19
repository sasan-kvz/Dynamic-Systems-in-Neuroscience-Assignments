%%hw2_3
clc; close all ; clear

ts = [0 20];
y0 = zeros(2,1);
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);

I=1;
[t,y] = ode45(@(t,y) nap_k(t,y,I), ts, [1;1],opts);

plot(t,y(:,1),'r') 
grid on ; xlabel('t [s]') ; ylabel('V [mv]'); title('I=1pA')


I=2;
[t,y] = ode45(@(t,y) nap_k(t,y,I), ts, [1;1],opts);

figure
plot(t,y(:,1),'r')
grid on ; xlabel('t [s]') ; ylabel('V [mv]'); title('I=2pA')


I=3;
[t,y] = ode45(@(t,y) nap_k(t,y,I), ts, [1;1],opts);

figure
plot(t,y(:,1),'r')
grid on ; xlabel('t [s]') ; ylabel('V [mv]'); title('I=3pA')


I=4;
[t,y] = ode45(@(t,y) nap_k(t,y,I), ts, [1;1],opts);

figure
plot(t,y(:,1),'r')
grid on ; xlabel('t [s]') ; ylabel('V [mv]'); title('I=4pA')


I=5;
[t,y] = ode45(@(t,y) nap_k(t,y,I), ts, [1;1],opts);

figure
plot(t,y(:,1),'r')
grid on ; xlabel('t [s]') ; ylabel('V [mv]'); title('I=5pA')

%% nap_k definition
function dydt = nap_k(t, y ,I)

dydt = zeros(2,1);
v = y(1);
n = y(2);

EL =-80; 
v_hn=-25;
%EL=-78;
%v_hn=-45;

mf = @(V) 1./(1+exp(((-20)-V)./15));
nf = @(V) 1./(1+exp((v_hn-V)./5));

if t<5
    I0 = 0;
else
    I0 = I;
end

dydt(1) = I0 - 20 * mf(v) *(v-60) - 10 * n* (v+90) - 8 * (v-EL);
dydt(2) = (nf(v)- n);

disp(t);
return
