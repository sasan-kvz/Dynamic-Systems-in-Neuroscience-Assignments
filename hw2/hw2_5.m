hw2_5
%% HW2_5a
close all; clear ; clc

c=281*10^-12; gl=30*10^-9; EL=-70.6*0.001; VTS=-50.4*0.001;delta=2*0.001; 
tw=144*0.001; a=4*10^-9; b=0.0805*10^-9; Vtr=20*0.001; h=0.01;

tspan=0:0.015:1000;
I=1e-6*stepfun(tspan,200);
v= zeros(1,length(I));
v2= zeros(1,length(I));


    for n = 1:length(I)-1
        
    v(n+1) = v(n) + h * ((1/c)*(- gl*(v(n)-EL) + gl*delta*exp((v(n)-VTS)/delta)+ I(n)- v2(n)));
    v2(n+1)=v2(n)+h*(1/tw*(a*(v(n)-EL)-v2(n)));
    if v(n+1)>Vtr
       v(n+1)=EL;
       v2(n+1)=v2(n)+b;
   end

    end

subplot 211
plot(tspan,v)
ylabel('V'); xlabel('t')

subplot 212
plot(tspan,I)
ylabel('I'); xlabel('t')

%% hw2_5b
close all; clear ; clc

c=281*10^-12; gl=30*10^-9; EL=-70.6*0.001; VTS=-50.4*0.001;delta=2*0.001; 
tw=144*0.001; a=4*10^-9; b=0.0805*10^-9; Vtr=20*0.001; h=0.01;

tspan=0:0.015:1000;
I= (2.5 * 1e-9)*square((2*pi/ 10)*tspan);
I(I<0) = 0;
I(tspan<200)=0;

v= zeros(1,length(I)); 
v2= zeros(1,length(I));

    for n = 1:length(I)-1
        
    v(n+1) = v(n) + h * ((1/c)*(- gl*(v(n)-EL) + gl*delta*exp((v(n)-VTS)/delta)+ I(n)- v2(n)));
    v2(n+1)=v2(n)+h*(1/tw*(a*(v(n)-EL)-v2(n)));
   if v(n+1)>Vtr
       v(n+1)=EL;
       v2(n+1)=v2(n)+b;
   end

    end
    

subplot 211
plot(tspan,v)
ylabel('V'); xlabel('t')

subplot 212
plot(tspan,I)
ylabel('I'); xlabel('t')

%% % HW2_5c
close all; clear ; clc

c=281*10^-12; gl=30*10^-9; EL=-70.6*0.001; VTS=-50.4*0.001;delta=2*0.001; 
tw=144*0.001; a=4*10^-9; b=0.0805*10^-9; Vtr=20*0.001; h=0.01;

tspan=-200:0.015:1200;
I=0.5*10^-9*(heaviside(tspan)-heaviside(tspan-200))+0.8*10^-9*(heaviside(tspan-500)-heaviside(tspan-1000));

v= zeros(1,length(I));
v2= zeros(1,length(I));

    for n = 1:length(I)-1
        
    v(n+1) = v(n) + h * ((1/c)*(- gl*(v(n)-EL) + gl*delta*exp((v(n)-VTS)/delta)+ I(n)- v2(n)));
    v2(n+1)=v2(n)+h*(1/tw*(a*(v(n)-EL)-v2(n)));
   if v(n+1)>Vtr
       v(n+1)=EL;
       v2(n+1)=v2(n)+b;
   end

    end
    
subplot 211
plot(tspan,v)
ylabel('V'); xlabel('t')

subplot 212
plot(tspan,I)
 ylabel('I'); xlabel('t')
