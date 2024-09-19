%%hw2_2
%% hw2_2a
clc; close all ; clear

t=0:0.01:250;
I=1.502*stepfun(t,100);
z= zeros(1,length(I));

n=1;
    while n <= length(I)-1
        
    z(n+1)=z(n)+(0.01)*((1/4)*(10*I(n)-(z(n)+75)));
    
   if z(n+1)>=-60
       z(n+1:n+1+100)=-85 ;
       n=n+100;
   else
       n=n+1;
   end
   
    
    end
    
plot(t,z)
title('Membrane Voltage'); ylabel('V'); xlabel('t')
%% hw2_2b
clc; close all ; clear

I = 1:1:60;
F = zeros(1, length(I));
for i = 1:length(I)-1
    [~, ~, f]= getf(I(i));
    F(i) = f;
end
[t, x1,f]= getf(1);

plot(I,F);
title('F-I Curve'); ylabel('F'); xlabel('I')

%% getf definition
function [t, z, f] = getf(scale)
    
    t=0:0.01:100;
    I=scale*stepfun(t,5);
   
    z= zeros(1,length(I));
    spikes = zeros(1,length(I));
    
    n=1;
    while n < length(I) -2
       z(n+1)=z(n)+(0.01)*((1/4)*(10*I(n)-(z(n)+75)));
       
       if z(n+1)>=-60
           z(n+1:n+1+100)=-85;
           n=n+100;
           spikes(n+1) = 1; 
       end
       n=n+1;
    end
    
    ix= find(spikes, 3);
    if (length(ix)>2)
        f = 1 / (t(ix(3)) - t(ix(2)));
    else
        f=0;
    end
end