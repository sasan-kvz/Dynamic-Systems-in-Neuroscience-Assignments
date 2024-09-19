%% hw4_3
clc; close all; clear; hold on;

%b=0.01;I=0;a=0.1;c=0.02;
%b=0.01;I=0;a=-0.1;c=0.02;
 b=0.01;I=5;a=-0.1;c=0.02;

s = @(t,x) [-x(1)*(x(1)-a)*(x(1)-1)-x(2)+I;b*(x(1)-c*x(2))];
vectorfield(s,-15:.2:15,-1:.2:1);
[t,xs] = ode45(s,[0 100],[1 .1]);
plot(xs(:,1),xs(:,2));
x=-5:.01:5;
plot(x,x./c,'r',x,I-x.*(x-a).*(x-1),'b')
hold off
axis([-3 3 -1 1])
xlabel('u(t)','FontSize',10)
ylabel('v(t)','FontSize',10)
hold off
%% vectorfield definition

function vectorfield(deqns,xval,yval,t)
if nargin==3;
    t=0;
end
m=length(xval);
n=length(yval);
x1=zeros(n,m);
y1=zeros(n,m);
for a=1:m
  for b=1:n
    pts = feval(deqns,t,[xval(a);yval(b)]);
    x1(b,a) = pts(1);
    y1(b,a) = pts(2);
  end
end
arrow=sqrt(x1.^2+y1.^2);
quiver(xval,yval,x1./arrow,y1./arrow,.5,'r');
axis tight;


