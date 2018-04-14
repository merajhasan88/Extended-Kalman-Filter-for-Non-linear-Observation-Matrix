clc
clear all

s0=[5;5;0;0]; %initial state, declared far from true one
M0=100*eye(4,4); %Large MSE
w=[0.01*rand;0.01*rand]; % For Noisy observations 
for n=0:1:99
rx=10-0.2.*n; %ideal
ry=-5+0.2.*n; %ideal
delta=1;
u2=0.0001; %noise variance
A=[1,0,delta,0;0,1,0,delta;0,0,1,0;0,0,0,1]; %transition matrix
Q=[0,0,0,0;0,0,0,0;0,0,u2,0;0,0,0,u2];
x=[sqrt(rx^2+ry^2);atan(ry/rx)]+w; %state matrix
H=[rx/sqrt(rx^2+ry^2), ry/sqrt(rx^2+ry^2),0,0;-ry/(rx^2+ry^2),rx/(rx^2+ry^2),0,0]; %observation matrix
C=[0.1,0;0,0.01];

%Extended Kalman Filter
s1=A*s0;
M1=(A*M0*A.')+Q;
K=(M1*H.')*(inv(C+(H*M1*H.')));
s2(:,n+1)=s1+K*(x-hofsofn(s1));
M2=(eye(4,4)-K*H*M1);

%Update
s0=s2(:,n+1);
M0=M2;
end

%ideal vs estimate plot
n1=0:1:99;
rx1=10-0.2.*n1; %ideal
ry1=-5+0.2.*n1; %ideal
hold on
plot(rx1,ry1)
plot(s2(1:1,:),s2(2:2,:))
xlim([-15,15])
ylim([-10,20])
xlabel('rx (ideal vs estimate)')
ylabel('ry (ideal vs estimate)')
hold off
legend('ideal','estimate')

function hsn = hofsofn(x)
hsn=[sqrt(x(1,1)^2+x(2,1)^2);atan(x(2,1)/x(1,1))];
end
