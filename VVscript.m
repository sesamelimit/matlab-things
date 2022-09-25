clc
clear all;
close all;

syms x1(t) x2(t) x3(t)

Y=[x1;x2;x3];

A=[-13.5,-12.5,-43;
    5.2,4.2,17.8;
    2.1,2.1,6.4];
S=[-5,-1,7;
    2,1,-3;
    1,0,-1;];
A_0=inv(S)*A*S;
B=[1;0;0];
K=[12.5;-4.5;-2.5];
B_0=inv(S)*B;
K_0=inv(S)*K;
m_1=-1.11;
m_2=5.30;

G=[-0.14;0;0];

C=Y(0)==[-5.36;4.62;2.52];

odes=diff(Y)==A*Y + B_0*m_1 + K_0*f(t);

[w1(t),w2(t),w3(t)]=dsolve(odes,C);

figure
plot3(w1(t),w2(t),w3(t));

function y=f(x)
y=1+2*sin(x+pi/3)+5*sin(2*x);
end


