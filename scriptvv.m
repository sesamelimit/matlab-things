%первый набросок который рисует одну петлю

clc
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

odes=diff(Y)==A_0*Y + B_0*m_1 + K_0*f(t);

[w1(t),w2(t),w3(t)]=dsolve(odes,C);

%не всегда правильно решает почему-то
t1=vpasolve(w1==-28.57,t,[0 2]);

figure
fplot3(w1,w2,w3,[0 double(t1)],'Color','k','LineWidth',1)
hold on
fplot3(w1(t1),w2(t1),w3(t1),'-o','Color','k','MarkerSize',6,'MarkerFaceColor','#709494')
grid on
axis square

m_2=5.30;
C2=Y(t1)==[-28.57;4.46;0.77];
odes2=diff(Y)==A_0*Y+B_0*m_2+K_0*f(t);
[v1(t),v2(t),v3(t)]=dsolve(odes2,C2);
t2=vpasolve(v1==-5.36,t,[double(t1) +inf]);
fplot3(v1,v2,v3,[double(t1) double(t2)],'Color','k','LineWidth',1)
fplot3(v1(t2),v2(t2),v3(t2),'-o','Color','k','MarkerSize',6,'MarkerFaceColor','#709494')

[Ys,Zs] = meshgrid(-1:0.2:7,-3:0.2:3);
Xs2=-28.57+0*Ys+0*Zs;
surf(Xs2,Ys,Zs, 'LineStyle', ':')
Xs=-5.36+0*Ys+0*Zs;
surf(Xs,Ys,Zs, 'LineStyle', ':')
colormap(white)
title('Решение ОДУ diff(Y)==A_0*Y+B_0*m_i+K_0*f(t)')
xlabel('x') 
ylabel('y') 
zlabel('z')
text(w1(t1),w2(t1),w3(t1),'t=1.6')
text(v1(t2),v2(t2),v3(t2),'t2=2\pi')
hold off

function y=f(x)
y=1+2*sin(x+pi/3)+5*sin(2*x);
end
