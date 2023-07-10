%% ME303 TERM PROJECT
%% QUESTION 2  
%% EULER'S METHOD

%% Start by closing open figures and clearing the memory and the termianl 
% outputs. 
clc
clear all
close all

%% Define initial and final time. Also define the step size and create
% a time matrix with that data. 
ti=0;
tf=17.06521656015796;
h=tf/6000;
T=ti:h:tf;

%% Define the sub-groups with the same length of T.
X1=zeros(size(T));
X2=zeros(size(T));
X3=zeros(size(T));
X4=zeros(size(T));


%% Specify their inital values.
X1(1)=0;
X2(1)=0.994;
X3(1)=-2.0015851063790825;
X4(1)=0;

%% The implemattion of the rk4 method.
tic
for J=1:length(T)-1
    
    x1_1 = X1_P(X1(J) , X2(J) , X3(J) , X4(J) );
    x2_1 = X2_P(X1(J) , X2(J) , X3(J) , X4(J) );
    x3_1 = X3_P(X1(J) , X2(J) , X3(J) , X4(J) );
    x4_1 = X4_P(X1(J) , X2(J) , X3(J) , X4(J) );

    
    x1_2 = X1_P(X1(J)+h/2*x1_1 , X2(J)+h/2*x2_1 , X3(J)+h/2*x3_1 , X4(J)+h/2*x4_1 );
    x2_2 = X2_P(X1(J)+h/2*x1_1 , X2(J)+h/2*x2_1 , X3(J)+h/2*x3_1 , X4(J)+h/2*x4_1 );
    x3_2 = X3_P(X1(J)+h/2*x1_1 , X2(J)+h/2*x2_1 , X3(J)+h/2*x3_1 , X4(J)+h/2*x4_1 );
    x4_2 = X4_P(X1(J)+h/2*x1_1 , X2(J)+h/2*x2_1 , X3(J)+h/2*x3_1 , X4(J)+h/2*x4_1 );

    
    x1_3 = X1_P(X1(J)+h/2*x1_2 , X2(J)+h/2*x2_2 , X3(J)+h/2*x3_2 , X4(J)+h/2*x4_2 );
    x2_3 = X2_P(X1(J)+h/2*x1_2 , X2(J)+h/2*x2_2 , X3(J)+h/2*x3_2 , X4(J)+h/2*x4_2 );
    x3_3 = X3_P(X1(J)+h/2*x1_2 , X2(J)+h/2*x2_2 , X3(J)+h/2*x3_2 , X4(J)+h/2*x4_2 );
    x4_3 = X4_P(X1(J)+h/2*x1_2 , X2(J)+h/2*x2_2 , X3(J)+h/2*x3_2 , X4(J)+h/2*x4_2 );

    
    x1_4 = X1_P(X1(J)+h*x1_3 , X2(J)+h*x2_3 , X3(J)+h*x3_3 , X4(J)+h*x4_3 );
    x2_4 = X2_P(X1(J)+h*x1_3 , X2(J)+h*x2_3 , X3(J)+h*x3_3 , X4(J)+h*x4_3 );
    x3_4 = X3_P(X1(J)+h*x1_3 , X2(J)+h*x2_3 , X3(J)+h*x3_3 , X4(J)+h*x4_3 );
    x4_4 = X4_P(X1(J)+h*x1_3 , X2(J)+h*x2_3 , X3(J)+h*x3_3 , X4(J)+h*x4_3 );

    
    X1(J+1) = X1(J) + (x1_1+2*x1_2+2*x1_3+x1_4)*h/6;
    X2(J+1) = X2(J) + (x2_1+2*x2_2+2*x2_3+x2_4)*h/6;
    X3(J+1) = X3(J) + (x3_1+2*x3_2+2*x3_3+x3_4)*h/6;
    X4(J+1) = X4(J) + (x4_1+2*x4_2+2*x4_3+x4_4)*h/6;
    
end
toc

%% Plot y1(t) and y2(t) 
plot(X2,X4)
hold on

%% Name the variables in the graph, name the x and y axis and name the 
% graph itself. Also turn on the grids for people to better see. 
title('y1(t) vs y2(t)')
xlabel('y1(t)')
ylabel('y2(t)')
grid on

%% Define the derivatives as a funcitons to pass on ode45.
function rvalx1 = X1_P(x1,x2,x3,x4)
mass1=0.012277471;
mass2=1-mass1;
d1=(((x2+mass1)^2+x4^2))^(3/2);
d2=(((x2-mass2)^2+x4^2))^(3/2);
rvalx1=x2 + 2*x3 - mass2*(x2+mass1)/d1 - mass1*(x2-mass2)/d2;
end

function rvalx2 = X2_P(x1,x2,x3,x4)
rvalx2=x1;
end

function rvalx3 = X3_P(x1,x2,x3,x4)
mass1=0.012277471;
mass2=1-mass1;
d1=(((x2+mass1)^2+x4^2))^(3/2);
d2=(((x2-mass2)^2+x4^2))^(3/2);
rvalx3=x4 - 2*x1 - mass2*x4/d1 - mass1*x4/d2;
end

function rvalx4=X4_P(x1,x2,x3,x4)
rvalx4=x3;
end