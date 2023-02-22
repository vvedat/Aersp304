%{
Part.d
1. reduce the second order to a 1st order ODE with 4 equations
using state variables.
2. define a function with time(t) and y(4 columns) representing 4 state 
variables. x, x_dot, y, y_dot
3. Call a ODE solver using function with given T and initial conditions
4. plot x and y
5. define the DCM that convert from b1,b2 to n1, n2
6. multiply x and y of solution by the DCM, to get X and Y
7. 
%}


load("EM_L2-304P1.mat");

tspan = [0 .2]; 
[t,y] = ode45(@ODEsystem,tspan,x0(1:4));
plot(t,y,'-o')
figure()
plot(y(:,1),y(:,3))

function dydt = ODEsystem(t,y)
load("EM_L2-304P1.mat");
dydt = zeros(4,1);


% y1 -> x
% y2 -> x_dot
% y3 -> y
% y4 -> y_dot
%  p1 = sqrt(y(1).^2+2.*MU1.*y(1)+MU1.^2+y(3).^2);
%  p2 = sqrt((1-MU1-y(1)).^2+y(3).^2)
% 
% Ux=((-1+MU1).*(y(1)+MU1)/p1.^3)-(MU1*(y(1)-1+MU1)./p2.^3)+y(1);
% Uy=-1*((1-MU1).*y(3)/p1.^3)-(MU1*y(3)/p2.^3)+y(3);



Ux=((-1+MU1).*(y(1)+MU1)/(y(1).^2+2.*MU1.*y(1)+MU1.^2+y(3).^2).^(3/2))-(MU1*(y(1)-1+MU1)./((1-MU1-y(1)).^2+y(3).^2).^(3/2))+y(1);
Uy=-1*((1-MU1).*y(3)/(y(1).^2+2.*MU1.*y(1)+MU1.^2+y(3).^2).^(3/2))-(MU1*y(3)/((1-MU1-y(1)).^2+y(3).^2).^(3/2))+y(3)

dydt(1) = y(2);
dydt(2) = Ux +2*y(4);
dydt(3) =  y(4);
dydt(4) =  Uy-2*y(2);

end