clc
clear

%interpolation
y= [0.046595548	0.039171938	0.036391507	0.035328126	0.034911479	0.034746169	0.034678138	0.034651474	0.034640264	0.034634682
];
x=1:10;

%cubic spline
%number of conditions conditions(4n)
total=(length(x)-1)*4;
A=zeros(total); % (x-xi)
B=zeros(total,1); %soln we need [a0_1,a1_1,a2_1,a3_1,,,a2_5,a3_5]
C=zeros(total,1);
n=length(x)-1;
% a0_i,a1_i,a2_i,a3_i, i=1...5

%Interpolation condition, p(x)=f(x), p1(x0) = f(x0); n+1 conditions
a=0;
for i =1:n
    xval=x(i+1);
    for j = 1:4
        a=a+1;
        A(i,a)=(x(i+1)-xval)^(j-1);
    end
    C(i)=y(i+1);
end
for j=1:4
    A(n+1,j)=(x(1)-x(2))^(j-1);
end
C(n+1)=y(1);
total_cond = n+1;

%Continuity condition, p_i(x_i)=p_i+1(x_i), n-1 conditions
a=0;
for i =1:n-1
    total_cond=total_cond+1;
    for j=1:4
        a=a+1;
        A(total_cond,a)= (x(i+1)-x(i+1))^(j-1);
        A(total_cond,a+4) = -(x(i+1)-x(i+2))^(j-1);
    end
    C(total_cond) = 0;
end

%Smoothness condition, p'_i(x_i) = p'_i+1(x_i) ,n-1 conditions
% P'_i(x) = a1+ 2*a2(x-xi)+3*a3(x-xi)^2
a=0;
for i =1:n-1
    total_cond=total_cond+1;
    for j=1:4
        a=a+1;
        if j==1
            A(total_cond,a)=0;
            A(total_cond,a+4) = 0;
        elseif j==2
            A(total_cond,a)=1;
            A(total_cond,a+4) = -1;
        else
            A(total_cond,a)= (x(i+1)-x(i+1))^(j-1);
            A(total_cond,a+4) = -(j-1)*(x(i+1)-x(i+2))^(j-2);
        end
    end
    C(total_cond) = 0;
end

%2nd derivative condition, P''_i(xi) = p''_i+1(xi), n-1 conditions
% p''_i(x) = 2*a2+6*a2(x-xi)
a=0;
for i = 1:n-1
    total_cond=total_cond+1;
    for j = 1:4
        a=a+1;
        if j==3
            A(total_cond,a)=2;
            A(total_cond,a+4)=-2;
        elseif j==4
            A(total_cond,a+4)=-6*(x(i+1)-x(i+2))^(j-3);
        end
    end
    C(total_cond)=0;
end

%natural spline condition, 2conditions
total_cond=total_cond+1;
C(total_cond)=0;
A(total_cond,3)=2;
A(total_cond,4)=6*(x(1)-x(2));

total_cond=total_cond+1;
C(total_cond)=0;
A(total_cond,end-1)=2;

B = inv(A)*C;
count=1;
for i = 1:length(x)-1
    for j=1:4
        soln(i,j)=B(count);
        count=count+1;
    end
end

%testing
testx = 1:0.01:10;
for i = 1:length(testx)
    xi = floor(testx(i)+1);
    xx = [1; (testx(i)-xi); (testx(i)-xi)^2; (testx(i)-xi)^3];
    if xi<11
        cspline(i) = soln(xi-1,:)*xx;
    else
        cspline(i) = soln(9,:)*xx;
    end
end
plot(testx,cspline)
hold on
scatter(x,y)
title('Cubic Spline')


%%
%Newton Polynomial
newton = zeros(n+1,n+1);
newton(1,:) = y;
for i =1:n
    for j = 1:n-i+1
        newton(i,j+1)
        newton(i,j)
        x(j+i)
        x(j)
        newton(i+1,j) = (newton(i,j+1)-newton(i,j))/(x(j+i)-x(j));
    end
end

newton_coefs= newton(:,1);
newton_x = x(1:9)';

for i = 1:length(testx)
    newton_soln = 0;
    xvar = testx(i)-newton_x;
    for j = 1:length(x)-1
        xcoef(j,1) = prod(xvar(1:j));
    end
    newton_soln = xcoef'*newton_coefs(2:end)+newton_coefs(1);
    newtonpoly(i) = newton_soln;
end
plot(testx,newtonpoly)
hold on
scatter(x,y)
title('Newton Polynomial')

%%
%Lagrange Polynomial


L=ones(length(x),length(testx));
for k = 1:length(x)
    for i = 1:length(x)
        if k~=i
            L(k,:)=L(k,:).*(testx-x(i))/(x(k)-x(i));
        end
    end
end

lagpoly = 0;
for i =1:length(x)
    lagpoly = lagpoly+y(i)*L(i,:);
end
plot(testx,lagpoly)
hold on
scatter(x,y)
title('Lagrange Polynomial')
%%
%Univariate Polynomial Regression
for i = 1:10
    UPRx(i,:) = [1,x(i),x(i)^2,x(i)^3,x(i)^4];
    %UPRx(i,:) = [1,x(i),x(i)^2,x(i)^3,x(i)^4,x(i)^5,x(i)^6];
end

beta = inv(UPRx'*UPRx)*UPRx'*y';
for i = 1:length(testx)
    xi = testx(i);
    UPRy(i) = [1,xi,xi^2,xi^3,xi^4]*beta;
    %UPRy(i) = [1,xi,xi^2,xi^3,xi^4,xi^5,xi^6]*beta;
end
plot(testx,UPRy)
hold on
scatter(x,y)
title('Univariate Polynomial Regression')


%%
% x = 1.5 @ index 51;
realsoln = 0.041638248064378;
model = [lagpoly(51),newtonpoly(51),cspline(51),UPRy(51)];
percenterror = abs(model-realsoln)/realsoln*100
plot(testx(1:300),cspline(1:300))
hold on
plot(testx(1:300),lagpoly(1:300))
plot(testx(1:300), newtonpoly(1:300))
plot(testx(1:300), UPRy(1:300))
scatter(1.5, realsoln)
legend('Cubic Spline','Lagrange Polynomial', 'Newton Polynomial','Univariate Polynomial Regression','Closed Form Solution')
