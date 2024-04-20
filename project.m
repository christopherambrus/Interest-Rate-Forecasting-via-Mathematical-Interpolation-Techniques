clear
clc

%maturities
maturities = [0.5, 1,2,3,5,7,10,20,30];

%CMT, 3m,6m,1yr,2yr,3yr,5yr,7yr,10yr,20,30
cm=readmatrix('cmt_rates.xls');
cm = cm(5:end,3:12);
constantmat = rmmissing(cm);

%fminsearch to optimize RMSE
x0=[5,5,5];
x1=0;
x2=1;
A=[];
B=[];
Aeq=[];
Beq=[];
lb=[];
ub=[];
%soln = fmincon(@objfunc,x0,A,B,Aeq,Beq,lb,ub);
soln = [0.0181325124554064,0.523594883310337,0.000237301709060087];
RMSE = objfunc(soln);

%%
%Yield Curve 
%paths from 1-30yrs, 10,000 paths. 
alpha = soln(1);
beta = soln(2);
sigma = soln(3);
dt=1;
r0 =0.05;
ratepath = zeros(10000,30/dt);
ratepath(:,1)=r0;

%short rate simulation using monte carlo
for i =2:30/dt
    for j = 1:10000
        delta_term = (alpha-beta*ratepath(j,i-1))*dt + sigma*randn(1)*sqrt(dt);
        ratepath(j,i) = ratepath(j,i-1) + delta_term;
    end
end
avgrate = mean(ratepath);

%ytm, yield curve
for i =1:length(avgrate)
    ytm(i) = vdt(avgrate(i),alpha,beta,sigma,i/dt);
end




%%

function rtnext = rt(rate,alpha,theta,sigma,dt)
    term1 = alpha*(theta-rate)*dt;
    term2 = sigma*sqrt(dt)*randn(1);
    rtnext = rate+term1+term2;
end

function [ytm] = vdt(rate,alp,bet,sig,T)
    term1 = (  (sig^2)/(2*(bet^2))   -(alp/bet)   )*T;
    term2 = (alp/(bet^2)-(sig^2)/(bet^3))*(1-exp(-bet*T));
    term3 = (sig^2)/(4*(bet^3))*(1-exp(-2*bet*T));
    At = exp(term1+term2+term3);
    Bt = 1/bet*(1-exp(-bet*T));
    vasicekdt = At*exp(-Bt*rate);
    ytm = -log(At)/T+Bt/T*rate;
end

function [At,Bt] = shortrate(alp,bet,sig,T)
    term1 = (  (sig^2)/(2*(bet^2))   -(alp/bet)   )*T;
    term2 = (alp/(bet^2)-(sig^2)/(bet^3))*(1-exp(-bet*T));
    term3 = (sig^2)/(4*(bet^3))*(1-exp(-2*bet*T));
    At = exp(term1+term2+term3);
    Bt = 1/bet*(1-exp(-bet*T));
end

function par = parrate(time,DISCOUNT,i)
    dis= DISCOUNT;
    nomi= 100-100*dis;
    denom= sum(discount(i,1:2*time));
    par= 2*nomi/denom;
end

function rmse = objfunc(x)
    maturities2= [0.5, 1,2,3,5,7,10,20,30];

    cm=readmatrix('cmt_rates.xls');
    cm = cm(5:end,3:12);
    constantmat = rmmissing(cm);



    % obtain short rate r0 for all of them
    [At,Bt] = shortrate(x(1),x(2),x(3),0.25);

    % step3) short rate for 0.25yr at each date.
    for j=1:551
        term1 = constantmat(j,1)+log(At)/0.25;
        shortR(j,1) = term1*0.25/Bt/100;
    end

    % step4) discount rate from 0.5yr ~ 30yr for all date.
    for i = 1:551
        for j = 1:60
        [Att,Btt]=shortrate(x(1),x(2),x(3),j/2);
        discount(i,j) = Att*exp(-Btt*shortR(i));
        end
    end
    
    % step5) use Discount factors to obtain model par rates. 
    for i = 1:551
        for j = 1:length(maturities2)
            time = maturities2(j);
            dis = discount(i,time*2);
            nomi= 100-100*dis;
            denom= sum(discount(i,1:2*time));
            modelpar(i,j) = 2*nomi/denom;
        end
    end
    
    %error between actual CMT and model par rates
    errormatrix = constantmat(:,2:end)-modelpar;
    %RMSE
    sqrt(sum(errormatrix.^2,'all')/(numel(errormatrix)))
    rmse = sqrt(sum(sum(errormatrix.^2))/(numel(errormatrix)));

end
