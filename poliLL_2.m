function f=poliLL_2(x)

global rf1 rf2 eqc ls mu0 C months den mutemp prob1 ls_new nu p pm mum m1 m2 mu nu_form % temp1 temp2 temp3 %nm nm0 prob1
global final_ls dres d_sq % NEEDED FOR PLOTTING
global w p mu den sum_temp
global D n C

%Initial values
mu1 = exp(x(1));
mu2 = exp(x(2));
mu3 = exp(x(3));
nu0 = x(4);

%define nu form
if nu_form == 0
    nu = nu0*ones(size(rf1,1),size(rf1,2));
else
    if nu_form == 1
        nu1 = x(5);
        nu = nu0 + nu1 * rf1;
    else
        if nu_form == 2
            nu2 = x(5);
            nu = nu0 + nu2 * rf2;
        else
            nu1 = x(5);
            nu2 = x(6);
            nu3 = (7);
            nu = nu0 + nu1 * rf1 + nu2 * rf2 + nu3 * eqc;
        end
    end
end
p = exp(-nu)./(1+exp(-nu));

%calculation of the overall effect of rain and quakes per day and location
logmu0 = mu1 * rf1 + mu2 * rf2 + mu3 * eqc;%matrix (x,t)
mutemp = (1-p).*exp(logmu0); %matrix (x,t)
mudash = mutemp';
lsdash = ls';
mu0 = sum(lsdash)./sum(mudash); 
%calculation of mu
mu = (mu0'*ones(1,length(mutemp))).*mutemp;

%REDISTRIBUTING LANDSLIDES
%min/max month number
%months
m1 = min(months);
m2 = max(months);
%initialize denominator (sum per month)
den = zeros(size(ls,1),size(ls,2));
% den = [months; den];

%calculation of the denominator sum((1-p)*mu)
for j = m1:m2
    mon_temp = months==j;
    p_temp = p(:,mon_temp); %select prob of zip model for month j
    mu_temp = mu(:,mon_temp); %select mu for month j
    sum_temp = sum((1-p_temp).*mu_temp,2); %calculation of the denominator for month j. result is a vector 1xnumber of municipalities
    sum_temp = repmat(sum_temp,1,size(p_temp,2)); %replicate vector to actual size of month
    den(:,mon_temp) = sum_temp; %save result into denominator matrix
end
% den(1,:)=[]; %get rid of month count
den = den + (den == 0)*1e-12;    % AVOID DIVIDE BY ZERO ERROR BELOW
w = ((1-p).*mu)./den; %calculation of the probabilities

%adding number of months to landslide matrix
lsm = ls; 
n1_r =  double.empty(1,0); %initialize a vector for non-accurate landslides on day 1 to be reallocated
ls_new =  double.empty(1,0);
for j = m1:m2
    mon_temp = months==j;
    ls_temp = lsm(:,mon_temp);%select landslides for month j
    daysofmonth = ones(size(ls_temp));
    n1_nact = ls_temp(:,1)*(1-prob1); %proportion of non accurate fraction to move
    n1_r = [n1_r diag(n1_nact)*daysofmonth]; %save non accurate once in a matrix;
    ls_temp(:,1) = ls_temp(:,1)*prob1; % replace the day-1 ls with the proportion of accurate first day landslides
    lsm(:,mon_temp) = ls_temp; %save the landslides into the original landslide matrix
end
% lsm(1,:)=[];% lsm is now a matrix with accurate ls day1 and other days
n1_nac = n1_r.*w; % multiply non-accurate ls to w probabilities
n = lsm + n1_nac; % new matrix of ls
final_ls = n;

%LIKELIHOOD
B = n; %create an identity matrix for n>0
B(B>0)=1;
D = B;  % for residuals
A = 1-B; %create an identity matrix for n=0
p_n0 = p + (1-p).*exp(-mu);  %probability when 0 landslides
logp_n0 = log(p_n0);
p_n1 = (1-p).*(mu.^n).*exp(-mu);%./factorial(n);
logp_n1 = log(1-p) + n.*log(mu) -mu;
 
A = A.* logp_n0;
B = B.* logp_n1;
C = A + B;
n_temp = n + (n == 0);
D = D.*(-n + n.*log(n_temp) - log(gamma(n+1)));

%dres = -2*sum(C - B.*log(gamma(n+1)));

d_sq = 2*sum(D - C);
%dres = sqrt(d_sq).*(((sum(D) - sum(C))>0) - ((sum(D) - sum(C))<0));

%d_sq = 2*(D - C);
%dres = sqrt(d_sq).*(((n - mu)>0) - ((n - mu)<0));
d_sq = sum(d_sq);
dres = sqrt(d_sq).*(((sum(n) - sum(mu))>0) - ((sum(n) - sum(mu))<0));

lambda = sum(sum(C));

f = -lambda;