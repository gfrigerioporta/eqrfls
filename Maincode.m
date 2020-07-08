% landslide.m
% ===========
global rf eqm eqd ls_loc rf1 rf2 eqc ls mu0 p mu eqtime dom nm lsz prob1 prob2 months nu nu_form  %%%%% LAST THREE ARE ADDED %%%%
global final_ls dres 

% MODEL (0, 1 or 2)
mult_models = [2]; %[1 2 3];
nu_forms = [3]; %[0 1 2 3];   % 0 = nu0 only, 1 = nu0 + nu1*rf1, 2 = nu0 + nu2*rf2, 3 = nu0 + nu1*rf1 + nu2*rf2 + nu3*eq 

% FLAGS
transform_rf1 = 0;
transform_rf2 = 0;

% PARAMETERS 
k = 150;
om = 0.98;
eq_alpha = 1.5;
eq_beta = 1;
max_eq_dist = 400;

% CONSTANTS
convno = 1;
noruns = 3;
min_ls = 10;  % only include municipalities with this many ls

% LOAD DATA
load('eqm.txt') % matrix of earthquake magnitude per day and location
load('eqd.txt') % matrix of earthquake distances per day and location
load('eqtime.txt') % earthquake times
load('rf.txt')% matrix of precipitation per day and location
load('ls.txt'); % matrix of landslides per day and location
load('dom.txt'); %day/month matrix

%restructuring the datasets
eqtime = eqtime-366+1;
ls_temp = ls;
months = ls_temp(1,2:end);
ls_temp(1,:)=[];
ls_temp(:,1)=[];
ls_mun = [ls(2:end,1) sum(ls_temp')'];
ls_loc = ls_mun(:,1);
ls_mun(:,1) = [];
rf_temp = rf;
rf_temp(:,1) = [];

e_loc = eqm(:,1);    % replacement code to match earthquake locations to landslide locations
ls80 = sum(ls(2:end,2:368),2);
ls_loc = ls_loc(ls_mun-ls80>=min_ls,:);   % get labels for locations with minimum number of ls
temp = 0*e_loc;
for i = 1:length(temp)
    for j = 1:length(ls_loc)
        if e_loc(i) == ls_loc(j)
            temp(i) = 1;             % temp = 1 if that location is included, 0 otherwise
        end
    end
end
eqm = eqm(temp==1,:);   % subset earthquake locations to match ls 
eqd = eqd(temp==1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rf = rf(ls_mun-ls80>=min_ls,:);% assumes ls and rf are sorted the same way by location label
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ls = ls_temp(ls_mun-ls80>=min_ls,:);   % use ls_temp, which is only ls numbers
ls = [ls_loc ls];            % to allow location to be stripped out later
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Building the earthquake component
nodays = size(ls_temp,2);

eqd(:,1)=[]; %get rid of them in the eqd matrix
eqm(:,1)=[]; %get rid of them in the eqm matrix
rf(:,1)=[];
ls(:,1)=[];

eqd = [eqtime'; eqd];
eqm = [eqtime'; eqm];

for i = 1:size(eqd,2)
    if mean(eqd(:,i))>max_eq_dist
        eqd(:,i) = zeros(size(eqd,1),1);
        eqm(:,i) = zeros(size(eqm,1),1);
    end
end

toDelete       = (sum(eqd(2:end,:), 1) == 0);
eqd(:,toDelete) = [];
eqm(:,toDelete) = [];
eqtime(toDelete) = [];
eqd(1,:)=[];
eqm(1,:)=[];

eqc = eq_maker2(eqd,eqm,eqtime,eq_alpha,eq_beta,nodays);

%Rainfall component
days = (0:k-1);
omdays = sort(om.^days);

rf(rf>900)=0;
rf(:,13149:end)=rf(:,13149:end)*10;


rf1 = zeros(size(rf,1),size(rf,2));
for i=1:size(rf,1)
    for j=2:size(rf,2)
         rf1(i,j) = (rf(i,j-1)+rf(i,j))/2;
    end
end

%%%%TRANSFORMATION%%%%
if transform_rf1
    tmp1 = rf1(rf1>0);
    mtmp1 = mean(mean(tmp1));
    stmp1 = std2(tmp1);
    rf1a = 0.5*(1 + erf((rf1-mtmp1)./(stmp1*sqrt(2))));
    rf1b = rf1.*(rf1==0) + rf1a.*(rf1>0);
    rf1b(rf1b>0) = (rf1b(rf1b>0)-min(min(rf1a)))/(max(max(rf1a))-min(min(rf1a)));
    rf1=rf1b;
end

%calculation of longterm
rf2 = zeros(size(rf,1),size(rf,2));                                         % size(rf)
for i=1:size(rf,1)
     for j=(k+2):size(rf,2)
     rf2(i,j) = sum(omdays.*rf(i,j-(k+2)+1:j-2))/k;                           % ERROR FIXED
     end
end

if transform_rf2
    tmp2 = rf2(rf2>0);
    mtmp2 = mean(mean(tmp2));    %%%% FIXED FROM tmp1
    stmp2 = std2(tmp2);          %%%% FIXED FROM tmp1
    rf2a = 0.5*(1 + erf((rf2-mtmp2)./(stmp2*sqrt(2))));
    rf2b = rf2.*(rf2==0) + rf2a.*(rf2>0);
    rf2b(rf2b>0) = (rf2b(rf2b>0)-min(min(rf2a)))/(max(max(rf2a))-min(min(rf2a)));
    rf2=rf2b;
end


rf1(:,1:366) = [];
rf2(:,1:366) = [];
ls(:,1:366) = [];
dom(:,1:366) = [];
eqc = eqc(:,1:length(dom));
months(:,1:366) = [];
nodays = size(months,2);

eqc(1,:)=[];
dom(1,:)=[];
eqc = eqc/10;
dom2 = dom;
dom(dom==1)=0;
dom(dom~=0)=1;

%generate number of days per month
nom = zeros(420,1);
j=0;
for i=2:size(dom2,2)
    if dom2(i)==1
        j=j+1;
        nom(j)=dom2(i-1);
    end
end
nom(end) = dom2(end);

nm = zeros(nodays,420);
j=0;                                              
for i=1:size(nom)
    nm(j+1:j+nom(i),i) = ones(nom(i),1);
    j=j+nom(i);
end

% calculate probability of type 1 and 2 events
type1days = length(nom);
type2days = sum(nom) - type1days;
lst = sum(ls);
num1 = sum(lst(dom==0));
num2 = sum(lst(dom==1));
prob1 = (type1days/(type1days+type2days))/(num1/(num1+num2));
lsz = ls;
lsz(lsz~=0) = 1;

%Parameter estimation
for hidx = 1:length(mult_models)
    h = mult_models(hidx) %for each type of mu model
    for kidx = 1:length(nu_forms) %for each type of nu model
        k = nu_forms(kidx)
        mult_model = h;
        nu_form = k;   % 0 = nu0 only, 1 = nu0 + nu1*rf1, 2 = nu0 + nu2*rf2, 3 = nu0 + nu1*rf1 + nu2*rf2 + nu3*eq 
        results = zeros(4,noruns); %matrix to save results
        allx = 1:0;
        allL = 1:0;
        for i = 1:noruns
            i
            no_conv = 0;  % number of convergent fit
            while no_conv < convno
                mu1_0 = -log(rand(1));
                mu2_0 = -log(rand(1));
                mu3_0 = -log(rand(1));
                nu0_0 = 5*randn(1);

                if nu_form == 0
                    if mult_model == 2
                        [x,fval,exitflag] = fminsearch(@poliLL_2,[log(mu1_0),log(mu2_0),log(mu3_0),nu0_0]);
                    else
                        if mult_model == 1
                            [x,fval,exitflag] = fminsearch(@poliLL_1,[log(mu1_0),log(mu2_0),log(mu3_0),nu0_0]);
                           
                        else
                            if mult_model == 3
                                [x,fval,exitflag] = fminsearch(@poliLL_3,[log(mu1_0),log(mu2_0),log(mu3_0),nu0_0]);
                            end
                        end
                    end
                else
                    if nu_form == 1
                        nu1_0 = 5*randn(1);
                        if mult_model == 2
                            [x,fval,exitflag] = fminsearch(@poliLL_2,[log(mu1_0),log(mu2_0),log(mu3_0),nu0_0,nu1_0]);
                        else
                            if mult_model == 1
                                [x,fval,exitflag] = fminsearch(@poliLL_1,[log(mu1_0),log(mu2_0),log(mu3_0),nu0_0,nu1_0]);
                            else
                                if mult_model == 3
                                    [x,fval,exitflag] = fminsearch(@poliLL_3,[log(mu1_0),log(mu2_0),log(mu3_0),nu0_0,nu1_0]);
                                end
                            end
                        end
                    else
                        if nu_form == 2
                            nu2_0 = 5*randn(1);
                            if mult_model == 2
                                [x,fval,exitflag] = fminsearch(@poliLL_2,[log(mu1_0),log(mu2_0),log(mu3_0),nu0_0,nu2_0]);
                                if fval < 39971
                                    save 'mu2nu2.mat';
                                    no_conv = convno+1;
                                end
                            else
                                if mult_model == 1
                                    [x,fval,exitflag] = fminsearch(@poliLL_1,[log(mu1_0),log(mu2_0),log(mu3_0),nu0_0,nu2_0]);
                                else
                                    if mult_model == 3
                                        [x,fval,exitflag] = fminsearch(@poliLL_3,[log(mu1_0),log(mu2_0),log(mu3_0),nu0_0,nu2_0]);
                                    end
                                end
                            end
                        else
                            nu1_0 = 5*randn(1);
                            nu2_0 = 5*randn(1);
                            nu3_0 = 5*randn(1);
                            if mult_model == 2
                                [x,fval,exitflag] = fminsearch(@poliLL_2MB,[log(mu1_0),log(mu2_0),log(mu3_0),nu0_0,nu1_0,nu2_0,nu3_0]);
                            else
                                if mult_model == 1
                                    [x,fval,exitflag] = fminsearch(@poliLL_1,[log(mu1_0),log(mu2_0),log(mu3_0),nu0_0,nu1_0,nu2_0,nu3_0]);
                                else
                                    if mult_model == 3
                                        [x,fval,exitflag] = fminsearch(@poliLL_3,[log(mu1_0),log(mu2_0),log(mu3_0),nu0_0,nu1_0,nu2_0,nu3_0]);
                                    end
                                end
       
                            end
                        end
                    end
                end

            if exitflag > 0 % convergent fit
                allx = [allx; x];
                allL = [allL; -fval];
                no_conv = no_conv + 1
                if no_conv == 1 % first convergent fit
                    best_fval = fval + 1;
                end
                if fval < best_fval % if new solution better than previous best, update
                    best_fval = fval
                    best_x = x;  % log(parameters) for best fit
                    best_ls = final_ls;    %% for plotting
                end
            end
        end

        if nu_form == 0
            results(1+h,i) = exp(best_x(1));
            results(2+h,i) = exp(best_x(2));
            results(3+h,i) = exp(best_x(3));
            results(4+h,i) = best_x(4);
            results(5+h,i) = -best_fval;
        else
            if nu_form == 1
                results(1+h,i) = exp(best_x(1));
                results(2+h,i) = exp(best_x(2));
                results(3+h,i) = exp(best_x(3));
                results(4+h,i) = best_x(4);
                results(5+h,i) = -exp(best_x(5));
                results(6+h,i) = -best_fval;
            else
                if nu_form == 2
                    results(1+h,i) = exp(best_x(1));
                    results(2+h,i) = exp(best_x(2));
                    results(3+h,i) = exp(best_x(3));
                    results(4+h,i) = best_x(4);
                    results(5+h,i) = -exp(best_x(5));
                    results(6+h,i) = -best_fval;
                else
                    results(1+h,i) = exp(best_x(1));
                    results(2+h,i) = exp(best_x(2));
                    results(3+h,i) = exp(best_x(3));
                    results(4+h,i) = best_x(4);
                    results(5+h,i) = -exp(best_x(5));
                    results(6+h,i) = -exp(best_x(6));
                    results(7+h,i) = -exp(best_x(7));
                    results(8+h,i) = -best_fval;
                end
            end
        end
        end
    end
end

%%%%% PLOTTING ROUTINE %%%%

%calculation of expected landslides
calc_mu_nu;
final_ls_round = round(final_ls);
lsd = sum(final_ls_round); %sum of daily landslides
lsm = sum(final_ls_round'); %sum of landslides per municipality
lsm_nr = sum(final_ls');
dd = 1:nodays;

%plot observed landslides vs mu0
a=figure(1);set(a, 'Visible', 'on')
 plot(lsm_nr,mu0,'o')
 xlabel('Cumulative number of observed landslides')
 ylabel('Parameter mu0 for each location')
ddy = 1981 + (dd-1)/365.25;
%plot observed vs expected landslides expobs_nu3_mu2
b=figure(10);set(b, 'Visible', 'on')
semilogy(ddy,lsd,'.k',ddy,mu,'k'); xlabel('Days'); ylabel('No. landslides'); legend('Observed','Expected');
ylim([min(mu) max(max(lsd),max(mu))+1]);xlim([min(ddy) max(ddy)]);


meanrf1 = mean(rf1);
meanrf2 = mean(rf2);
dd1 = ddy(meanrf1 > 0);
dd2 = ddy(meanrf2 > 0);
meanrf1 = meanrf1(meanrf1>0);
meanrf2 = meanrf2(meanrf2>0);

ax = gca;
ax.YRuler.Exponent = 1;
yt = get(gca,'YTick');
set(gca,'YTickLabel', sprintf('%.0f',yt))


c=figure('units','normalized','outerposition',[0 0 1 1]);set(c, 'Visible', 'on')
subplot(3,1,3); semilogy(ddy,mean(eqc),'.k'); xlabel('Year'); ylabel('Earthquake component'); 
subplot(3,1,2); semilogy(ddy,lsd,'.k',ddy,mu,'k'); xlabel('Year'); ylabel('No. landslides'); legend('Observed','Expected'); %ylim([min(mu) 15])
hAx = subplot(3,1,1);[hAx,H1,H2]= plotyy(dd1,meanrf1,dd2,meanrf2); xlabel('Year');ylabel(hAx(1),'Short term rainfall') % left y-axis 
set(H2,'LineStyle', 'none','marker','.','color', [0.5 0.5 0.5]);set(hAx(2),'YLim',[-2.5 2.5],'YTick',[-2.5:2.5:2.5],'YColor', [0.5 0.5 0.5]);ylabel(hAx(2),'Long term rainfall');
set(H1,'linewidth',1,'color', [0 0 0]);set(hAx(1),'YLim',[0 max(meanrf1)],'YTick',[0:25:50],'YColor', [0 0 0])


c1=figure('units','normalized','outerposition',[0 0 1 1]);set(c1, 'Visible', 'on')
subplot(3,1,3); semilogy(ddy,mean(eqc),'.k'); xlabel('Year'); ylabel('Earthquake component'); xlim([2016 max(ddy)]);
subplot(3,1,2); semilogy(ddy,lsd,'.k',ddy,mu,'k'); xlabel('Year'); ylabel('No. landslides'); legend('Observed','Expected');xlim([2016 max(ddy)]);
%subplot(3,1,1); plot(dd1,meanrf1/max(meanrf1),'.k',dd2,meanrf2/max(meanrf2),'-k'); xlabel('Year'); ylabel('Rainfall components'); legend('RF_1','RF_2'); 
hAx = subplot(3,1,1);[hAx,H1,H2]= plotyy(dd1,meanrf1,dd2,meanrf2); xlabel('Year');ylabel(hAx(1),'Short term rainfall') % left y-axis 
set(H2,'LineStyle', 'none','marker','.','color', [0.5 0.5 0.5]);set(hAx(2),'YLim',[-2.5 2.5],'YTick',[-2.5:2.5:2.5],'YColor', [0.5 0.5 0.5]);ylabel(hAx(2),'Long term rainfall');xlim([2016 max(ddy)]);
set(H1,'linewidth',1,'color', [0 0 0]);set(hAx(1),'YLim',[0 50],'YTick',[0:25:50],'YColor', [0 0 0])

%logical vectors to separate components in days with or without landslides
lsr = lsd~=0;
lsf = lsd==0;
%splitting components in landslide days/no landslide days
rf1_r = m2rv(rf1(:,lsr));
rf2_r = m2rv(rf2(:,lsr));
eqc_r = m2rv(eqc(:,lsr));

rf1_f = m2rv(rf1(:,lsf));
rf2_f = m2rv(rf2(:,lsf));
eqc_f = m2rv(eqc(:,lsf));

%Plot of components:
e=figure(100);set(e, 'Visible', 'on')
subplot(1,3,1);semilogy(rf1_f,rf2_f,'s','color', [0.5 0.5 0.5], 'MarkerSize',5); xlabel('Short term rainfall'); ylabel('Long term rainfall'); ylim([min(rf2_f) max(rf2_f)]); % title('2 days vs 14 days rain components');
hold on
subplot(1,3,1);semilogy(rf1_r,rf2_r,'.k'); xlabel('Short term rainfall'); ylabel('Long term rainfall');  ylim([min(rf2_r) max(rf2_r)]);;% title('2 days vs 14 days rain components');
hold off
subplot(1,3,2);semilogy(rf1_f,eqc_f,'s','color', [0.5 0.5 0.5], 'MarkerSize',5); xlabel('Short term rainfall'); ylabel('Earthquake'); ylim([min(eqc_f) max(eqc_f)]);% title('2 days rain components vs earthquake component');
hold on
subplot(1,3,2);semilogy(rf1_r,eqc_r,'.k'); xlabel('Short term rainfall'); ylabel('Earthquake'); ylim([min(eqc_r) max(eqc_r)]);% title('2 days rain components vs earthquake component');
hold off
subplot(1,3,3);semilogy(rf2_f,eqc_f,'s','color', [0.5 0.5 0.5], 'MarkerSize',5); xlabel('Long term rainfall'); ylabel('Earthquake'); ylim([min(eqc_f) max(eqc_f)]);% title('14 days rain components vs earthquake component');
hold on
subplot(1,3,3);semilogy(rf2_r,eqc_r,'.k'); xlabel('Long term rainfall'); ylabel('Earthquake'); ylim([min(eqc_r) max(eqc_r)]);% title('14 days rain components vs earthquake component');
hold off

%components plot but dividing each component by its minimum
rf1_ra = rf1_r/min(min(rf1_r(rf1_r~=0)));
rf2_ra = rf2_r/min(min(rf2_r(rf2_r~=0)));
eqc_ra = eqc_r/min(min(eqc_r(eqc_r~=0)));

rf1_fa = rf1_f/min(min(rf1_f(rf1_f~=0)));
rf2_fa = rf2_f/min(min(rf2_f(rf2_f~=0)));
eqc_fa = eqc_f/min(min(eqc_f(eqc_f~=0)));

f=figure(104);set(f, 'Visible', 'on')
subplot(1,3,1);semilogy(rf1_fa,rf2_fa,'s','color', [0.5 0.5 0.5], 'MarkerSize',5); xlabel('Short term Rainfall'); ylabel('Long term Rainfall'); % title('2 days vs 14 days rain components');
hold on
subplot(1,3,1);semilogy(rf1_ra,rf2_ra,'.k'); xlabel('Short term rainfall'); ylabel('Long term Rainfall'); % title('2 days vs 14 days rain components');
hold off
subplot(1,3,2);semilogy(rf1_fa,eqc_fa,'s','color', [0.5 0.5 0.5], 'MarkerSize',5); xlabel('Short term Rainfall'); ylabel('Earthquake'); % title('2 days rain components vs earthquake component');
hold on
subplot(1,3,2);semilogy(rf1_ra,eqc_ra,'.k'); xlabel('Short term rainfall'); ylabel('Earthquake'); % title('2 days rain components vs earthquake component');
hold off
subplot(1,3,3);semilogy(rf2_fa,eqc_fa,'s','color', [0.5 0.5 0.5], 'MarkerSize',5); xlabel('Long term Rainfall'); ylabel('Earthquake'); % title('14 days rain components vs earthquake component');
hold on
subplot(1,3,3);semilogy(rf2_ra,eqc_ra,'.k'); xlabel('Long term Rainfall'); ylabel('Earthquake'); % title('14 days rain components vs earthquake component');
hold off

%components plot but multiplying each component by the parameter
rf1_rm = mu1*rf1(:,lsr);
rf2_rm = mu2*rf2(:,lsr);
eqc_rm = mu3*eqc(:,lsr);

rf1_fm = mu1*rf1(:,lsf);
rf2_fm = mu2*rf2(:,lsf);
eqc_fm = mu3*eqc(:,lsf);

h=figure(102);set(h, 'Visible', 'on')
subplot(1,3,1);semilogy(rf1_fm,rf2_fm,'s','color', [0.5 0.5 0.5], 'MarkerSize',5); xlabel('Short term Rainfall'); ylabel('Long term Rainfall'); ylim([min(min(min(rf2_fm)),min(min(rf2_rm))) max(max(max(rf2_fm)),max(max(rf2_rm)))]);% title('2 days vs 14 days rain components');
hold on
subplot(1,3,1);semilogy(rf1_rm,rf2_rm,'.k'); xlabel('Short term Rainfall'); ylabel('Long term Rainfall'); % title('2 days vs 14 days rain components');
hold off
subplot(1,3,2);semilogy(rf1_fm,eqc_fm,'s','color', [0.5 0.5 0.5], 'MarkerSize',5); xlabel('Short term Rainfall'); ylabel('Earthquake'); ylim([min(min(min(eqc_fm)),min(min(eqc_rm))) max(max(max(eqc_fm)),max(max(eqc_rm)))]);% title('2 days rain components vs earthquake component');
hold on
subplot(1,3,2);semilogy(rf1_rm,eqc_rm,'.k'); xlabel('Short term Rainfall'); ylabel('Earthquake'); % title('2 days rain components vs earthquake component');
hold off
subplot(1,3,3);semilogy(rf2_fm,eqc_fm,'s','color', [0.5 0.5 0.5], 'MarkerSize',5); xlabel('Long term Rainfall'); ylabel('Earthquake'); ylim([min(min(min(eqc_fm)),min(min(eqc_rm))) max(max(max(eqc_fm)),max(max(eqc_rm)))]);% title('14 days rain components vs earthquake component');
hold on
subplot(1,3,3);semilogy(rf2_rm,eqc_rm,'.k'); xlabel('Long term Rainfall'); ylabel('Earthquake'); % title('14 days rain components vs earthquake component');
hold off
