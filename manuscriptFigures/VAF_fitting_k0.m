%% VAF model fitting depending on k0

% Model 1: Binomial
% Model 2: Beta-Binomial

% Binomial: f_BI(0,n,VAF) = (1-VAF)^n;
% Beta-Binomial: f_BB(0|n,VAF) = B(a,n-b)/B(a,b) 
% = B(VAF*n+1,n(2-VAF)+1)/B(n*VAF+1,n(1-VAF)+1)
% where, VAF: % of variant; n: number of sequencing

% k0 is how many detect is needed to confirm existance
k0 = 1;

% True Data
TRUE_WASHU = [0.9942642161	0.9908672697	0.9912493949	0.9851209685	0.9632734802	0.8586419378];
TRUE_UW = [0.9900086345	0.9846174932	0.9765317505	0.9060026874	0.7090390686	0.5106605148];
TRUE_NYGC = [0.9921672629	0.9870988118	0.9879028015	0.9696736911	0.8148399579	0.6129644245];
TRUE_BCM = [0.9936012088	0.9888197647	0.9890341184	0.9769205108	0.8997997352	0.738161729];
TRUE_SIM1 = [0.988065869	0.9820181868	0.9832012587	0.971953927	0.8963035878	0.6970009713];
TRUE_SIM2 = [0.991350068	0.986524366	0.988865736	0.984967641	0.985132888	0.962594099];
% TRUE_DATA = [TRUE_WASHU;TRUE_BCM;TRUE_NYGC;TRUE_UW;TRUE_SIM1;TRUE_SIM2];
TRUE_DATA = [TRUE_WASHU;TRUE_BCM;TRUE_NYGC;TRUE_UW];

% Recall Rate: RR
VAF = [0.1	0.05	0.02	0.01	0.005	0.0025];
n_WASHU = 508; n_UW = 226; n_NYGC = 271; n_BCM = 464; n_SIM1 = 485; n_SIM2 = 4000;
% n_total = [n_WASHU, n_BCM, n_NYGC, n_UW, n_SIM1, n_SIM2];
n_total = [n_WASHU, n_BCM, n_NYGC, n_UW];

RR_WASHU_BI = NaN(length(n_total),length(VAF));
RR_WASHU_BB = NaN(length(n_total),length(VAF));

for j = 1:length(n_total)
    p_nondetect = 0;
    for k = 0:(k0-1)
        p_nondetect = p_nondetect + ...
            nchoosek(n_total(j),k) .* VAF.^k .* (1-VAF).^(n_total(j)-k);
    end
    RR_WASHU_BI(j,:) = 1-p_nondetect;
    for i = 1:length(VAF)
        p_nondetect = 0;
        for k = 0:(k0-1)
            p_nondetect = p_nondetect + ...
                nchoosek(n_total(j),k) .* betabinomial(k,n_total(j),VAF(i));
        end

        RR_WASHU_BB(j,i) = 1-p_nondetect;
    end
end

% Recall Rate with sequencing error
pE = 0.001; % Illumina Nova seq error rate = 0.1%
RR_WASHU_BI_ER = NaN(length(n_total),length(VAF));
RR_WASHU_BB_ER = NaN(length(n_total),length(VAF));
for j = 1:length(n_total)
    n_temp = n_total(j);

    % Wrong, old one
    % p_nondetect = 0;
    % for k = 0:(k0-1)
    %     error_noVAF_rate = 0;
    %     for l = k:n_temp
    %         error_noVAF_rate = error_noVAF_rate + (1-pE)^(n_temp-l) * (pE)^l * ...
    %             nchoosek(l,k) * (2/3)^(l-k) * (1/3)^k;
    %     end
    %     p_nondetect = p_nondetect + ...
    %             nchoosek(n_total(j),k) .* VAF.^k .* (1-VAF).^(n_total(j)-k) * ...
    %             error_noVAF_rate;
    % end
    % RR_WASHU_BI_ER(j,:) = 1-p_nondetect;

    p_nondetect = 0;
    for k = 0:k0-1 %k0:n_temp
        p_nondetect = p_nondetect + ...
            nchoosek(n_temp,k)*(VAF*(1-pE)+(1-VAF)*pE/3).^k.*...
            (1-(VAF*(1-pE)+(1-VAF)*pE/3)).^(n_temp-k);
    end
    RR_WASHU_BI_ER(j,:) = 1-p_nondetect;

    for i = 1:length(VAF)

        % Wrong, old one
        % p_nondetect = 0;
        % for k = 0:(k0-1)
        %     error_noVAF_rate = 0;
        %     for l = k:n_temp
        %         error_noVAF_rate = error_noVAF_rate + (1-pE)^(n_temp-l) * (pE)^l * ...
        %             nchoosek(l,k) * (2/3)^(l-k) * (1/3)^k;
        %     end
        %     p_nondetect = p_nondetect + ...
        %             nchoosek(n_total(j),k) .* betabinomial(k,n_total(j),VAF(i)) * ...
        %             error_noVAF_rate;
        % end
        % RR_WASHU_BB_ER(j,i) = 1-p_nondetect;

        p_nondetect = 0;
        for k = 0:k0-1 %k0:n_temp
        p_nondetect = p_nondetect + nchoosek(n_temp,k) * ...
            betabinomial(k,n_total(j),(VAF(i)*(1-pE)+(1-VAF(i))*pE/3));
        end
        RR_WASHU_BB_ER(j,i) = 1-p_nondetect;
        
    end
end

% Recall Rate when sequencing error is Poisson dist.
pE = 0.001; % Illumina Nova seq error rate = 0.1%
RR_WASHU_BI_ER_PO = NaN(length(n_total),length(VAF));
RR_WASHU_BB_ER_PO = NaN(length(n_total),length(VAF));
for j = 1:length(n_total)
    n_temp = n_total(j);

    % % Wrong, old one
    % x = 0:n_temp;
    % y = poisspdf(x,n_temp*pE);
    % p_nondetect = 0;
    % for k = 0:(k0-1)
    %     error_noVAF_rate = 0;
    %     for l = k:n_temp
    %         error_noVAF_rate = error_noVAF_rate + y(l+1) * ...
    %             nchoosek(l,k) * (2/3)^(l-k) * (1/3)^k;
    %     end
    %     p_nondetect = p_nondetect + ...
    %             nchoosek(n_total(j),k) .* VAF.^k .* (1-VAF).^(n_total(j)-k) * ...
    %             error_noVAF_rate;
    % end
    % RR_WASHU_BI_ER_PO(j,:) = 1-p_nondetect;

    p_nondetect = 0;
    for k = 0:(k0-1)
        poi = 0;
        for l = 0:(k0-k-1)
            lambda = (1-VAF)*pE/3;
            poi = poi + exp(-lambda).*lambda.^l./factorial(l);
        end
        q = VAF.*(1-pE);
        p_nondetect = p_nondetect + nchoosek(n_total(j),k) .* q.^k .* (1-q).^(n_temp-k) .* poi;
    end
    RR_WASHU_BI_ER_PO(j,:) = 1-p_nondetect;


    for i = 1:length(VAF)

        % % Wrong, old one
        % p_nondetect = 0;
        % for k = 0:(k0-1)
        %     error_noVAF_rate = 0;
        %     for l = k:n_temp
        %         error_noVAF_rate = error_noVAF_rate + y(l+1) * ...
        %             nchoosek(l,k) * (2/3)^(l-k) * (1/3)^k;
        %     end
        %     p_nondetect = p_nondetect + ...
        %             nchoosek(n_total(j),k) .* betabinomial(k,n_total(j),VAF(i)) * ...
        %             error_noVAF_rate;
        % end

        p_nondetect = 0;
        for k = 0:k0-1 
            poi = 0;
            for l = 0:(k0-k-1)
                lambda = (1-VAF(i))*pE/3;
                poi = poi + exp(-lambda).*lambda.^l./factorial(l);
            end
            p_nondetect = p_nondetect + nchoosek(n_temp,k) * ...
            betabinomial(k,n_total(j),(VAF(i)*(1-pE))) .* poi;
        end

        RR_WASHU_BB_ER_PO(j,i) = 1-p_nondetect;
    end
end

% Figure 1; Binomial model fitting
figure, hold on; 
b = bar(TRUE_DATA','FaceColor','flat');
xticks([1:length(VAF)]);
xticklabels({'0.1','0.05','0.02','0.01','0.005','0.0025'});
xticklabels({'10','5','2','1','0.5','0.25'});
% for i = 1:4
%     b(i).CData = [1 0 0];
% end
% b(5).CData = [0.8 0.8 0.8];
% b(6).CData = [0.8 0.8 0.8];
x = 1:length(VAF);
for i = 1:length(n_total)
    plot(x+(i-2.5)*0.18,RR_WASHU_BI_ER(i,:),'rs','MarkerSize',7);
end
legend({'Binomial w/SQ-ER'})
for i = 1:length(n_total)
    plot(x+(i-2.5)*0.18,RR_WASHU_BI(i,:),'ks','MarkerSize',7);
end
% legend({'WASHU','BCM','NYGC','UW','SIM485','SIM4000',...
%     'Binomial w/SQ-ER','','','','','',...
%     'Binomial','','',''},'Location','southwest')
legend({'WASHU','BCM','NYGC','UW',...
    'Binomial w/SQ-ER','','','','','',...
    'Binomial','','',''},'Location','southwest')
ylim([0 1]);
title('Variant Detect Rate for different VAF(%)')
set(gca, 'XDir','reverse')
hold off;

% Figure 2; comparing between models (VAF=10%,0.25%)
figure, yl = [0 1];

subplot(1,2,1), hold on;
bar(TRUE_DATA(:,1));
plot(RR_WASHU_BI(:,1),'o-');
plot(RR_WASHU_BB(:,1),'o-');
plot(RR_WASHU_BI_ER(:,1),'o-');
plot(RR_WASHU_BB_ER(:,1),'o-');
plot(RR_WASHU_BI_ER_PO(:,1),'o-');
plot(RR_WASHU_BB_ER_PO(:,1),'o-');
xticks([1:length(n_total)]);
% xticklabels({'WASHU','BCM','NYGC','UW','SIM485','SIM4000'})
xticklabels({'WASHU','BCM','NYGC','UW'})
legend({'True set','Binomial','Beta-Binomial',...
    'Bin w/SQ-ER','Beta-Bin w/SQ-ER','Bin w/Poi-SQ-ER',...
    'Beta-Bin w/Poi-SQ-ER'},'Location','southeast');
title('Recall rate comparison between models (VAF=10%)')
ylim(yl);
hold off;

subplot(1,2,2), hold on;
bar(TRUE_DATA(:,end));
plot(RR_WASHU_BI(:,end),'o-');
plot(RR_WASHU_BB(:,end),'o-');
plot(RR_WASHU_BI_ER(:,end),'o-');
plot(RR_WASHU_BB_ER(:,end),'o-');
plot(RR_WASHU_BI_ER_PO(:,end),'o-');
plot(RR_WASHU_BB_ER_PO(:,end),'o-');
xticks([1:length(n_total)]);
% xticklabels({'WASHU','BCM','NYGC','UW','SIM485','SIM4000'})
xticklabels({'WASHU','BCM','NYGC','UW'})
legend({'True set','Binomial','Beta-Binomial',...
    'Bin w/SQ-ER','Beta-Bin w/SQ-ER','Bin w/Poi-SQ-ER',...
    'Beta-Bin w/Poi-SQ-ER'},'Location','southeast');
title('Recall rate comparison between models (VAF=0.25%)')
ylim(yl);
hold off;


% figure, hold on;
% bar(TRUE_DATA');
% plot(RR_WASHU_BI','o-');
% ylim([0.3 1])
% title('Binomial estimate');
% 
% figure, hold on;
% bar(TRUE_DATA');
% plot(RR_WASHU_BB','o-');
% ylim([0.3 1])
% title('Beta-Binomial estimate');



%% Example data setup (replace these with your actual data):
% Suppose you have 5 sets of data, each with 6 points
N = 4;  % number of data sets
M = 6;  % number of data points in each set

%% Preallocate arrays for SSE and RMSE
SSE_model1 = zeros(N,1);
SSE_model2 = zeros(N,1);
RMSE_model1 = zeros(N,1);
RMSE_model2 = zeros(N,1);
RMSE_model3 = zeros(N,1);
RMSE_model4 = zeros(N,1);
RMSE_model5 = zeros(N,1);
RMSE_model6 = zeros(N,1);

%% Compute SSE and RMSE for each data set
for i = 1:N
    % Extract real data, model1 predictions, model2 predictions
    y_true = TRUE_DATA(i, :);
    y_pred1 = RR_WASHU_BI(i, :);
    y_pred2 = RR_WASHU_BB(i, :);
    y_pred3 = RR_WASHU_BI_ER(i, :);
    y_pred4 = RR_WASHU_BB_ER(i, :);
    y_pred5 = RR_WASHU_BI_ER_PO(i, :);
    y_pred6 = RR_WASHU_BB_ER_PO(i, :);
    
    % Sum of Squared Errors for model1 and model2
    y_mean = mean(y_true);
    SStot = sum((y_true-y_mean).^2);

    SSE_model1(i) = sum((y_true - y_pred1).^2);
    SSE_model2(i) = sum((y_true - y_pred2).^2);
    SSE_model3(i) = sum((y_true - y_pred3).^2);
    SSE_model4(i) = sum((y_true - y_pred4).^2);
    SSE_model5(i) = sum((y_true - y_pred5).^2);
    SSE_model6(i) = sum((y_true - y_pred6).^2);

    R2_model3(i) = 1-(SSE_model3(i)/SStot);
    
    % Root Mean Square Error for model1 and model2
    RMSE_model1(i) = sqrt(mean((y_true - y_pred1).^2));
    RMSE_model2(i) = sqrt(mean((y_true - y_pred2).^2));
    RMSE_model3(i) = sqrt(mean((y_true - y_pred3).^2));
    RMSE_model4(i) = sqrt(mean((y_true - y_pred4).^2));
    RMSE_model5(i) = sqrt(mean((y_true - y_pred5).^2));
    RMSE_model6(i) = sqrt(mean((y_true - y_pred6).^2));
end

%% Aggregate metrics across all data sets
% You can either look at sums/means across all data sets or look at them individually.

total_SSE_model1 = sum(SSE_model1);
total_SSE_model2 = sum(SSE_model2);

mean_SSE_model1  = mean(SSE_model1);
mean_SSE_model2  = mean(SSE_model2);

total_RMSE_model1 = sum(RMSE_model1);
total_RMSE_model2 = sum(RMSE_model2);

mean_RMSE_model1  = mean(RMSE_model1);
mean_RMSE_model2  = mean(RMSE_model2);
mean_RMSE_model3  = mean(RMSE_model3);
mean_RMSE_model4  = mean(RMSE_model4);
mean_RMSE_model5  = mean(RMSE_model5);
mean_RMSE_model6  = mean(RMSE_model6);

mean_RMSE = [mean_RMSE_model1,mean_RMSE_model2,mean_RMSE_model3,...
    mean_RMSE_model4,mean_RMSE_model5,mean_RMSE_model6];


%% Figure 4; Heatmap of VAF(%) x N(#Seq)
figure, hold on;
bar(mean_RMSE,'FaceColor',[0.7 0.7 0.7]);
xticks([1:6]);
xticklabels({'Binomial','Beta-Binomial','Bin w/SQ-ER','Beta-Bin w/SQ-ER',...
    'Bin w/Poi-SQ-ER','Beta-Bin w/Poi-SQ-ER'});
% yticks([]);
ylabel('RMSE')
title('RMSE for six models (lower the better)')
R2 = mean(R2_model3);

hold off;


%% Figure 3; RMSE
VAF = [0.1 0.05 0.025 0.01 0.005 0.0025 0.001];
n_seq = [30 60 100 200 300 400 500 600 1000];
pE = 0.001; % Assuming sequencing error rate is 0.1%
RR = NaN(length(n_seq), length(VAF));

for j = 1:length(n_seq)
    n_temp = n_seq(j);
    p_detect = 0;
    for k = k0:n_temp
        p_detect = p_detect + ...
            nchoosek(n_temp,k)*(VAF*(1-pE)+(1-VAF)*pE/3).^k.*...
            (1-(VAF*(1-pE)+(1-VAF)*pE/3)).^(n_temp-k);
    end
    RR(j,:) = p_detect;
end
RR = RR*100; % convert unit to percent(%)

figure,
h = heatmap(RR');
% h.CellLabelFormat = '%.2e';
title(['Variant detection rate (%; required supporting read = ' num2str(k0), ')'] ) 
h.YLabel = 'VAF (%)';
h.YDisplayLabels = {'10', '5', '2.5', '1', '0.5','0.25', '0.1'};
h.XLabel = 'Coverage (X)';
h.XDisplayLabels = ({'30', '60', '100', '200','300' ,'400', '500', '600', '1000'});
J = customcolormap_preset('red-white-blue');
h.Colormap = J;

%% Display results
fprintf('=== Model Comparison Results ===\n');

fprintf('Model 1: Total SSE  = %.4f, Mean SSE  = %.4f\n', total_SSE_model1, mean_SSE_model1);
fprintf('Model 1: Total RMSE = %.4f, Mean RMSE = %.4f\n', total_RMSE_model1, mean_RMSE_model1);

fprintf('Model 2: Total SSE  = %.4f, Mean SSE  = %.4f\n', total_SSE_model2, mean_SSE_model2);
fprintf('Model 2: Total RMSE = %.4f, Mean RMSE = %.4f\n', total_RMSE_model2, mean_RMSE_model2);

if mean_SSE_model1 < mean_SSE_model2
    disp('=> Model 1 has a lower SSE and may be the better fit.');
else
    disp('=> Model 2 has a lower SSE and may be the better fit.');
end

if mean_RMSE_model1 < mean_RMSE_model2
    disp('=> Model 1 has a lower RMSE and may be the better fit.');
else
    disp('=> Model 2 has a lower RMSE and may be the better fit.');
end





function pdf = betabinomial(k,n,p)
    pdf = beta(k+n*p+1,-k+n*(2-p)+1)/beta(n*p+1,n*(1-p)+1);
end

function pdf = poisson(k,lambda)
    pdf = lambda^k*exp(-lambda)/factorial(k);
end

