# pop_code_aud_ctx
Analysis code from the manuscript "Hierarchical Differences in Population Coding Within Auditory Cortex." by Downer et al, 2017 Journal of Neurophysiology

 clear

% data must be pre-organized and loaded

% data are organized already, load them (the data comes from lines 0-350 of noise_correlation_population_model_ml_indmeans:

% There needs to be a number of matrices of size n*m where n is # of 
% neurons and m is # of distinct stimuli. In the following code, these
% variables store mean spike count ('active_mean_rate'), variance 
% ('active_var_rate'), standard deviation ('active_std_rate')

% There also need to be designations of the slope for each neuron, in other
% words does it encode AM with and increasing or decreasing function? This
% is a 1*n vector where n is # neurons, and each element is the slope of
% one neuron. Here, it is 'a_SU_slope_ml_indmeans'



% The model below has the following properties:

% correlations range from ~0:~0.3
corrs = 0.005:0.1:0.205;

% vary number of neurons in population
pop_size = [10 20 50 100 200];
pop_size = 100;

simulated_noise_correlation_ml = [];
active_performance_ml = [];
active_Iole_ml = [];

for pop = 1:length(pop_size)
    num_neurons = pop_size(pop)
    
    % 400 trials
    n_trials = 400;
    
    % type of fit: cubic ('poly3')
    
    raw_standard_fano = [];
    raw_test_fano = [];
    normd_fano = [];
    slope_vector = [];
    rtuning_vector = [];
    rtuning_mean = [];
    rtuning_median = [];
    
    % The simulation is run separately for each stimulus, in this case AM
    % depths:
    amd = [0 0.06 0.16 0.28 0.4 0.6 0.8 1];  
    
    % run through 20 iterations (20 distinct realizations of the modeled
    % populations) to average over
    for iteration = 1:20
        units_to_use = [];
        for i = 1:num_neurons % draw from random from the recorded neurons
            u = randperm(length(active_var_rate)); 
            units_to_use(i) = u(1);
        end
        
        
        max_means = max(active_mean_rate(units_to_use,:),[],2);
        best_stim = [];
        normd_rates = [];
        normd_means = [];
        mnormd = [];
        for i = 1:num_neurons
            for stim = 1:8
                normd_rates{i,stim} = trialwise_rates_active_ml{units_to_use(i),stim}(:)/max_means(i);
                normd_means(i,stim) = mean(normd_rates{i,stim});
            end
            if size(find(active_mean_rate(units_to_use(i),:)==max_means(i)))
                x = find(active_mean_rate(units_to_use(i),:)==max_means(i));
                best_stim(i) = x(1);
            else
                best_stim(i) = find(active_mean_rate(units_to_use(i),:)==max_means(i));
            end
            [r,mnormd(i),b] = regression(amd,normd_means(i,:));
        end
        
        slope_vector(iteration,:) = mnormd;
        rtuning_dist = corr(normd_means');
        rtuning_dist(logical(eye(size(rtuning_dist)))) = NaN;
        rtuning_vector(iteration,:) = rtuning_dist(:);
        rtuning_mean(iteration) = nanmean(rtuning_vector(iteration,:));
        rtuning_median(iteration) = nanmedian(rtuning_vector(iteration,:));
        
        ordervector = sortrows([mnormd' (1:length(mnormd))']);
        ordervector = ordervector(:,2);
        
        x = num_neurons;
        f100 = fit(slope_vector(iteration,ordervector)',normd_means(ordervector,8),'poly3');
        f080 = fit(slope_vector(iteration,ordervector)',normd_means(ordervector,7),'poly3');
        f060 = fit(slope_vector(iteration,ordervector)',normd_means(ordervector,6),'poly3');
        f040 = fit(slope_vector(iteration,ordervector)',normd_means(ordervector,5),'poly3');
        f028 = fit(slope_vector(iteration,ordervector)',normd_means(ordervector,4),'poly3');
        f016 = fit(slope_vector(iteration,ordervector)',normd_means(ordervector,3),'poly3');
        f006 = fit(slope_vector(iteration,ordervector)',normd_means(ordervector,2),'poly3');
        f000 = fit(slope_vector(iteration,ordervector)',normd_means(ordervector,1),'poly3');
        
        
        fit_100 = f100.p1*(slope_vector(iteration,ordervector)').^3 + f100.p2*(slope_vector(iteration,ordervector)').^2 + f100.p3.*(slope_vector(iteration,ordervector)') + f100.p4;
        fit_080 = f080.p1*(slope_vector(iteration,ordervector)').^3 + f080.p2*(slope_vector(iteration,ordervector)').^2 + f080.p3.*(slope_vector(iteration,ordervector)') + f080.p4;
        fit_060 = f060.p1*(slope_vector(iteration,ordervector)').^3 + f060.p2*(slope_vector(iteration,ordervector)').^2 + f060.p3.*(slope_vector(iteration,ordervector)') + f060.p4;
        fit_040 = f040.p1*(slope_vector(iteration,ordervector)').^3 + f040.p2*(slope_vector(iteration,ordervector)').^2 + f040.p3.*(slope_vector(iteration,ordervector)') + f040.p4;
        fit_028 = f028.p1*(slope_vector(iteration,ordervector)').^3 + f028.p2*(slope_vector(iteration,ordervector)').^2 + f028.p3.*(slope_vector(iteration,ordervector)') + f028.p4;
        fit_016 = f016.p1*(slope_vector(iteration,ordervector)').^3 + f016.p2*(slope_vector(iteration,ordervector)').^2 + f016.p3.*(slope_vector(iteration,ordervector)') + f016.p4;
        fit_006 = f006.p1*(slope_vector(iteration,ordervector)').^3 + f006.p2*(slope_vector(iteration,ordervector)').^2 + f006.p3.*(slope_vector(iteration,ordervector)') + f006.p4;
        fit_000 = f000.p1*(slope_vector(iteration,ordervector)').^3 + f000.p2*(slope_vector(iteration,ordervector)').^2 + f000.p3.*(slope_vector(iteration,ordervector)') + f000.p4;
        
        fits{1} = fit_000;fits{2} = fit_006;fits{3} = fit_016;fits{4} = fit_028;fits{5} = fit_040;fits{6} = fit_060;fits{7} = fit_080;fits{8} = fit_100;
        
        % check: does it look OK?
        
        figure
        plot(slope_vector(iteration,ordervector)',fit_100,'g','LineWidth',6)
        hold on
        plot(slope_vector(iteration,ordervector)',fit_000,'k','LineWidth',6)
        scatter(slope_vector(iteration,ordervector)',normd_means(ordervector,8),'g')
        scatter(slope_vector(iteration,ordervector)',normd_means(ordervector,1),'k')
        legend('100%','0%')
        set(gca,'FontSize',24)
        % If it looks promising, run a simulation:
        
        % Vary rnoise according to variable 'corrs'
        
        for stim = 1:8; % 8 distinct stimuli, one simulation for each
            for rep = 1:length(corrs)
                d_nc = corrs(rep); % desired noise correlation
                
                C = eye(num_neurons); % correlation matrix
                C(C==0) = d_nc;
                q_n=num_neurons;
                
                % these next 4 lines are to create 'u' from Shadlen 96 paper:
                q_u1= 1 / (d_nc*q_n^.5);
                q_u2= ( 2/q_n + d_nc - 2*d_nc/q_n - 2/q_n*( (1-d_nc)*(1-d_nc+d_nc*q_n) )^.5 )^.5;
                q_u3= (1 + ( (1-d_nc)*(1-d_nc+d_nc*q_n) )^.5 );
                q_u=q_u1*q_u2*q_u3; % q_u is 'u' from Shadlen 96
                % this next bit is 'v' from Shadlen 96
                q_v= 1/q_n^.5 * q_u2;
                
                q_mat=ones(q_n)*q_v;
                q_mat(logical(eye(size(q_mat))))=q_u; % q_mat is 'Q' from Shadlen 96
                
                % C = Q*Q';
                
                % make a matrix, where each row is a simulated neuron and each column is a trial
                resp_input=randn(q_n,n_trials); % same as 'z' from Shadlen 96
                % corr_rates will be 'n_trials' trials be z scored firing rates for 'num_neurons' neurons. It's a num_neurons x n_trials matrix
                corr_rates = q_mat*resp_input;
                
                % now convert the corr_rates to rates and variances based on the data:
                
                test_sigma_sq = [];test_mu = [];standard_sigma_sq = [];standard_mu = [];
                
                test_sigma_sq = active_var_rate(units_to_use,stim);
                standard_sigma_sq = active_var_rate(units_to_use,1);
                test_mu = active_mean_rate(units_to_use,stim);
                standard_mu = active_mean_rate(units_to_use,1);
                active_var_rate(active_var_rate==0) = 0.001;
                active_mean_rate(active_mean_rate==0) = 0.001;
                
                raw_standard_fano(:,stim,iteration) = active_var_rate(units_to_use,1)./active_mean_rate(units_to_use,1);
                raw_test_fano(:,stim,iteration) = active_var_rate(units_to_use,stim)./active_mean_rate(units_to_use,stim);
                
                % first, ensure that each distribution here TRULY has a mean of
                % 0:
                
                corr_rates = bsxfun(@plus,mean(corr_rates,2).*-1,corr_rates);
                
                % Preserve each unit's fano factor as well:
                
                ideal_test_mean = active_mean_rate(units_to_use,stim)./max_means;
                pool_test_mean = active_mean_rate(units_to_use,stim);
                ideal_standard_mean = active_mean_rate(units_to_use,1)./max_means;
                pool_standard_mean = active_mean_rate(units_to_use,1);
                
                ideal_test_fano = raw_test_fano(:,stim,iteration);
                ideal_standard_fano = raw_standard_fano(:,stim,iteration);
                % now just some algebra to get ideal variance:
                ideal_test_var = ideal_test_mean.*ideal_test_fano;
                pool_test_var = pool_test_mean.*ideal_test_fano;
                ideal_standard_var = ideal_standard_mean.*ideal_standard_fano;
                pool_standard_var = pool_standard_mean.*ideal_standard_fano;
                
                % multiply by variances:
                
                test_rates = bsxfun(@times,ideal_test_var.^0.5,corr_rates);
                pool_test_rates = bsxfun(@times,pool_test_var.^0.5,corr_rates);
                standard_rates = bsxfun(@times,ideal_standard_var.^0.5,corr_rates);
                pool_standard_rates = bsxfun(@times,pool_standard_var.^0.5,corr_rates);
                
                % then add normalized or raw means:
                
                test_rates = bsxfun(@plus,ideal_test_mean,test_rates);
                pool_test_rates = bsxfun(@plus,pool_test_mean,pool_test_rates);
                standard_rates = bsxfun(@plus,ideal_standard_mean,standard_rates);
                pool_standard_rates = bsxfun(@plus,pool_standard_mean,pool_standard_rates);
                                
                x = mean(pool_test_rates)-mean(mean(pool_test_rates));
                xx = mean(pool_standard_rates)-mean(mean(pool_test_rates));
                y = mean(pool_test_rates)-mean(mean(pool_standard_rates));
                yy = mean(pool_standard_rates)-mean(mean(pool_standard_rates));
                
                poolhit = sum(abs(x)<abs(y))/n_trials;
                poolfalse = sum(abs(xx)<abs(yy))/n_trials;
                
                % Now, do the above separately for pos & neg slope...
                
                % pos slope first:
                
                x = mean(pool_test_rates(a_SU_slope_ml_indmeans(units_to_use)>0,:))-mean(mean(pool_test_rates(a_SU_slope_ml_indmeans(units_to_use)>0,:)));
                xx = mean(pool_standard_rates(a_SU_slope_ml_indmeans(units_to_use)>0,:))-mean(mean(pool_test_rates(a_SU_slope_ml_indmeans(units_to_use)>0,:)));
                y = mean(pool_test_rates(a_SU_slope_ml_indmeans(units_to_use)>0,:))-mean(mean(pool_standard_rates(a_SU_slope_ml_indmeans(units_to_use)>0,:)));
                yy = mean(pool_standard_rates(a_SU_slope_ml_indmeans(units_to_use)>0,:))-mean(mean(pool_standard_rates(a_SU_slope_ml_indmeans(units_to_use)>0,:)));
                
                pospoolhit = sum(abs(x)<abs(y))/n_trials;
                pospoolhit = pospoolhit*(sum(a_SU_slope_ml_indmeans(units_to_use)>0))/num_neurons;
                pospoolfalse = sum(abs(xx)<abs(yy))/n_trials;
                pospoolfalse = pospoolfalse*(sum(a_SU_slope_ml_indmeans(units_to_use)>0))/num_neurons;
                
                x2 = mean(pool_test_rates(a_SU_slope_ml_indmeans(units_to_use)<0,:))-mean(mean(pool_test_rates(a_SU_slope_ml_indmeans(units_to_use)<0,:)));
                xx2 = mean(pool_standard_rates(a_SU_slope_ml_indmeans(units_to_use)<0,:))-mean(mean(pool_test_rates(a_SU_slope_ml_indmeans(units_to_use)<0,:)));
                y2 = mean(pool_test_rates(a_SU_slope_ml_indmeans(units_to_use)<0,:))-mean(mean(pool_standard_rates(a_SU_slope_ml_indmeans(units_to_use)<0,:)));
                yy2 = mean(pool_standard_rates(a_SU_slope_ml_indmeans(units_to_use)<0,:))-mean(mean(pool_standard_rates(a_SU_slope_ml_indmeans(units_to_use)<0,:)));
                
                negpoolhit = sum(abs(x2)<abs(y2))/n_trials;
                negpoolhit = negpoolhit*(sum(a_SU_slope_ml_indmeans(units_to_use)<0))/num_neurons;
                negpoolfalse = sum(abs(xx2)<abs(yy2))/n_trials;
                negpoolfalse = negpoolfalse*(sum(a_SU_slope_ml_indmeans(units_to_use)<0))/num_neurons;
                
                
                normd_fano(:,stim,iteration) = var(test_rates,[],2)./mean(test_rates,2);
                % now calculate performance on each trial
                false = 0;
                hit = 0;
                

                % now calculate 4 error values per "trial": 2 matrices of "trials" (standard and test) and 2 fits (standard and test)
                
                test_v_test = sum((bsxfun(@minus,test_rates(ordervector,:),fits{stim}(:))).^2);
                test_v_stan = sum((bsxfun(@minus,test_rates(ordervector,:),fits{1}(:))).^2);
                % any time the error of the test_v_test is lower than the
                % test_v_stan, it's a hit:
                
                hit_rate = sum([test_v_test'<test_v_stan']);
                
                stan_v_test = sum((bsxfun(@minus,standard_rates(ordervector,:),fits{stim}(:))).^2);
                stan_v_stan = sum((bsxfun(@minus,standard_rates(ordervector,:),fits{1}(:))).^2);
                fa_rate = sum([stan_v_test'<stan_v_stan']);

                % now use those values to calculate performance:
                
                
                simulated_noise_correlation_ml(rep,stim,iteration,pop) = (nansum(nansum(corr(test_rates')))-num_neurons)/((num_neurons^2)-num_neurons);
                active_performance_ml(rep,stim,iteration,pop) = dPrime2(hit_rate/n_trials,fa_rate/n_trials);
                active_pooling_perfomance_ml(rep,stim,iteration,pop) = dPrime2(poolhit,poolfalse);
                active_sub_pooling_perfomance_ml(rep,stim,iteration,pop) = dPrime2(pospoolhit+negpoolhit,pospoolfalse+negpoolfalse);
                
                % calculate Iole:
                s = [0 1];
                mui = [mean(standard_rates,2) mean(test_rates,2)];
                L = [];
                for i = 1:length(mui); % I'll be obtaining 100 stimulus covariance values
                    cvrn = cov(s,mui(i,:));
                    L(i) = cvrn(2,1);
                end
                L = L';
                Cn = cov(test_rates');
                Cmean = cov(mui');
                R = Cn+Cmean;
                A = (R^-1)*L;
                active_Iole_ml(rep,stim,iteration,pop) = L'*A;
                
            end
        end
        
    end
    colors = [0 0 0; 0.7 0.5 0.5; 0.5 0 1; 0 0 1];

    
    figure
    hold on
    % colors = [0 0 0;1 0 0;1 0.5 0;0.5 0.75 0;0.25 1 0.25;0 0.75 0.25;...
    %     0 0.5 0.75;0 0.25 0.75; 0 0 0.75; 0 0 1; 0.4 0.4 0.4]
    amdepth = [0 0.06 0.16 0.28 0.4 0.6 0.8 1];
    for i = 1:length(corrs)
        plot(amdepth,nanmean(active_performance_ml(i,:,:),3),'color',colors(i,:),'LineWidth',8)
    end
    set(gca,'FontSize',24,'Xtick',[0 0.5 1],'XtickLabel',{'0' '50' '100'})
    hline(1)
    legend(' rsc~0',' rsc~0.1',' rsc~0.2')
    ylim([0 3.5])
    
    
    figure
    hold on
    % colors = [0 0 0;1 0 0;1 0.5 0;0.5 0.75 0;0.25 1 0.25;0 0.75 0.25;...
    %     0 0.5 0.75;0 0.25 0.75; 0 0 0.75; 0 0 1; 0.4 0.4 0.4]
    for i = 1:length(corrs)
        plot(nanmean(active_Iole_ml(i,:,:),3),'color',colors(i,:))
    end
    hline(0)
    
    figure
    hold on
    % colors = [0 0 0;1 0 0;1 0.5 0;0.5 0.75 0;0.25 1 0.25;0 0.75 0.25;...
    %     0 0.5 0.75;0 0.25 0.75; 0 0 0.75; 0 0 1; 0.4 0.4 0.4]
    for i = 1:length(corrs)
        plot(nanmean(active_pooling_perfomance_ml(i,:,:),3),'color',colors(i,:))
    end
    hline(1)
    
        
    figure
    hold on
    % colors = [0 0 0;1 0 0;1 0.5 0;0.5 0.75 0;0.25 1 0.25;0 0.75 0.25;...
    %     0 0.5 0.75;0 0.25 0.75; 0 0 0.75; 0 0 1; 0.4 0.4 0.4]
    for i = 1:length(corrs)
        plot(nanmean(active_sub_pooling_perfomance_ml(i,:,:),3),'color',colors(i,:))
    end
    hline(1)
    
end




