clear
clc
Unmixedness = 0.05; % Unmixedness parameter
NReactors = 10; %Number of Reactor in the normal distribution
mus = [0.497, 0.541, 0.808, 0.899];
Conditions = [7, 30, 85, 100];

problem = rem(NReactors, 2);
if problem == 1
    problem = 'odd';
else
    problem = 'even';
end

%%%%%%%%Input%%%%%%%%%%
while Unmixedness <= 0.3 %Running the loop for a range of unmixedness 
    % values
    X = linspace(0, 2, 500);  

    k = 2;
    for j = 1:(length(mus)) %Running the loop for each thrust condition         
        standard = 'y'; 

        mu = mus(j);
        Sigma = mu*Unmixedness;

    %     figure('name', '1')
    %     if standard == 'n'
    %         Y = normpdf(X,mu,Sigma);
    %         plot(X,Y);      
    %     else
    %         p = (X-mu)./Sigma;
    %         Y = normpdf(p);
    %         plot(p,Y);
    %         axis([-5 5 0 0.4])      
    %     end
    %     grid on

    %     figure('name', '2')
    %     if standard == 'n'
    %         Y = normcdf(X,mu,Sigma);    
    %         plot(X,Y);
    %     else
    %         p = (X-mu)./Sigma;
    %         Y = normcdf(p);
    %         plot(p,Y);
    %         axis([-5 5 0 1])
    %     end
    %     grid on

        g = normcdf([-3 3]); %Cumulative distribution function is from -3 
        % to 3 since this range includes 99.9% of the flow/
        Reactors = 6/(NReactors); %Reactors can be thought of as Z in
        % the context of normalized equivalence ratio (Phi)
        Reac = Reactors/2; %Reac can also be thought of as Z
        i = 1;
        h = 1:i; %Array holding the amount of flow in reactors on either 
        % sides of the mean and at equal distance from the mean reactor.
        % This amount is found from the cumulative distribution function.        
        Percentages = zeros([NReactors,1,1]); %Array holding the amount of 
        % flow for each reactor at a particular phi. This amount is 
        % found by dividing h by 2.
        Phis = zeros([NReactors,1,1]); %Array holding the phi value for 
        % each reactor
        
        for i=1:NReactors
            if strcmp(problem,'odd')
                if i == 1
                    mid = (NReactors+1)/2; %Calculation of the middleth 
                    % reactor for an odd number of total reactors                    
                    g = normcdf([-Reac Reac]);
                    h(i) = g(2)-g(1);
                    Phis(mid) = mu;
                    Percentages(mid) = h(i);
                elseif i > 1 && i < (mid+1)                    
                    g = normcdf([-(Reac+Reactors*(i-1)) (Reac+Reactors*(i-1))]);
                    h(i) = (g(2)-g(1))-f;                    
                    PhiPlus = Reactors*(i-1);
                    Phis((mid-1)+i) = Sigma*PhiPlus+mu;
                    if i < mid
                        Percentages((mid-1)+i) = h(i)/2;
                    elseif i == mid
                        Percentages((mid-1)+i) = (1-c)/2;
                    end
                end
                
            else
                if i == 1                    
                    mid = (NReactors)/2;
                    g = normcdf([-Reactors Reactors]);
                    h(i) = g(2)-g(1);
                    PhiPlus = Reac*i;                     
                    Phis(mid) = Sigma*(-PhiPlus)+mu;                  
                    Phis(mid+1) = Sigma*(PhiPlus)+mu;                    
                    Percentages(mid) = h(i)/2;
                    Percentages(mid+1) = h(i)/2;  
                elseif i > 1 && i < (mid+1)
                    g = normcdf([-(Reactors*i) (Reactors*i)]);
                    h(i) = (g(2)-g(1))-f;
                    PhiPlus = Reac+Reactors*(i-1);
                    Phis(mid+i) = Sigma*PhiPlus+mu;
                    if i < mid
                        Percentages(mid+i) = h(i)/2;
                    elseif i == mid
                        Percentages(mid+i) = (1-c)/2;
                    end
                end
            end
            
            if i > 1 && i < (mid+1)
                Phis((mid+1)-i) = Sigma*(-PhiPlus)+mu;
                if i < mid
                    Percentages((mid+1)-i) = h(i)/2;
                elseif i == mid
                    Percentages((mid+1)-i) = (1-c)/2;
                end
                c = sum(Percentages);
            end
            
            f = sum(h);
            
            if (f > 1)
                break
            end            
        end
        
        PhisAndFlows(1,1) = Unmixedness; %Array that will eventually be 
        % outputted in to a csv file for each unmixedness value.
        PhisAndFlows(1,k) = Conditions(j);
        PhisAndFlows(2:length(Percentages)+1,1) = Percentages;
        PhisAndFlows(2:length(Phis)+1,k) = Phis;

        k = k+1;
    end
    s = pwd;
    s1 = '\PhisAndFlows_';
    s2 = num2str(Unmixedness);
    s3 = '_';
    s4 = num2str(NReactors);
    s5 = '.csv';
    s1 = [s s1 s2 s3 s4 s5];
    dlmwrite(s1, PhisAndFlows, 'delimiter', ',');
    PhisAndFlows
    Unmixedness = Unmixedness + 0.05;
end


