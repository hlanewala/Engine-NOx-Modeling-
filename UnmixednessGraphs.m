clear
clc
a = 'C:\Users\user\Documents\Master thesis\Cantera Simulation\JT15D\single toroid\Propane\Procedure 6\Data_Reactors.csv';
d = csvread(a);
Row_len = length(d);
Col_len = length(d(1,:));

x = d(2:5, 1);

NOx = d(2:Row_len, [2,4,6,8,10,12]);

CO = d(2:Row_len, [3,5,7,9,11,13]);

PSR_Start_Number = 10;
PSR_End_Number = 19;
PSR_Number = PSR_Start_Number;

Unmixedness = 0.05;

for t = 1:length(NOx(1,:))
    
    legend_numbers = PSR_End_Number - PSR_Start_Number + 1;
    n = 1;
    i = 1;
    legend_names = cell(1,legend_numbers);
    PSR_Number = PSR_Start_Number;
    
    figure('name', ['Unmixedness - ' num2str(Unmixedness)]);
    while 1>0    
        subplot(1,2,1);
        plot(x, NOx(i:i+3,t), '-o');
        title(['Unmixedness - ' num2str(Unmixedness)]);
        xlabel('Thrust Level (%)');
        ylabel('EINO_x (g/kg fuel)');
        legend_names{n} = [num2str(PSR_Number) ' PZ Reactors'];    
        hold all
        if i+5 > length(d)
            hold off 
            legend(legend_names)
        end

        subplot(1,2,2);
        plot(x, CO(i:i+3,t), '-o')
        title(['Unmixedness - ' num2str(Unmixedness)]);
        xlabel('Thrust Level (%)');
        ylabel('EICO (g/kg fuel)');
        legend_names{n} = [num2str(PSR_Number) ' PZ Reactors'];    
        hold all

        n = n+1;

        PSR_Number = PSR_Number + 1;
        i = i+5;
        if i > length(d)        
            break
        end
    end
    Unmixedness = Unmixedness+0.05;
    legend(legend_names)
end
