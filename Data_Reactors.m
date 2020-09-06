clear
clc
Data = zeros(1000,16);
UM = 0.05;
NReactors = 10;

for i=1:6
    Data(1,i*2) = i*UM;
end

RowCounter = 1;
CurrentRow = 2;
while 1 > 0    
    a = 'C:\Users\user\Documents\Master thesis\Cantera Simulation\PW305\Propane\Procedure 6 - New Inputs\Parametric Study Tests\';
    a1 = num2str(NReactors);
    a2 = ' Reactors';
    a = [a a1 a2];
    try
        cd(a)
        RowCounter = RowCounter+5;
    catch err
        break
    end  
    
    Unmixedness = 0.05;
    CurrentColumn = 2;
    
    while Unmixedness <= 0.3        
        a1 = '\Unmixedness - ';
        a2 = num2str(Unmixedness);
        a3 = '\Emissions.csv';
        b = [a a1 a2 a3];
        try            
            d = xlsread(b);
        catch err1
            DataString = 'There seems to be no data for ';
            DataString1 = num2str(NReactors);
            DataString2 = ' Reactors';
            DataString = [DataString DataString1 DataString2];            
            disp(DataString);
            CurrentColumn = 14;
            break
        end
        Data(CurrentRow:CurrentRow+(length(d(:,1))-1), 1) = d(1:length(d(:,1)),1);
        Data(CurrentRow:CurrentRow+(length(d(:,1))-1), CurrentColumn:CurrentColumn+1) = d(1:length(d(:,1)),2:3);
        Unmixedness = Unmixedness+0.05; 
        CurrentColumn = CurrentColumn+2;    
    end
    Data(CurrentRow,CurrentColumn+2) = NReactors;
    CurrentRow = CurrentRow + 5;
    NReactors = NReactors + 1;
    
end
Data = Data(1:RowCounter, 1:16);

s = 'C:\Users\user\Documents\Master thesis\Cantera Simulation\PW305\Propane\Procedure 6 - New Inputs\Parametric Study Tests';
s1 = '\Data_Reactors.csv';
s1 = [s s1];
dlmwrite(s1, Data, 'delimiter', ',');

