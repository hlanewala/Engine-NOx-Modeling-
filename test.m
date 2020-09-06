clear
clc

NReactors = 16;
a = 'C:\Users\user\Documents\Master thesis\Cantera Simulation\PW305\Propane\Procedure 6 - New Inputs\Parametric Study Tests\';
a1 = num2str(NReactors);
a2 = ' Reactors';
a = [a a1 a2];
try
    cd(a)
catch err
    break
end  

Unmixedness = 0.05;

a1 = '\Unmixedness - ';
a2 = num2str(Unmixedness);
a3 = '\Emissions.csv';
b = [a a1 a2 a3];
d = xlsread(b);
disp(d);

