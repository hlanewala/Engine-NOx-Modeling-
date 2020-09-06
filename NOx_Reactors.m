clear
clc
a = pwd;
a1 = '\Data_Reactors.csv';
a = [a a1];
d = csvread(a);
Row_len = length(d(:,1));
Col_len = length(d(1,:));
s = 0;
for g = 1:6
    s = s + 0.05;
    UM = 2*g;
    x = zeros(20,1);
    Idle = zeros(20,1);
    Landing = zeros(20,1);
    ClimbOut = zeros(20,1);
    TakeOff = zeros(20,1);
    i = 2;
    Counter = 0;
    while i <= Row_len
       Counter = Counter+1;
       x(Counter) = d(i,16);
       Idle(Counter) = d(i,UM);
       Landing(Counter) = d(i+1,UM);
       ClimbOut(Counter) = d(i+2,UM);
       TakeOff(Counter) = d(i+3,UM);
       i = i+5;  
    end
    x = x(1:Counter);
    Idle = Idle(1:Counter);
    Landing = Landing(1:Counter);
    ClimbOut = ClimbOut(1:Counter);
    TakeOff = TakeOff(1:Counter);

    figure('name', ['Unmixedness - ' num2str(s)])
    plot(x,Idle,'-+',x,Landing,'-*',x,ClimbOut,'-s', x,TakeOff, '-o')
    title('Effect of Primary Zone size on NOx @ 12.5% reacting dilution zone air & 100% reacting cooling flow');
    xlabel('Primary Zone length');
    ylabel('EINO_x (g/kg fuel');
    legend('7%', '30%', '85%', '100%');
end









