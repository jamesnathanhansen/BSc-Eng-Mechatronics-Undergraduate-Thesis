t = [0:1:4999];
x = [0:2.5/999:2.5];
y = ones(1,3000)*2.5;
z = [2.5:-2.5/999:0];
run = [x,y,z];
run = awgn(run,20);

figure(1),clf
plot(t,run)
grid on
xlabel('Number of Samples')
ylabel('Mass (kg)')