
arterial_dimension = readtable('artery_size.csv');
rho = 1.05; %blood density [g/cm^3]
rho = 1000*rho; %change to [kg/m^3]
mu = 0.04; %blood dynamic viscosity [poise] %0.04 initial
mu = 0.1*mu; %change to [Pa s]

Length = arterial_dimension(:,'Length');
Radius = arterial_dimension(:,'Radius');
Strength = arterial_dimension(:,'millionEDyn_cm');
Thickness = arterial_dimension(:,'Thickness');

Length = table2array(Length);
Radius = table2array (Radius);
Strength =table2array (Strength);
Thickness =table2array(Thickness);

% Parameters for different segments: length(l), radius(r), thickness(t),and arterial strength(E).
% Since some segments are combinations of several segments, r,t, and E are averaged out.

l = zeros(1,9);
l(1)= sum(Length([1,2]));
l(2)= Length(3);
l(3)= Length(4);
l(4)= sum(Length( [5 6]));
l(5)= sum(Length ([7 8]));
l(6)= sum(Length ([9 10 11]));
l(7)= sum(Length ([12 13]));
l(8)= sum(Length ([14 15]));
l(9)= Length (16);
l = l*1e-2; %Unit conversion


r = zeros(1,9);
r(1)= sum(Radius([1,2]).*Length([1,2]))/l(1);
r(2)= Radius(3);
r(3)= Radius(4);
r(4)= sum(Radius ([5,6]).*Length([5,6]))/l(4) ; 
r(5)= sum(Radius ([7,8]).*Length([7,8]))/l(5) ;
r(6)= sum(Radius ([9,10,11]).*Length([9,10,11]))/l(6);
r(7)= sum(Radius ([12,13]).*Length([12,13]))/l(7);
r(8)= sum(Radius([14,15]).*Length([14,15]))/l(8);
r(9)= Radius(16);
r = r*1e-2;

t = zeros(1,9);
t(1)= sum(Thickness([1,2]).*Length([1,2]))/l(1) ;
t(2)= Thickness(3);
t(3)= Thickness(4);
t(4)= sum(Thickness ([5,6]).*Length([5,6]))/l(4) ; 
t(5)= sum(Thickness ([7,8]).*Length([7,8]))/l(5) ;
t(6)= sum(Thickness ([9,10,11]).*Length([9,10,11]))/l(6);
t(7)= sum(Thickness ([12,13]).*Length([12,13]))/l(7);
t(8)= sum(Thickness([14,15]).*Length([14,15]))/l(8);
t(9)= Thickness (16);
t = t*1e-2;

E = zeros(1,9);
E(1)= sum(Strength([1,2]).*Length([1,2]))/l(1) ;  
E(2)= Strength(3);
E(3)= Strength(4);
E(4)= sum(Strength ([5,6]).*Length([5,6]))/l(4) ; 
E(5)= sum(Strength ([7,8]).*Length([7,8]))/l(5) ;
E(6)= sum(Strength ([9,10,11]).*Length([9,10,11]))/l(6);
E(7)= sum(Strength ([12,13]).*Length([12,13]))/l(7);
E(8)= sum(Strength([14,15]).*Length([14,15]))/l(8);
E(9)= Strength(16);
E = E*1e-3; % TRANSFER TO N/m

R_cal = l.*2*(2+9)*pi*mu./((pi*r.^2).^2);           %Resistance
L_cal = l.*rho./(pi*r.^2);                          %Inductance
c_cal = l.*2.*(pi*r.^2).^(3/2)./(4/3*sqrt(pi)*E.*t);%Capacitance
R_t_cal = 0.2*mean(R_cal);                          %Terminal Resistance
clear rho mu E r l h


save('initial_model_parameters.mat')
%%
%This part is used to tune the resistance, inductance, and capacitance.
%The purpose of Tuning is to investigate the influences of RLC changes.
%Tuning RLC helps to simulate the unique peek in the reference.

R=R_cal ;R_t = R_t_cal;
L = L_cal;
c=c_cal;
R_ratio = 100;
R=R_cal * R_ratio;
R_t = R_t_cal *R_ratio;
L= L_cal*1.3; 
c= c_cal*0.5; 

load '2019-07-15.mat'
ICP=d.icp.wave;
pvdf=d.pvdf_test.wave;
co = d.co.wave;
abp= d.abp.wave;
cvp = d.cvp.wave;
C1 =-1/(c(5)+c(6)+c(7)+c(8)+c(9));
C2 = -1/(c(1)+c(2));
C3 = -1/(c(3)+c(4));

%The state-space system:
%X_dot = Ax+B, Y=Cx+D
%  x here is blood flow rate and blood pressure in each segment

A= [0 0 0 0 0 0 0 C1 C1 C1 C1 C1;...  
  0 0 0 C2 C2 0 0 -C2 0 0 0 0;...
  0 0 0 0 0 C3 C3 0 -C3 0 0 0;...
  0 1/L(1) 0 -(R(1)+R_t)/L(1) 0 0 0 0 0 0 0 0;...
  0 1/L(2) 0  0 -(R(2)+R_t)/L(2) 0 0 0 0 0 0 0;...
  0 0 1/L(3) 0 0 -(R(3)+R_t)/L(3) 0 0 0 0 0 0;...
  0 0 1/L(3) 0 0 0 -(R(4)+R_t)/L(4) 0 0 0 0 0;...
  1/L(5) -1/L(5) 0 0 0 0 0 -R(5)/L(5) 0 0 0 0;...
  1/L(6) 0 -1/L(6) 0 0 0 0 0 -R(6)/L(6) 0 0 0;...
  1/L(7) 0 0 0 0 0 0 0 -R(7)/L(7) 0 0 0;...
  1/L(8) 0 0 0 0 0 0 0 0 0 -R(8)/L(8) 0;...
  1/L(9) 0 0 0 0 0 0 0 0 0 0 -R(9)/L(9)];
B= [-C1; zeros(11,1)];
C=[0 1/2 1/2 0 0 0 0 0 0 0 0 0];
sys = ss(A,B,C,0);
b = 1080000:1280000; %b is the steady-state period
abp_use = abp(b);
co_use = co(b);      % the data sets which I tested 871000:900000 ; 1080000-1280000;
ICP_use = ICP(b);    % I found that this function has overshoot and it takes a long to stablize
abp_used = interp(abp_use,100);  % the response
ICP_used = interp(ICP_use,100);
time = 0:1/200:(length(abp_used)-1)/200; %The sample rate is 200
Y = lsim(sys, abp_used, time);     %Simulated ICP: the output of the system

%%
%This block is used to evaluate the performance of the system, since the system output has different numerical unit
%I calculated the error rate here.

block = 3556000:9564000;   %Select a short period of time to calculate the required parameters
delta_t=0.2;               %This is the delay time comes from the system
ICP_sim = Y(block);                     %arterial pressure in steady-state
ICP_r = ICP_used(block);                %used block here
ICP_std = std(ICP_r,"omitnan");         %reference standard deviation
sim_ICP_std = std(ICP_sim,"omitnan");   %output standard deviation
ICP_mean = mean(ICP_r,"omitnan");       %reference average
sim_ICP_mean = mean(ICP_sim,"omitnan"); %output average
Ratio = ICP_std/artery_std;             %which will be used to correct the numerical difference between the reference and output

%The output(Y) is simply a sine wave that started from 0, however, the reference is not a sine wave
%To best simulate the ICP waveform, the output(Y) needs to subtract the mean at first and then multiply the Ratio to
%make sure the output and reference have the same numerical level. At last, add ICP_mean
%The main purpose of this is to see if my model is able to simulate the ICP fluctuation which is the main problem in the emergency

sim_ICP_ratioed = Ratio .* ( Y(block-delta_t*20000)-sim_ICP_mean )+ ICP_mean;
Error = (sim_ICP_ratioed - ICP_r).^2;
Error_mean = mean(Error);


%%
%The section below is used to plot

figure(8)
subplot(3,2,1)
plot(ICP_r);%140000:350000
hold on 
plot(204000:504000,ICP_r(204000:504000))
hold on 
plot(5804000:6008000,ICP_r(5804000:6008000))
title('ICP segment');
legend ('Unstable period','Segement1','Segement1')
subplot(3,2,2)
plot(artery_r);     %artery_ratioed 140000:350000 Ratioed simulated waveform segment
hold on 
plot(204000:504000,artery_r(204000:504000))
hold on 
plot(5804000:6008000,artery_r(5804000:6008000))
title('Simulated ICP segment ');
legend ('Unstable period','Segement1','Segement1')

subplot(3,2,3)
plot(ICP_r(204000:504000)); 
title('ICP segment(1)');
subplot(3,2,4)
plot(artery_ratioed(204000:504000));
title('Ratioed simulated waveform segment(1) ');

subplot(3,2,5)
plot(ICP_r(5804000:6008000)); 
title('ICP segment(2)');
subplot(3,2,6)
plot(artery_ratioed(5804000:6008000));
title('Ratioed simulated waveform segment(2) ');

figure (9)
plot (ICP_r(4000:204000));
hold on 
plot(artery_ratioed(4000:204000));

legend ('Reference','Ratioed output')
