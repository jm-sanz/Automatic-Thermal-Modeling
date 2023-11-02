% Date: 31.12.2023
% Author: José Miguel Sanz Alcaine (jm_sanz@unizar.es)
% Code for replicating the results from the work: Estimation of Semiconductor Power Losses Through Automatic Thermal Modeling 
% Sampling time is Ts= 100ms
clear all, close all, clc

% Load experiments data (MATLAB structures format)
load("data\test_sat.mat", "test_sat");
load("data\test_lin.mat", "test_lin");
load("data\test_driver.mat", "test_driver");
load("data\test_inductor.mat", "test_inductor");

% Concatenate structures
test_array = [test_sat, test_lin, test_driver, test_inductor];
number_test = length(test_array);
r_acc = zeros(1,number_test);

n = 6; % Number of power variables
m = 8; % Number of temperature variables

% Format data into Powers and Temperature matrices
for i=1:number_test
    if(i==1)
        r_acc(i)=test_array(i).size;
    else
        r_acc(i)=r_acc(i-1) + test_array(i).size;    
    end
end

Powers = zeros(n,r_acc(end));
Temperatures = zeros(m,r_acc(end));

for i=1:number_test

    test_data = test_array(i);
    if(i == 1)
        initial_ix = 1;
    else
        initial_ix = r_acc(i-1) + 1;
    end

    final_ix = r_acc(i);


    Temperatures(1,initial_ix:final_ix) = test_data.state.temp_hs - test_data.state.temp_amb;
    Temperatures(2,initial_ix:final_ix) = test_data.state.temp_ls- test_data.state.temp_amb;
    Temperatures(3,initial_ix:final_ix) = test_data.state.temp_drv_hs - test_data.state.temp_amb;
    Temperatures(4,initial_ix:final_ix) = test_data.state.temp_drv_ls - test_data.state.temp_amb;
    Temperatures(5,initial_ix:final_ix) = test_data.state.temp_vin - test_data.state.temp_amb;
    Temperatures(6,initial_ix:final_ix) = test_data.state.temp_vsw - test_data.state.temp_amb;
    Temperatures(7,initial_ix:final_ix) = test_data.state.temp_gnd - test_data.state.temp_amb;
    Temperatures(8,initial_ix:final_ix) = test_data.state.temp_inductor - test_data.state.temp_amb;


       if(i ~= 1)
        delta = (Temperatures(1:m,initial_ix - 1) - Temperatures(1:m,initial_ix));
        Temperatures(1:m,initial_ix:final_ix) = Temperatures(1:m,initial_ix:final_ix) + delta;
       end

    Powers(1,initial_ix:final_ix) = test_data.input.power_transistorh;
    Powers(2,initial_ix:final_ix) = test_data.input.power_transistorl;
    Powers(3,initial_ix:final_ix) = test_data.input.power_drvh;
    Powers(4,initial_ix:final_ix) = test_data.input.power_drvl;
    Powers(5,initial_ix:final_ix) = test_data.input.power_trackh;
    Powers(6,initial_ix:final_ix) = test_data.input.power_inductor;
end

%% TRAIN AND TEST WITH THE SAME DATA
% State-Space Definitions
u_k = (Temperatures(:,1:end-1));
u_k_1 = Temperatures(:,2:end);
x_k = (Powers(:,1:end-1));
z_k = [u_k; x_k];

% Least Square by means of pseudoinverse
W = u_k_1*pinv(z_k);
A = W(:,1:m);
B = W(:,m+1:end);

% Least Square by means of pseudoinverse with regularization
% epsilon = 1e2;
% W = u_k_1*z_k'*pinv(z_k*z_k' + epsilon*eye(n+m));
% A = W(:,1:m);
% B = W(:,m+1:end);

% Open Loop Simulation of Temperature evolution with the same data
u_sim = zeros(8,length(Temperatures));
u_sim(:,1) = Temperatures(1:m,1);

for i=1:length(Temperatures)-1
    u_sim(:,i+1) = A*u_sim(:,i) + B*Powers(:,i);
end

% Power Estimation
u_k_filt = zeros(m, length(u_k));
u_k_1_filt = zeros(m, length(u_k_1));

% Filtering moving average 
k = 500; % Window Length
for i=1:m
    u_k_filt(i,:) = movmean(u_k(i,:),k);
    u_k_1_filt(i,:) = movmean(u_k_1(i,:),k);
end

% Calculation (Eq. 9)
x_sim = pinv(B)*(u_k_1_filt - A*u_k_filt);

% Plots
figure(1)
subplot(2,1,1)
plot(u_sim')
xlabel('Samples (n)')
ylabel('\Delta{\it T} (ºC)')
title('Results tested with training data')
subplot(2,1,2)
plot(x_sim')
xlabel('Samples (n)')
ylabel('Power (W)')
%% TRAIN AND TEST WITH DIFFERENT DATA 

% Selection of data for testing
sat_test_range = 2200000:2591880;
lin_test_range = 2872000:2949340;
drv_test_range = 3555000:3668000;
ind_test_range = 4044000:4136600;

% Saving test data
Temp_sat_test = Temperatures(:,sat_test_range);
Temp_lin_test = Temperatures(:,lin_test_range);
Temp_drv_test = Temperatures(:,drv_test_range);
Temp_ind_test = Temperatures(:,ind_test_range);

Powers_sat_test = Powers(:,sat_test_range);
Powers_lin_test = Powers(:,lin_test_range);
Powers_drv_test = Powers(:,drv_test_range);
Powers_ind_test = Powers(:,ind_test_range);

% Removing this data from training
Temperatures(:, [sat_test_range, lin_test_range, drv_test_range, ind_test_range]) = [];
Powers(:, [sat_test_range, lin_test_range, drv_test_range, ind_test_range]) = [];

% State-Space Definitions 
u_k = (Temperatures(:,1:end-1));
u_k_1 = Temperatures(:,2:end);
x_k = (Powers(:,1:end-1));
z_k = [u_k; x_k];

% Least Square by means of pseudoinverse
W = u_k_1*pinv(z_k);
A = W(:,1:m);
B = W(:,m+1:end);

% Least Square by means of pseudoinverse with regularization
% epsilon = 1e2;
% W = u_k_1*z_k'*pinv(z_k*z_k' + epsilon*eye(n+m));
% A = W(:,1:m);
% B = W(:,m+1:end);

% Open loop simulation of temperature evolution with test data

% Test vectors definitions
u_sim_sat = zeros(m,length(Temp_sat_test));
u_sim_sat(:,1) = Temp_sat_test(:,1);

u_sim_lin = zeros(m,length(Temp_lin_test));
u_sim_lin(:,1) = Temp_lin_test(:,1);

u_sim_drv = zeros(m,length(Temp_drv_test));
u_sim_drv(:,1) = Temp_drv_test(:,1);

u_sim_ind = zeros(m,length(Temp_ind_test));
u_sim_ind(:,1) = Temp_ind_test(:,1);

% Simulation of saturations tests
for i=1:length(Temp_sat_test)-1
    u_sim_sat(:,i+1) = A*u_sim_sat(:,i) + B*Powers_sat_test(:,i);
end

% Simulation of linear tests
for i=1:length(Temp_lin_test)-1
    u_sim_lin(:,i+1) = A*u_sim_lin(:,i) + B*Powers_lin_test(:,i);
end

% Simulation of driver tests
for i=1:length(Temp_drv_test)-1
   u_sim_drv(:,i+1) = A*u_sim_drv(:,i) + B*Powers_drv_test(:,i);
end

% Simulation of inductor tests
for i=1:length(Temp_ind_test)-1
    u_sim_ind(:,i+1) = A*u_sim_ind(:,i) + B*Powers_ind_test(:,i);
end

% Power estimation
u_k_sat = Temp_sat_test(:,1:end-1);
u_k_lin = Temp_lin_test(:,1:end-1);
u_k_drv = Temp_drv_test(:,1:end-1);
u_k_ind = Temp_ind_test(:,1:end-1);

u_k_1_sat = Temp_sat_test(:,2:end);
u_k_1_lin = Temp_lin_test(:,2:end);
u_k_1_drv = Temp_drv_test(:,2:end);
u_k_1_ind = Temp_ind_test(:,2:end);

u_k_sat_filt = zeros(m, length(u_k_sat));
u_k_lin_filt = zeros(m, length(u_k_lin));
u_k_drv_filt = zeros(m, length(u_k_drv));
u_k_ind_filt = zeros(m, length(u_k_ind));

u_k_1_sat_filt = zeros(m, length(u_k_sat));
u_k_1_lin_filt = zeros(m, length(u_k_lin));
u_k_1_drv_filt = zeros(m, length(u_k_drv));
u_k_1_ind_filt = zeros(m, length(u_k_ind));

% Filtering moving average 
k = 500; % Window Length
for i=1:m
    u_k_sat_filt(i,:) = movmean(u_k_sat(i,:),k);
    u_k_1_sat_filt(i,:) = movmean(u_k_1_sat(i,:),k);

    u_k_lin_filt(i,:) = movmean(u_k_lin(i,:),k);
    u_k_1_lin_filt(i,:) = movmean(u_k_1_lin(i,:),k);

    u_k_drv_filt(i,:) = movmean(u_k_drv(i,:),k);
    u_k_1_drv_filt(i,:) = movmean(u_k_1_drv(i,:),k);

    u_k_ind_filt(i,:) = movmean(u_k_ind(i,:),k);
    u_k_1_ind_filt(i,:) = movmean(u_k_1_ind(i,:),k);
end
% Calculation (Eq. 9)
x_sim_sat = pinv(B)*(u_k_1_sat_filt - A*u_k_sat_filt);
x_sim_lin = pinv(B)*(u_k_1_lin_filt - A*u_k_lin_filt);
x_sim_drv = pinv(B)*(u_k_1_drv_filt - A*u_k_drv_filt);
x_sim_ind = pinv(B)*(u_k_1_ind_filt - A*u_k_ind_filt);


% Plots (case of sat. as example)
figure(2)
subplot(2,1,1)
plot(u_sim_sat')
xlabel('Samples (n)')
ylabel('\Delta{\it T} (ºC)')
title('Results tested with data different from training data (Saturation case)')
subplot(2,1,2)
plot(x_sim_sat')
xlabel('Samples (n)')
ylabel('Power (W)')
