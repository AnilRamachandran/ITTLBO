function [mean_cal,sigma_cal]=find_mean_var(y_fun,str)
%%finding mean and var to update y values

mean_cal=mean(y_fun);
sigma_cal=std(y_fun);


