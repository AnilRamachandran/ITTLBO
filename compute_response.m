function [y_response]=compute_response(y_response,mean,var,funcstr)
%%update y values

y_response=(y_response-mean)./var;
