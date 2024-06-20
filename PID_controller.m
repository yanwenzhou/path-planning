function [steering_angle] = PID_controller(P,I,D,target_angle,current_angle)
%PID_CONTROL 此处显示有关此函数的摘要
%   此处显示详细说明
persistent self_integral;
if isempty(self_integral)
   self_integral=0;
end
persistent self_prev_error;
if isempty(self_prev_error)
   self_prev_error=0;
end
error = target_angle-current_angle; % 计算误差
proportional = P*error;             % 计算比例项
self_integral  = self_integral+error;                 % 积分项等于误差累加
integral = I * self_integral;            % Calculate the integral term
derivative = D * (error - self_prev_error); %计算微分
self_prev_error = error;            % 赋值上一个误差
% 计算转向角
steering_angle = proportional + integral + derivative;
end

