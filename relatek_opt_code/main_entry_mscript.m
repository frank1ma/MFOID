
%  Flat Output Identification - MIMO (FOID-M) 
%  Date : 02 - 19 - 2021
%  Frank (Shangjie) Ma
%  ----------------------------------------------------------------------

clear 
clc
swEPSfigure
swFigSize


% parallel cpu computing
% localpc = parcluster('local');
% parpool('local',localpc.NumWorkers)
%%
%Simulation Config.
%
% sampling
fs = 1000;  % Sample frequency
tval = 20; % time for simulation

% reference
dt = 1 / fs;         % step siz
t = (0:dt:tval-dt)'; % time vector
refin =  0.1*randn(2,size(t,1));   % random input
%load('my_input.mat')
% linear output matrix
C = [1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0];
% C = [1 0 0 0 0 0;
%     0 0 1 0 0 0;
%     0 0 0 0 1 0];
%%
% Model Selection for simulation

%[out,uin,t,B] = rdid_nonlnr_mdl_4d(fs,tval,C,refin);
[out,uin,t,B] = rdid_nonlnr_mdl_msd(fs,tval,C,refin);
% [out,uin,t] = rdid_lnr_mdl_ivp(fs,tval,C,refin);      % inverted pendulum
% [out,uin,t] = rdid_lnr_mdl_3dsys(fs,tval,C,refin);    % 3d system
% [out,uin,t] = rdid_lnr_mdl_3dsp(fs,tval,C,refin);     % special 3d system
% [out,uin,t] = rdid_lnr_mdl_flexible_link(fs,tval,C,refin);

%%
% Args Configuration for y1

% number of dimensions not required
% it is optional in function getrd to limit the highest relative degree 
n1 = 2;
col_opts = 3;
raw_data = out';
raw_input = uin;

% search settings for combination
srch_opts= {[-50,50],1/100};  % {[start end],search step}

% data interval args
% mode = 'auto' or 'manul'
% for 'auto', a interval should be specified 
% arg list: {{mode,[interval(%)]}}
crit_opts = {'manul',[]};
%crit_opts = {'auto',[0.001 0.06]};

% optional windowing function. Default is Kasiver with beta 35
% win_opts = kaiser(size(t,1)/2+1,35);

win_opts = kaiser(1000,35);

% optional filter applying to raw output and input. Defalut empty.
% see MATLAB Signal Processing Toolbox
filt_opts =[];
%
plot_opts = -0.1;

[rdval] = RelaTek_s(n1,t,raw_data,raw_input,srch_opts,crit_opts,win_opts,filt_opts,plot_opts);

%S = -inv([B(5,1) B(6,1);B(5,2) B(6,2)])*[B(4,1);B(4,2)];
%%
coeff_reg =[-4:1:6;-4:1:6];
size_coeff_reg = size(coeff_reg,2);
pair = zeros(4,size_coeff_reg);
interval=[];
manual_for_all_switch = 0;
% Args Configuration for y2
n2 = 4;
for i = 1:size_coeff_reg
    [rdval_1,ptgval_1,rdval_2,ptgval_2,combval,interval] = RelaTek_m(n2,t,raw_data,raw_input,srch_opts,crit_opts,win_opts,filt_opts,plot_opts,col_opts,coeff_reg(:,i),interval);
    if (~manual_for_all_switch) && (~isempty(interval))
        crit_opts = {'skip',[]};
    end
    %plot_result(rdval_1,ptgval_1,combval)
    %plot_result(rdval_2,ptgval_2,combval)
    pair(:,i) = [coeff_reg(1,i);plot_result(rdval_1,ptgval_1,combval);coeff_reg(2,i);plot_result(rdval_2,ptgval_2,combval)];
end

x=-10:0.1:10;
    
p1=polyfit(pair(1,:),pair(2,:),1);
plot(x,p1(1)*x+p1(2),'-r')
hold on
plot(pair(1,:),pair(2,:),'xr')
p2=polyfit(pair(3,:),pair(4,:),1);
plot(x,p2(1)*x+p2(2),'--b')
plot(pair(3,:),pair(4,:),'sb')
xlabel('$c_2$','Interpreter','Latex')
ylabel('$c_3$','Interpreter','Latex')
%plot(-1.0459,-0.63114,'mo')
syms t
eqn = p1(1)*t+p1(2)-p2(1)*t-p2(2)==0
V1=vpasolve(eqn,t,[-20 20])
V2=p1(1)*V1+p1(2)
%%
% Plots
%========================================================================
%plot_result(rdval_1,ptgval_1,combval)
S = -inv([B(5,1) B(6,1);B(5,2) B(6,2)])*[B(4,1);B(4,2)]