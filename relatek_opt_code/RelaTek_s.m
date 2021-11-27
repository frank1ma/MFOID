function [rdval] = RelaTek_s(dim,t,raw_output,raw_input,srch_opts,crit_opts,win_opts,filt_opts,plot_opts)
%RelaTeK - Tool for Relative Degree Detection
%   This function serves as a tool to detect the relative degree of linear
%   and nonlinear systems by evaluating the frequency response functions of
%   them, which provide an estimated magnitude - frequency plot of the
%   concerced system.
%
%   The power spectral density estimation is based on the welch's method
%   dim       - preserved for future use
%   t         - time sequence
%   raw_output - output collection
%   raw_input  - input signal
%   srch_opts  - search args
%   crit_opts  - criterion args
%   win_opts   - window method
%   filt_opts  - data filter
%
%   Reference:
%       [1] Documentation of pwelch, tfestimation, The MathWorks, Inc.
%
%   date : 06 - 30 - 2019
%   Frank S.Ma


narginchk(6,10)
nargoutchk(0,7)

if nargin < 7
    win_opts = [];
end

if nargin < 8
    filt_opts = [];
end

if isempty(win_opts)
    win_opts = kaiser(size(t,1)/2+1,55);
end

if isempty(filt_opts)
    filt_output_x = raw_output;
    filt_input_u  = raw_input;
else
    filt_output_x = filter(filt_opts,raw_output);   % filter design from designfilt -
    filt_input_u      = filter(filt_opts,raw_input);    % Signal Processing Toolbox
end

% sampling freq.
fs = 1 / (t(2)-t(1));

% threshold of minimum points of segments
min_num_seg = 0.05;

% pass searching parameters
% searching list
% size of search
[srch_range,srch_step] = srch_opts{:};
srch_list = srch_range(1):srch_step:srch_range(2);
srch_len = size(srch_list,2);

% number of outputs in combo
% number of data points
% array for relative degree
% vector def. for cpsd
num_of_points = size(t,1);
% segzp=cell(1,srch_len);
% segdelarr=cell(1,srch_len);
% segarr =cell(1,srch_len);
% max_seg_arr=zeros(srch_len,1);
% psd for u
[psd_u_1,f] = cpsd(filt_input_u(1,:),filt_input_u(1,:),win_opts,[],num_of_points,fs);
[psd_u_2,f] = cpsd(filt_input_u(2,:),filt_input_u(2,:),win_opts,[],num_of_points,fs);
% psd for data
[psd_data1,~] = cpsd(filt_input_u(1,:),filt_output_x(:,1),win_opts,[],num_of_points,fs);
[psd_data2,~] = cpsd(filt_input_u(2,:),filt_output_x(:,1),win_opts,[],num_of_points,fs);
% convert to real positive estf
estf_data1 = (psd_data1./psd_u_1);
estf_data2 = (psd_data2./psd_u_2);
% psd_data = abs(psd_data);
% psd_u    = abs(psd_u);

% Uncomment to show psd of two outputs
% figure(1)
% plot(log10(f),log10(abs(estf_data1)));
% %plot(log10(f),log10(abs(estf_data(:,1)+estf_data(:,2))));
% hold on
% grid on
% figure(2)
% plot(log10(f),log10(abs(estf_data2)));
% hold on
% grid on
% select interval for fitting
num_of_psd = size(psd_data1,1);
crit_flg = crit_opts{1};          % mode selection

switch crit_flg
    
    case 'auto'
        cutoff_start = crit_opts{2}(1);
        cutoff_end = crit_opts{2}(2);
        interval_start = ceil(cutoff_start * num_of_psd);
        interval_end   = floor(cutoff_end * num_of_psd);
        
    case 'manul'
        msgtoken = 0;
        opts1.Interpreter = 'tex';
        opts1.Default = 'Go select';
        ans_select = questdlg('\fontsize{12}Select interval for estimation', ...
            'Selection', ...
            'Go select','Cancel',opts1);
        
        switch ans_select
            
            case 'Go select'               % select intervals
                figure;
                plot(log10(f),log10(abs(estf_data1)));
                grid on;
                title('Please select two points')
                
                while msgtoken == 0        % if not confirmed
                    [x_cor,msgtoken]=manulSelect();
                end
                interval_start = ceil(10^x_cor(1)/(fs/2)* num_of_psd);
                interval_end   = floor(10^x_cor(2)/(fs/2)* num_of_psd);
                
            case 'Cancel'                 % no selection
                disp('No selection of interval. Program terminated.');
                rdval = 0;
                combval = 0;
                return;
        end
        
end
close(gcf);

faxis = log10(f(interval_start:interval_end,1));

% construct matrix for data  num_of_columns by 2
estf_data_in_1 = estf_data1(interval_start:interval_end,:);
estf_data_in_2 = estf_data2(interval_start:interval_end,:);
%plot(faxis,log10(abs((estf_data_in_1))));

% search begins here

%parallel computing

estf_feed_data_1 = log10(abs(estf_data_in_1));
estf_feed_data_2 = log10(abs(estf_data_in_2));
[rd1,~] = getrd(dim,estf_feed_data_1,faxis,min_num_seg,plot_opts,0);
[rd2,~] = getrd(dim,estf_feed_data_2,faxis,min_num_seg,plot_opts,0);

% return
rdval = [rd1 rd2];
end

function [x_cor,msgtoken]=manulSelect()
[x_cor,y_cor] = ginput(2);  % get two points
hold on
points_plot = plot(x_cor,y_cor,'rx'); % show selection
current_ylim = get(gca,'YLim');
line_plot1 = line([x_cor(1) x_cor(1)],[current_ylim(1) current_ylim(2)]);
line_plot2 = line([x_cor(2) x_cor(2)],[current_ylim(1) current_ylim(2)]);
hold off
opts2.Interpreter = 'tex';
opts2.Default = 'Yes';
ans_confirm = questdlg('\fontsize{12}Confirm the interval?', ...
    'Confirmation', ...
    'Yes','No',opts2);

switch ans_confirm
    
    case 'Yes'
        msgtoken = 1;                     % confirm sign
    case 'No'
        delete(points_plot)               % delete plot
        delete(line_plot1)
        delete(line_plot2)
        msgtoken = 0;
end
end



% [EOF]