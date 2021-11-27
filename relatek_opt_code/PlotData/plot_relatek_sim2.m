% Plot for Relatek
% Simulation 1 - Inverted Pendulum 
% Shangjie(Frank) Ma
% 10-24-2019

%
load('sim2_result')
load('sim2_combval')
load('sim2_rdval')
load('sim2_ptgval')

% format setting
swEPSfigure
swFigSize
TitleFlag = 0;   % 1 to show title, 0 hide.
GridFlag  = 1; % 1  to turn grid on, 0 off.

% find best ratio and its max relative degree in the charts
% do not change
maxrd = max(rdval);             % max rd
minrd = min(rdval);             
num_of_max_rd = find(rdval==maxrd); % possible ratio for max rd
 ptgcfdt= ptgval(num_of_max_rd)./sum(ptgval(num_of_max_rd));
 idx_of_max_ptgcfdt = find(ptgcfdt==max(ptgcfdt));
         
target_list = combval(num_of_max_rd(1) + idx_of_max_ptgcfdt - 1);
[~,target_idx] = min(abs(target_list));
target_ratio = target_list(target_idx);
target_ratio
maxrd

%% sim2 Closed-loop Response 
f1=figure(1);
% subplot 1
f1Ax1=subplot(211);
fSim1ClosedloopResponseX1=plot(t,out(1,:));
hold on;
plot(t(1),out(1,1),'rx')
legend('response','initial value','Orientation','horizontal')
hold off;
ylabel('$\theta_{\delta}(rad)$')
xlabel('time$(s)$')
set(f1Ax1,'YLim',[-2e-3 3.5e-3]);      % axes for subplot 1

if (TitleFlag)                   % title
    title('Simulation 1 Closed-loop Response-$\theta_{\delta}$');
end
if (GridFlag)                    % grid 
    grid on
end
% subplot 2 
f1Ax2=subplot(212);
fSim1ClosedloopResponseX3=plot(t,out(2,:));
hold on;
plot(t(1),out(2,1),'rx')
legend('response','initial value','Orientation','horizontal')
hold off;
ylabel('$\phi_{\delta}(rad)$')
xlabel('time$(s)$') 
set(f1Ax2,'YLim',[-5e-3 7e-3]);    % axes for subplot 2

if (TitleFlag)
    title('Simulation 1 Closed-loop Response-$\phi_{\delta}$');
end
if (GridFlag)
    grid on
end

%% sim2 relative degree global chart
f2=figure(2);
plot(combval,rdval)
ylabel('Relative Degree($r$)')
xlabel('Ratio($C_r$)')
% axes
f2ax = f2.CurrentAxes;
set(f2ax,'YLim',[0 4]);   
% title
if (TitleFlag)            
    title('Relative Degree');
end
% grid 
if (GridFlag)             
    grid on
end 
% arrow position
 arrow_start_cor = [target_ratio+0.2 maxrd+0.3];   
 arrow_end_cor =[target_ratio+0.02 maxrd+0.05];
 width_f2ax = f2ax.Position(3);
 height_f2ax = f2ax.Position(4);
 end_point_x = ((arrow_end_cor(1)-f2ax.XLim(1))/(f2ax.XLim(2)-f2ax.XLim(1)))*width_f2ax + f2ax.Position(1);
 end_point_y = ((arrow_end_cor(2)-f2ax.YLim(1))/(f2ax.YLim(2)-f2ax.YLim(1)))*height_f2ax + f2ax.Position(2);
 start_point_x = ((arrow_start_cor(1)-f2ax.XLim(1))/(f2ax.XLim(2)-f2ax.XLim(1)))*width_f2ax + f2ax.Position(1);
 start_point_y = ((arrow_start_cor(2)-f2ax.YLim(1))/(f2ax.YLim(2)-f2ax.YLim(1)))*height_f2ax + f2ax.Position(2);
 % text
f2text = ['$r^{*}=',num2str(maxrd),'$' newline '$C_r^{*}=' num2str(target_ratio) '$'];  
% draw arrow
annotation('arrow',[ start_point_x,end_point_x],[start_point_y,end_point_y]);           
% add text
annotation('textbox',[ start_point_x,start_point_y,0.1,0.1],'String',f2text,'LineStyle','none','Interpreter','latex','FitBoxToText','on'); 

%% sim2 relative degree zoom-in chart
f3 = figure(3);
zoom_in_range = 20;
size_max_rd = size(num_of_max_rd,1);
% left
yyaxis left
stairs(combval(num_of_max_rd(1)-zoom_in_range:num_of_max_rd(size_max_rd)+zoom_in_range),rdval(num_of_max_rd(1)-zoom_in_range:num_of_max_rd(size_max_rd)+zoom_in_range),...
    '-','LineWidth',1.5)       
ylabel('Relative Degree')
xlabel('Ratio($C_r$)')
hold on
plot(target_ratio,maxrd,'ro','MarkerSize',6)
hold off
% right
yyaxis right
bar(combval(num_of_max_rd),ptgcfdt)
ylabel('Relative Data Usage')
if (TitleFlag)
    title('Relative Degree and Data Usage Zoomed In');
end
if (GridFlag)
    grid on
end
% axes
f3ax = f3.CurrentAxes;
f3ax.YAxis(1).Limits=[0 4];
f3ax.YAxis(2).Limits=[0 0.04];
f3ax.XAxis.Limits=[-0.0167 0.018];
%legend
legend('Relative Degree','Optimal point','Relative Data Usage','Orientation','vertical')
% arrow position
 arrow_start_cor = [target_ratio-0.002 maxrd+0.3];   
 arrow_end_cor =[target_ratio maxrd+0.08];
 width_f3ax = f3ax.Position(3);
 height_f3ax = f3ax.Position(4);
 end_point_x = (arrow_end_cor(1)-f3ax.XLim(1))/(f3ax.XLim(2)-f3ax.XLim(1))*width_f3ax + f3ax.Position(1);
 end_point_y = ((arrow_end_cor(2)-f3ax.YAxis(1).Limits(1))/(f3ax.YAxis(1).Limits(2)-f3ax.YAxis(1).Limits(1)))*height_f3ax + f3ax.Position(2);
 start_point_x = ((arrow_start_cor(1)-f3ax.XLim(1))/(f3ax.XLim(2)-f3ax.XLim(1)))*width_f3ax + f3ax.Position(1);
 start_point_y = ((arrow_start_cor(2)-f3ax.YAxis(1).Limits(1))/(f3ax.YAxis(1).Limits(2)-f3ax.YAxis(1).Limits(1)))*height_f3ax + f3ax.Position(2);
 % text
f2text = ['$r^{*}=',num2str(maxrd),'$' newline '$C_r^{*}=' num2str(target_ratio) '$'];  
% draw arrow
annotation('arrow',[ start_point_x,end_point_x],[start_point_y,end_point_y]);           
% add text
annotation('textbox',[ start_point_x-0.2,start_point_y+0.05,0.1,0.1],'String',f2text,'LineStyle','none','Interpreter','latex','FitBoxToText','on'); 