function [target_ratio]=plot_result(rdval,ptgval,combval)
maxrd = max(rdval);             % max rd
minrd = min(rdval);             
num_of_max_rd = find(rdval==maxrd); % possible ratio for max rd
if size(combval,1) ~= 1         % manul mode
    if (maxrd~=minrd)&&(size(num_of_max_rd,1)>1) 
        
%        confidence rate
        ptgcfdt= ptgval(num_of_max_rd)./sum(ptgval(num_of_max_rd));
        idx_of_max_ptgcfdt = find(ptgcfdt==max(ptgcfdt));
         
        target_list = combval(num_of_max_rd(1) + idx_of_max_ptgcfdt - 1);
        [~,target_idx] = min(abs(target_list));
        target_ratio = target_list(target_idx);
%         figure(2)
%         yyaxis left
%         plot(combval,rdval)
%         ylabel('Relative Degree')
%         ylim([0 5])
%         xlabel('Ratio')
%         yyaxis right
%         bar(combval(num_of_max_rd),ptgval(num_of_max_rd))
%         ylim([0 0.3])
%         ylabel('Confidence Rate')
%         %title('Relative Degree and Confidence Rate')
%         save()
%         figure(3)
%         yyaxis left
%         plot(combval(num_of_max_rd(1)-2:num_of_max_rd(end)),rdval(num_of_max_rd(1)-2:num_of_max_rd(end)))
% 
%         ylabel('Relative Degree')
%         xlabel('Ratio')
%         yyaxis right
%         bar(combval(num_of_max_rd),ptgval(num_of_max_rd))
%         ylabel('Confidence Rate')
        %title('Relative Degree and Confidence Rate(Zoomed in)')
       
        %target_ratio
      % maxrd
        
    else                                         % if only one max rd 
%         figure
%         plot(combval,rdval)
%         ylabel('Relative Degree')
%         xlabel('Ratio')
%         title('Relative Degree')
% %         
         target_ratio = combval(num_of_max_rd);
       % maxrd
        
    end
end