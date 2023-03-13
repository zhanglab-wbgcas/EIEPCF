function roc_data = ROC(graph_data, flag)
    % Method:    
    % Draw ROC(Receiver Operating Characteristics,受试者工作特征) curve.
	
	% Input:
    % graph_data:一个包含画图数据的二维矩阵，总共两列。第一列为金标准，1代表有调控关系，0代表无调控关系；第二列为与金标准相对应的预测的调控关系
    % flag:是否画ROC曲线图的标志值，1代表画图，0代表不画图，默认值为1
   
    % Output:
    % roc_data:一个结构体，记录了画ROC曲线时自变量X的所有取值点、与自变量X相对应的所有函数值Y,以及ROC曲线下的面积

    % Author: PengHuixiang
    % Version data: 2022-03-23

    
    %% Check input parameters

	% 输入的参数最少1个，最多2个，否则报错
    narginchk(1,2);
    
    % 检查输入的画图数据矩阵的值是否均为有限的数值
    if ~all(isnumeric(graph_data(:))) || ~all(isfinite(graph_data(:)))
        error('Warning: All of the values in graph_data must be numeric and finite')
    end
   
    % 检查输入的画图数据矩阵的金标准列是否包含所有的相互作用边和非相互作用边
    check_gold = logical(graph_data(:,1));
    if all(check_gold==1)
        error('Warning: The gold_standard has only interactive edges!')
    end
    if all(check_gold==0)
        error('Warning: The gold_standard has only non-interactive edges!')
    end
	
    % 检查输入的画图标志值flag是否为1或0
    if nargin > 2 && (~isa(flag,'numeric') || (isa(flag,'numeric') && ~ismember(flag,[0,1])))
        error('Input parameter flag must be 1 or 0.');
    end
    
    
	%% Set and display parameters
    if nargin < 2
        flag = 1;
    end
    
    fprintf('Size of the graph data matrix: (%d,%d) \n',size(graph_data,1),size(graph_data,2));
    fprintf('Flag value of whether to draw ROC curve: %d \n\n',flag);
	
    
    %% Get the ROC curve plotting coordinates and the area under the curve
    
    % 将画图数据矩阵按推断出来的调控关系列的调控强度升序排列
    sorted_preValues = sortrows(graph_data,2);
    % 获得画图数据矩阵调控关系列中的唯一值
    unique_preValues = unique(sorted_preValues(:,2));
    % 获得唯一值的数量
    uniValues_number = length(unique_preValues);
    
    % 定义一个二维矩阵，用于记录不同时刻的假阳性率和真阳性率
    FPR_and_TPR = zeros(uniValues_number,2);
   
    % 求不同时刻的假阳性率和真阳性率
    for i = 1:uniValues_number
        
        % 求此时的真阳性数、假阳性数、假阴性数和真阴性数
        TP = length(graph_data(graph_data(:,1)==1 & graph_data(:,2)>unique_preValues(i)));
        FP = length(graph_data(graph_data(:,1)==0 & graph_data(:,2)>unique_preValues(i)));
        FN = length(graph_data(graph_data(:,1)==1 & graph_data(:,2)<=unique_preValues(i)));
        TN = length(graph_data(graph_data(:,1)==0 & graph_data(:,2)<=unique_preValues(i)));
       
        % 记录此时的特异性和敏感性
        % FPR_and_TPR(i,:) = [TN/(TN+FP), TP/(TP+FN)];
        FPR_and_TPR(i,:) = [FP/(FP+TN), TP/(TP+FN)]; 
        
    end
    
    % 获得画ROC曲线的坐标点
    % xroc = flipud([1; 1-FPR_and_TPR(:,1); 0]);
    xroc = flipud([1; FPR_and_TPR(:,1); 0]);
    yroc = flipud([1; FPR_and_TPR(:,2); 0]);

	% 采用梯形法求积分，即ROC曲线下面积
    auc = trapz(xroc, yroc);
    
    % 如果调用此方法时接收返回值，则给返回值赋值，否则不给返回值赋值
    if nargout
        roc_data.xroc = xroc;
        roc_data.yroc = yroc;
        roc_data.auc = auc;
    end
    
	if flag
		% 绘制ROC曲线图
		plot(xroc,yroc,'r.-');
		hold off;
		xlabel('False positive rate (1-Specificity)');
		ylabel('True positive rate (Sensitivity)');
		title('ROC curve');
		axis square;
	end
    
    
end

