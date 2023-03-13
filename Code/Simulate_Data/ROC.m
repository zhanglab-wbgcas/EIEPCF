function roc_data = ROC(graph_data, flag)
    % Method:    
    % Draw ROC(Receiver Operating Characteristics,�����߹�������) curve.
	
	% Input:
    % graph_data:һ��������ͼ���ݵĶ�ά�����ܹ����С���һ��Ϊ���׼��1�����е��ع�ϵ��0�����޵��ع�ϵ���ڶ���Ϊ����׼���Ӧ��Ԥ��ĵ��ع�ϵ
    % flag:�Ƿ�ROC����ͼ�ı�־ֵ��1����ͼ��0������ͼ��Ĭ��ֵΪ1
   
    % Output:
    % roc_data:һ���ṹ�壬��¼�˻�ROC����ʱ�Ա���X������ȡֵ�㡢���Ա���X���Ӧ�����к���ֵY,�Լ�ROC�����µ����

    % Author: PengHuixiang
    % Version data: 2022-03-23

    
    %% Check input parameters

	% ����Ĳ�������1�������2�������򱨴�
    narginchk(1,2);
    
    % �������Ļ�ͼ���ݾ����ֵ�Ƿ��Ϊ���޵���ֵ
    if ~all(isnumeric(graph_data(:))) || ~all(isfinite(graph_data(:)))
        error('Warning: All of the values in graph_data must be numeric and finite')
    end
   
    % �������Ļ�ͼ���ݾ���Ľ��׼���Ƿ�������е��໥���ñߺͷ��໥���ñ�
    check_gold = logical(graph_data(:,1));
    if all(check_gold==1)
        error('Warning: The gold_standard has only interactive edges!')
    end
    if all(check_gold==0)
        error('Warning: The gold_standard has only non-interactive edges!')
    end
	
    % �������Ļ�ͼ��־ֵflag�Ƿ�Ϊ1��0
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
    
    % ����ͼ���ݾ����ƶϳ����ĵ��ع�ϵ�еĵ���ǿ����������
    sorted_preValues = sortrows(graph_data,2);
    % ��û�ͼ���ݾ�����ع�ϵ���е�Ψһֵ
    unique_preValues = unique(sorted_preValues(:,2));
    % ���Ψһֵ������
    uniValues_number = length(unique_preValues);
    
    % ����һ����ά�������ڼ�¼��ͬʱ�̵ļ������ʺ���������
    FPR_and_TPR = zeros(uniValues_number,2);
   
    % ��ͬʱ�̵ļ������ʺ���������
    for i = 1:uniValues_number
        
        % ���ʱ����������������������������������������
        TP = length(graph_data(graph_data(:,1)==1 & graph_data(:,2)>unique_preValues(i)));
        FP = length(graph_data(graph_data(:,1)==0 & graph_data(:,2)>unique_preValues(i)));
        FN = length(graph_data(graph_data(:,1)==1 & graph_data(:,2)<=unique_preValues(i)));
        TN = length(graph_data(graph_data(:,1)==0 & graph_data(:,2)<=unique_preValues(i)));
       
        % ��¼��ʱ�������Ժ�������
        % FPR_and_TPR(i,:) = [TN/(TN+FP), TP/(TP+FN)];
        FPR_and_TPR(i,:) = [FP/(FP+TN), TP/(TP+FN)]; 
        
    end
    
    % ��û�ROC���ߵ������
    % xroc = flipud([1; 1-FPR_and_TPR(:,1); 0]);
    xroc = flipud([1; FPR_and_TPR(:,1); 0]);
    yroc = flipud([1; FPR_and_TPR(:,2); 0]);

	% �������η�����֣���ROC���������
    auc = trapz(xroc, yroc);
    
    % ������ô˷���ʱ���շ���ֵ���������ֵ��ֵ�����򲻸�����ֵ��ֵ
    if nargout
        roc_data.xroc = xroc;
        roc_data.yroc = yroc;
        roc_data.auc = auc;
    end
    
	if flag
		% ����ROC����ͼ
		plot(xroc,yroc,'r.-');
		hold off;
		xlabel('False positive rate (1-Specificity)');
		ylabel('True positive rate (Sensitivity)');
		title('ROC curve');
		axis square;
	end
    
    
end

