%run vb_CCDT_connectomic_plots first

subj = 'HUP154';
saveon = 1;
ploton = 0;
subploton = 1;
cv=1;
cvMethod = '10-fold'; %'10-fold' vs. 'N-1'

%% Classifier Training Set-up
if ~cv
    cvMethod = 'ModelFit';
end

all500Trials = (vRT>0&vDT<550);
all1500Trials = (vRT>0&vDT>1450);
trials(1,:) = all1500Trials;
trials(2,:) = ifast1500;
trials(3,:) = islow1500;
trials(4,:) = all500Trials;
trials(5,:) = ifast500;
trials(6,:) = islow500;
j=0;
k=0;
for i=1:length(vRT)
    if trials(1,i)==1
        j=j+1;
        if trials(1,i)+trials(2,i)==2
            ccFast1500(j,1)=1;
        else
            ccFast1500(j,1)=0;
        end
        if trials(1,i)+trials(3,i)==2
            ccSlow1500(j,1)=1;
        else
            ccSlow1500(j,1)=0;
        end
    elseif trials(4,i)==1
        k=k+1;
        if trials(4,i)+trials(5,i)==2
            ccFast500(k,1)=1;
        else
            ccFast500(k,1)=0;
        end
        if trials(4,i)+trials(6,i)==2
            ccSlow500(k,1)=1;
        else
            ccSlow500(k,1)=0;
        end
    end
    
end

ccMid500 = (~ccFast500&~ccSlow500);
ccSlow500 = logical(ccSlow500);
ccFast500 = logical(ccFast500);
ccMid500 = logical(ccMid500);

ccMid1500 = (~ccFast1500&~ccSlow1500);
ccSlow1500 = logical(ccSlow1500);
ccFast1500 = logical(ccFast1500);

ccFastAvg1500 = (vRT(all1500Trials)<median(vRT(all1500Trials)));
ccSlowAvg1500 = (vRT(all1500Trials)>=median(vRT(all1500Trials)));

ccFastAvg500 = (vRT(all500Trials)<median(vRT(all500Trials)));
ccSlowAvg500 = (vRT(all500Trials)>=median(vRT(all500Trials)));



%%
ccTrialsRT_500 = vRT(all500Trials);
ccTrialsRT_500_Bi = ccTrialsRT_500;
ccTrialsRT_500_Bi(ccMid500) = [];
ccTrialsStr_500_HF = mean(str_precue(all500Trials,:,2),2);
ccTrialsStr_500_HF_Bi = ccTrialsStr_500_HF;
ccTrialsStr_500_HF_Bi(ccMid500) = [];
ccTrialsQexp_500_HF = mean(nodalQexp(all500Trials,:,2),2);
ccTrialsQexp_500_HF_Bi = ccTrialsQexp_500_HF;
ccTrialsQexp_500_HF_Bi(ccMid500) = [];
ccTrialsStr_500_LF = mean(str_precue(all500Trials,:,1),2);
ccTrialsStr_500_LF_Bi = ccTrialsStr_500_LF;
ccTrialsStr_500_LF_Bi(ccMid500) = [];
ccTrialsQexp_500_LF = mean(nodalQexp(all500Trials,:,1),2);
ccTrialsQexp_500_LF_Bi = ccTrialsQexp_500_LF;
ccTrialsQexp_500_LF_Bi(ccMid500) = [];
ccTrialsStr_500_LG = mean(str_precue(all500Trials,:,3),2);
ccTrialsStr_500_LG_Bi = ccTrialsStr_500_LG;
ccTrialsStr_500_LG_Bi(ccMid500) = [];
ccTrialsQexp_500_LG = mean(nodalQexp(all500Trials,:,3),2);
ccTrialsQexp_500_LG_Bi = ccTrialsQexp_500_LF;
ccTrialsQexp_500_LG_Bi(ccMid500) = [];

classNames500 = cell(sum(all500Trials),1);
classNames500(ccSlow500) = {'SlowRT'};
classNames500(ccFast500) = {'FastRT'};
classNames500(ccMid500) = {'MidRT'};
classNames500Binary = classNames500;
classNames500Binary(ccMid500) = [];

ccTrialsRT_1500 = vRT(all1500Trials);
ccTrialsRT_1500_Bi = ccTrialsRT_1500;
ccTrialsRT_1500_Bi(ccMid1500) = [];
ccTrialsStr_1500_HF = mean(str_precue(all1500Trials,:,2),2);
ccTrialsStr_1500_HF_Bi = ccTrialsStr_1500_HF;
ccTrialsStr_1500_HF_Bi(ccMid1500) = [];
ccTrialsQexp_1500_HF = mean(nodalQexp(all1500Trials,:,2),2);
ccTrialsQexp_1500_HF_Bi = ccTrialsQexp_1500_HF;
ccTrialsQexp_1500_HF_Bi(ccMid1500) = [];
ccTrialsStr_1500_LF = mean(str_precue(all1500Trials,:,1),2);
ccTrialsStr_1500_LF_Bi = ccTrialsStr_1500_LF;
ccTrialsStr_1500_LF_Bi(ccMid1500) = [];
ccTrialsQexp_1500_LF = mean(nodalQexp(all1500Trials,:,1),2);
ccTrialsQexp_1500_LF_Bi = ccTrialsQexp_1500_LF;
ccTrialsQexp_1500_LF_Bi(ccMid1500) = [];
ccTrialsStr_1500_LG = mean(str_precue(all1500Trials,:,3),2);
ccTrialsStr_1500_LG_Bi = ccTrialsStr_1500_LG;
ccTrialsStr_1500_LG_Bi(ccMid1500) = [];
ccTrialsQexp_1500_LG = mean(nodalQexp(all1500Trials,:,3),2);
ccTrialsQexp_1500_LG_Bi = ccTrialsQexp_1500_LG;
ccTrialsQexp_1500_LG_Bi(ccMid1500) = [];


classNames1500 = cell(sum(all1500Trials),1);
classNames1500(ccSlow1500) = {'SlowRT'};
classNames1500(ccFast1500) = {'FastRT'};
classNames1500(ccMid1500) = {'MidRT'};
classNames1500Binary = classNames1500;
classNames1500Binary(ccMid1500) = [];

% Use above or below average instead of fastest and slowest third
cN_1500 = cell(sum(all1500Trials),1);
cN_1500(ccFastAvg1500) = {'FastRT'};
cN_1500(ccSlowAvg1500) = {'SlowRT'};

cN_500 = cell(sum(all500Trials),1);
cN_500(ccFastAvg500) = {'FastRT'};
cN_500(ccSlowAvg500) = {'SlowRT'};

%% Classifier training
%Above vs. below avg
sixFdat_1500 = [ccTrialsStr_1500_HF_Bi, ccTrialsQexp_1500_HF_Bi, ccTrialsStr_1500_LF_Bi, ccTrialsQexp_1500_LF_Bi, ccTrialsStr_1500_LG_Bi, ccTrialsQexp_1500_LG_Bi];
sixFdat_500 = [ccTrialsStr_500_HF_Bi, ccTrialsQexp_500_HF_Bi, ccTrialsStr_500_LF_Bi, ccTrialsQexp_500_LF_Bi, ccTrialsStr_500_LG_Bi, ccTrialsQexp_500_LG_Bi];
sixFdat_1500cont = [ccTrialsStr_1500_HF, ccTrialsQexp_1500_HF, ccTrialsStr_1500_LF, ccTrialsQexp_1500_LF, ccTrialsStr_1500_LG, ccTrialsQexp_1500_LG];
sixFdat_500cont = [ccTrialsStr_500_HF, ccTrialsQexp_500_HF, ccTrialsStr_500_LF, ccTrialsQexp_500_LF, ccTrialsStr_500_LG, ccTrialsQexp_500_LG];

fourFdat_1500 = [ccTrialsStr_1500_HF_Bi, ccTrialsQexp_1500_HF_Bi, ccTrialsStr_1500_LF_Bi, ccTrialsQexp_1500_LF_Bi];
fourFdat_500 = [ccTrialsStr_500_HF_Bi, ccTrialsQexp_500_HF_Bi, ccTrialsStr_500_LF_Bi, ccTrialsQexp_500_LF_Bi];
fourFdat_1500cont = [ccTrialsStr_1500_HF, ccTrialsQexp_1500_HF, ccTrialsStr_1500_LF, ccTrialsQexp_1500_LF];
fourFdat_500cont = [ccTrialsStr_500_HF, ccTrialsQexp_500_HF, ccTrialsStr_500_LF, ccTrialsQexp_500_LF];

twoFdat_1500 = [ccTrialsStr_1500_HF_Bi, ccTrialsQexp_1500_HF_Bi];
twoFdat_500 = [ccTrialsStr_500_HF_Bi, ccTrialsQexp_500_HF_Bi];
twoFdat_1500cont = [ccTrialsStr_1500_HF, ccTrialsQexp_1500_HF];
twoFdat_500cont = [ccTrialsStr_500_HF, ccTrialsQexp_500_HF];

classifierStruct = struct;

for xi = 1:3
    switch xi
        case 1
            currFdat_1500 = twoFdat_1500;
            currFdat_500 = twoFdat_500;
            currFdat_1500cont = twoFdat_1500cont;
            currFdat_500cont = twoFdat_500cont;
            fS = {'Str_HF, Qexp_HF'};
            
        case 2
            currFdat_1500 = fourFdat_1500;
            currFdat_500 = fourFdat_500;
            currFdat_1500cont = fourFdat_1500cont;
            currFdat_500cont = fourFdat_500cont;
            fS = {'Str_HF, Qexp_HF, Str_LF, Qexp_LF'};
            
        case 3
            
            currFdat_1500 = sixFdat_1500;
            currFdat_500 = sixFdat_500;
            currFdat_1500cont = sixFdat_1500cont;
            currFdat_500cont = sixFdat_500cont;
            fS = {'Str_HF, Qexp_HF, Str_LF, Qexp_LF, Str_LG, Qexp_LG'};
            
    end
    
    for xii = 1:2
        switch xii
            case 1
                currClassNames_1500 = classNames1500Binary;
                currClassNames_500 = classNames500Binary;
                cFdat_1500 = currFdat_1500;
                cFdat_500 = currFdat_500;
                cRT_1500 = ccTrialsRT_1500_Bi;
                cRT_500 = ccTrialsRT_500_Bi;
                cutoff = 'third';
                
                
            case 2
                currClassNames_1500 = cN_1500;
                currClassNames_500 = cN_500;
                cFdat_1500 = currFdat_1500cont;
                cFdat_500 = currFdat_500cont;
                cRT_1500 = ccTrialsRT_1500;
                cRT_500 = ccTrialsRT_500;
                cutoff = 'half';
                
                
        end
        
        resp_1500 = strcmp(currClassNames_1500,'FastRT');
        resp_500 = strcmp(currClassNames_500,'FastRT');
        
        %SVM
        mdlSVM = fitcsvm(cFdat_1500,currClassNames_1500);
        if ~cv
            [~,score_svm] = resubPredict(mdlSVM);
        elseif cv&&strcmpi(cvMethod,'10-fold')
            mdlSVMcv = fitPosterior(mdlSVM, 'kfold', 10);
            [~,score_svm] = resubPredict(mdlSVMcv);
        elseif cv&&strcmpi(cvMethod,'N-1')
            mdlSVMcv = fitPosterior(mdlSVM, 'leaveout', 'on');
            [~,score_svm] = resubPredict(mdlSVMcv);
        end
        [Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(currClassNames_1500,score_svm(:,1),'FastRT', 'NBoot', 1000, 'boottype', 'cper');
        
        mdlSVM_500 = fitcsvm(cFdat_500,currClassNames_500);
        if ~cv
            [~,score_svm_500] = resubPredict(mdlSVM_500);
        elseif cv&&strcmpi(cvMethod,'10-fold')
            mdlSVM_500cv = fitPosterior(mdlSVM_500, 'kfold', 10);
            [~,score_svm_500] = resubPredict(mdlSVM_500cv);
        elseif cv&&strcmpi(cvMethod,'N-1')
            mdlSVM_500cv = fitPosterior(mdlSVM_500, 'leaveout', 'on');
            [~,score_svm_500] = resubPredict(mdlSVM_500cv);
        end
        [Xsvm_500,Ysvm_500,Tsvm_500,AUCsvm_500] = perfcurve(currClassNames_500,score_svm_500(:,1),'FastRT', 'NBoot', 1000, 'boottype', 'cper');
        
        %Logistic regression
        
        mdlGLM = fitglm(cFdat_1500,resp_1500, 'distribution', 'binomial');
        if ~cv
            score_glm = mdlGLM.Fitted.Probability;
        elseif cv&&strcmpi(cvMethod,'10-fold')
            mdlGLMcv = fitclinear(cFdat_1500,currClassNames_1500,'learner','logistic','kfold',10);
            [~,score_glm] = kfoldPredict(mdlGLMcv);
        elseif cv&&strcmpi(cvMethod, 'N-1')
            mdlGLMcv = fitclinear(cFdat_1500,currClassNames_1500,'learner','logistic','leaveout','on');
            [~,score_glm] = kfoldPredict(mdlGLMcv);
        end
        [Xglm,Yglm,Tglm,AUCglm] = perfcurve(currClassNames_1500,score_glm(:,1),'FastRT', 'NBoot', 1000, 'boottype', 'cper');
        %
        %     [B,dev,stats] = glmfit(cFdat_1500,resp_1500,'binomial');
        %     [Bc,devc,statsc] = glmfit(currFdat_1500cont,vRT(all1500Trials));
        %
        mdlGLM_500 = fitglm(cFdat_500,resp_500, 'distribution', 'binomial');
        if ~cv
            score_glm_500 = mdlGLM_500.Fitted.Probability;
        elseif cv&&strcmpi(cvMethod,'10-fold')
            mdlGLM_500cv = fitclinear(cFdat_500,currClassNames_500,'learner','logistic','kfold',10);
            [~,score_glm_500] = kfoldPredict(mdlGLM_500cv);
        elseif cv&&strcmpi(cvMethod, 'N-1')
            mdlGLM_500cv = fitclinear(cFdat_500,currClassNames_500,'learner','logistic','leaveout','on');
            [~,score_glm_500] = kfoldPredict(mdlGLM_500cv);
        end
        [Xglm_500,Yglm_500,Tglm_500,AUCglm_500] = perfcurve(currClassNames_500,score_glm_500(:,1),'FastRT', 'NBoot', 1000, 'boottype', 'cper');
        
        %     [B_500,dev_500,stats_500] = glmfit(currFdat_500,resp_500,'binomial');
        %     [B_500c,dev_500c,stats_500c] = glmfit(currFdat_500cont,vRT(all500Trials));
        %
        
        %Decision Tree
        mdlTree = fitctree(cFdat_1500,currClassNames_1500);
        if ~cv
            [~,score_tree] = resubPredict(mdlTree);
        elseif cv&&strcmpi(cvMethod,'10-fold')
            mdlTreecv = crossval(mdlTree, 'kfold', 10);
            [~,score_tree] = kfoldPredict(mdlTreecv);
        elseif cv&&strcmpi(cvMethod,'N-1')
            mdlTreecv = crossval(mdlTree, 'leaveout', 'on');
            [~,score_tree] = kfoldPredict(mdlTreecv);
        end
        [Xtree,Ytree,Ttree,AUCtree] = perfcurve(currClassNames_1500,score_tree(:,1),'FastRT', 'NBoot', 1000, 'boottype', 'cper');
        
        mdlTree_500 = fitctree(cFdat_500,currClassNames_500);
        if ~cv
            [~,score_tree_500] = resubPredict(mdlTree_500);
        elseif cv&&strcmpi(cvMethod,'10-fold')
            mdlTree_500cv = crossval(mdlTree_500, 'kfold', 10);
            [~,score_tree_500] = kfoldPredict(mdlTree_500cv);
        elseif cv&&strcmpi(cvMethod,'N-1')
            mdlTree_500cv = crossval(mdlTree_500, 'leaveout', 'on');
            [~,score_tree_500] = kfoldPredict(mdlTree_500cv);
        end
        [Xtree_500,Ytree_500,Ttree_500,AUCtree_500] = perfcurve(currClassNames_500,score_tree_500(:,1),'FastRT', 'NBoot', 1000, 'boottype', 'cper');
        
        % Decision tree for multiclass. Can also use fitcecoc for SVM based
        % multiclass model support
        % mdlTree_MC = fitctree([ccTrialsStr_HF, ccTrialsQexp_HF, ccTrialsStr_LF],classNames1500);
        % [~,score_treeMC] = resubPredict(mdlTree_MC);
        % diffscore = score_treeMC(:,1) - max(score_treeMC(:,2),score_treeMC(:,3));
        % [XtreeMC,YtreeMC,TtreeMC,~,OPTROCPTtreeMC,subytreeMC,subnamestreeMC] = perfcurve(classNames1500,diffscore,'FastRT');
        
        %% Plots
        if ploton
            h1=figure;
            plot(Xsvm,Ysvm)
            hold on
            plot(Xglm,Yglm, 'color', [0.91 0.41 0.17])
            plot(Xtree,Ytree, 'color', [0.5 0 0.9])
            plot([0 1],[0 1], 'k--')
            box off
            xlabel('False Positive Rate')
            ylabel('True Positive Rate')
            suptitle(['ROC for Binary Classifiers: 1500ms Trials; Feature Space =  ' num2str(xi*2)]);
            if cv
                title(['fastest vs slowest: ' cutoff, '; CV =  ', cvMethod]);
            else
                title(['fastest vs slowest: ' cutoff, '; ModelFit']);
            end
            legend(['SVM (auc=' num2str(AUCsvm, '%4.2f') ')'], ['Logit GLM (auc=' num2str(AUCglm, '%4.2f') ')'], ['Decision Tree (auc=' num2str(AUCtree, '%4.2f') ')'], 'location', 'best')
            set(gcf,'color','w')
            
            h2=figure;
            plot(Xsvm_500,Ysvm_500)
            hold on
            plot(Xglm_500,Yglm_500, 'color', [0.91 0.41 0.17])
            plot(Xtree_500,Ytree_500, 'color', [0.5 0 0.9])
            plot([0 1],[0 1],'k--')
            box off
            xlabel('False Positive Rate')
            ylabel('True Positive Rate')
            suptitle(['ROC for Binary Classifiers: 500ms Trials; Feature Space =  ' num2str(xi*2)]);
            if cv
                title(['fastest vs slowest: ' cutoff, '; CV =  ', cvMethod]);
            else
                title(['fastest vs slowest: ' cutoff, '; ModelFit']);
            end
            legend(['SVM (auc=' num2str(AUCsvm_500, '%4.2f') ')'], ['Logit GLM (auc=' num2str(AUCglm_500, '%4.2f') ')'], ['Decision Tree (auc=' num2str(AUCtree_500, '%4.2f') ')'], 'location', 'best')
            set(gcf,'color','w')
        end
        
        tbl_glm = table(score_glm, currClassNames_1500, cRT_1500, 'VariableName', {'FastRT_prob', 'True_Class', 'Observed_RT'});
        %     writetable(tbl_glm, sprintf('%s_LogisticRegBinary.xls',subj));
        tbl_glm_500 = table(score_glm_500, currClassNames_500, cRT_500, 'VariableName', {'FastRT_prob', 'True_Class', 'Observed_RT'});
        %     writetable(tbl_glm, sprintf('%s_LogisticRegBinary_500.xls',subj));
        
        tbl_svm = table(score_svm(:,1), score_svm(:,2), currClassNames_1500, cRT_1500, 'VariableName', {'FastRT_prob', 'SlowRT_prob', 'True_Class', 'Observed_RT'});
        %     writetable(tbl_tree, sprintf('%s_DecisionTreeBinary.xls',subj));
        tbl_svm_500 = table(score_svm_500(:,1), score_svm_500(:,2), currClassNames_500, cRT_500, 'VariableName', {'FastRT_prob', 'SlowRT_prob', 'True_Class', 'Observed_RT'});
        %     writetable(tbl_tree_500, sprintf('%s_DecisionTreeBinary_500.xls',subj));
        
        
        tbl_tree = table(score_tree(:,1), score_tree(:,2), currClassNames_1500, cRT_1500, 'VariableName', {'FastRT_prob', 'SlowRT_prob', 'True_Class', 'Observed_RT'});
        %     writetable(tbl_tree, sprintf('%s_DecisionTreeBinary.xls',subj));
        tbl_tree_500 = table(score_tree_500(:,1), score_tree_500(:,2), currClassNames_500, cRT_500, 'VariableName', {'FastRT_prob', 'SlowRT_prob', 'True_Class', 'Observed_RT'});
        %     writetable(tbl_tree_500, sprintf('%s_DecisionTreeBinary_500.xls',subj));
        
        % create classifierStruct
        classifierStruct(xi).features = fS;
        classifierStruct(xi).cutoff(xii).cutoff = cutoff;
        if ~cv
            classifierStruct(xi).glm_1500(xii).model = mdlGLM;
            classifierStruct(xi).svm_1500(xii).model = mdlSVM;
            classifierStruct(xi).tree_1500(xii).model = mdlTree;
            classifierStruct(xi).glm_500(xii).model = mdlGLM_500;
            classifierStruct(xi).svm_500(xii).model = mdlSVM_500;
            classifierStruct(xi).tree_500(xii).model = mdlTree_500;
            classifierStruct(xi).cvMethod = 'None';
        else
            classifierStruct(xi).cvMethod = cvMethod;
            classifierStruct(xi).glm_1500(xii).model = mdlGLMcv;
            classifierStruct(xi).svm_1500(xii).model = mdlSVMcv;
            classifierStruct(xi).tree_1500(xii).model = mdlTreecv;
            classifierStruct(xi).glm_500(xii).model = mdlGLM_500cv;
            classifierStruct(xi).svm_500(xii).model = mdlSVM_500cv;
            classifierStruct(xi).tree_500(xii).model = mdlTree_500cv;
        end
        classifierStruct(xi).glm_1500(xii).score = score_glm;
        classifierStruct(xi).svm_1500(xii).score= score_svm;
        classifierStruct(xi).tree_1500(xii).score = score_tree;
        classifierStruct(xi).glm_1500(xii).perf.auc = AUCglm;
        classifierStruct(xi).svm_1500(xii).perf.auc = AUCsvm;
        classifierStruct(xi).tree_1500(xii).perf.auc = AUCtree;
        classifierStruct(xi).glm_1500(xii).perf.X = Xglm;
        classifierStruct(xi).svm_1500(xii).perf.X = Xsvm;
        classifierStruct(xi).tree_1500(xii).perf.X = Xtree;
        classifierStruct(xi).glm_1500(xii).perf.Y = Yglm;
        classifierStruct(xi).svm_1500(xii).perf.Y = Ysvm;
        classifierStruct(xi).tree_1500(xii).perf.Y = Ytree;
        classifierStruct(xi).glm_1500(xii).table = sortrows(tbl_glm,1,'descend');
        classifierStruct(xi).svm_1500(xii).table = sortrows(tbl_svm,1,'descend');
        classifierStruct(xi).tree_1500(xii).table = sortrows(tbl_tree,1,'descend');
        
        classifierStruct(xi).glm_500(xii).score = score_glm_500;
        classifierStruct(xi).svm_500(xii).score = score_svm_500;
        classifierStruct(xi).tree_500(xii).score = score_tree_500;
        classifierStruct(xi).glm_500(xii).perf.auc = AUCglm_500;
        classifierStruct(xi).svm_500(xii).perf.auc = AUCsvm_500;
        classifierStruct(xi).tree_500(xii).perf.auc = AUCtree_500;
        classifierStruct(xi).glm_500(xii).perf.X = Xglm_500;
        classifierStruct(xi).svm_500(xii).perf.X = Xsvm_500;
        classifierStruct(xi).tree_500(xii).perf.X = Xtree_500;
        classifierStruct(xi).glm_500(xii).perf.Y = Yglm_500;
        classifierStruct(xi).svm_500(xii).perf.Y = Ysvm_500;
        classifierStruct(xi).tree_500(xii).perf.Y = Ytree_500;
        classifierStruct(xi).glm_500(xii).table = sortrows(tbl_glm_500,1,'descend');
        classifierStruct(xi).svm_500(xii).table = sortrows(tbl_svm_500,1,'descend');
        classifierStruct(xi).tree_500(xii).table = sortrows(tbl_tree_500,1,'descend');
        
        
        if saveon&&ploton
            disp('saving figures...')
            figure(h1);
            print([sprintf('%s_',subj) 'classifierROC.ps'],'-dpsc2','-append');
            figure(h2);
            print([sprintf('%s_',subj) 'classifierROC.ps'],'-dpsc2','-append');
        end
        clear mdl* score_* tbl_* X* Y* AUC* h1 h2 predict_*
    end
end

%% overall subplot
% feature space plot
% figure;
% plot(ccTrialsStr_1500_HF_Bi, ccTrialsQexp_1500_HF_Bi,'k.')
% hold on
% plot(ccTrialsStr_1500_HF_Bi(resp_1500), ccTrialsQexp_1500_HF_Bi(resp_1500), 'go')
% plot(ccTrialsStr_1500_HF_Bi(~resp_1500), ccTrialsQexp_1500_HF_Bi(~resp_1500), 'mo')
% box off
% set(gcf,'color','w')

if subploton
    h3=figure;
    xy=0;
    for ii=1:3
        for jj=1:2
            xy=xy+1;
            subplot(3,2,xy)
            hold on
%             if classifierStruct(ii).svm_1500(jj).perf.auc < 0.5
%                 plot(1-classifierStruct(ii).svm_1500(jj).perf.X,1-classifierStruct(ii).svm_1500(jj).perf.Y)
%             elseif classifierStruct(ii).glm_1500(jj).perf.auc < 0.5
%                 plot(1-classifierStruct(ii).glm_1500(jj).perf.X,1-classifierStruct(ii).glm_1500(jj).perf.Y, 'color', [0.91 0.41 0.17])
%             elseif classifierStruct(ii).tree_1500(jj).perf.auc < 0.5
%                 plot(1-classifierStruct(ii).tree_1500(jj).perf.X,1-classifierStruct(ii).tree_1500(jj).perf.Y, 'color', [0.5 0 0.9])
%             else
                plot(classifierStruct(ii).svm_1500(jj).perf.X(:,1),classifierStruct(ii).svm_1500(jj).perf.Y(:,1))
                plot(classifierStruct(ii).glm_1500(jj).perf.X(:,1),classifierStruct(ii).glm_1500(jj).perf.Y(:,1), 'color', [0.91 0.41 0.17])
                plot(classifierStruct(ii).tree_1500(jj).perf.X(:,1),classifierStruct(ii).tree_1500(jj).perf.Y(:,1), 'color', [0.5 0 0.9])
%             end
            plot([0 1],[0 1],'k--')
            box off
            xticks([0 1])
            yticks([0 1])
            title(['fastest vs. slowest: ' classifierStruct(ii).cutoff(jj).cutoff, '; featSpace =  ', num2str(ii*2)]);
            legend(['SVM (auc=' num2str(classifierStruct(ii).svm_1500(jj).perf.auc(1), '%4.2f') ')'], ['Logit GLM (auc=' num2str(classifierStruct(ii).glm_1500(jj).perf.auc(1), '%4.2f') ')'], ['Decision Tree (auc=' num2str(classifierStruct(ii).tree_1500(jj).perf.auc(1), '%4.2f') ')'], 'location', 'best')
            set(gcf,'color','w')
            if xy==6
                xlabel('False Positive Rate')
                ylabel('True Positive Rate')
            end
            
        end
    end
    if cv
        suptitle(['ROC for Binary Classifiers: 1500ms Trials; CV = ' cvMethod]);
    else
        suptitle('ROC for Binary Classifiers: 1500ms Trials; ModelFit');
    end
end


if subploton
    h4=figure;
    xy=0;
    for ii=1:3
        for jj=1:2
            xy=xy+1;
            subplot(3,2,xy)
            hold on
%             if classifierStruct(ii).svm_500(jj).perf.auc < 0.5
%                 plot(1-classifierStruct(ii).svm_500(jj).perf.X,1-classifierStruct(ii).svm_500(jj).perf.Y)
%             elseif classifierStruct(ii).glm_500(jj).perf.auc < 0.5
%                 plot(1-classifierStruct(ii).glm_500(jj).perf.X,1-classifierStruct(ii).glm_500(jj).perf.Y, 'color', [0.91 0.41 0.17])
%             elseif classifierStruct(ii).tree_500(jj).perf.auc < 0.5
%                 plot(1-classifierStruct(ii).tree_500(jj).perf.X,1-classifierStruct(ii).tree_500(jj).perf.Y, 'color', [0.5 0 0.9])
%             else
                plot(classifierStruct(ii).svm_500(jj).perf.X(:,1),classifierStruct(ii).svm_500(jj).perf.Y(:,1))
                plot(classifierStruct(ii).glm_500(jj).perf.X(:,1),classifierStruct(ii).glm_500(jj).perf.Y(:,1), 'color', [0.91 0.41 0.17])
                plot(classifierStruct(ii).tree_500(jj).perf.X(:,1),classifierStruct(ii).tree_500(jj).perf.Y(:,1), 'color', [0.5 0 0.9])
%             end
            plot([0 1],[0 1],'k--')
            box off
            xticks([0 1])
            yticks([0 1])
            title(['fastest vs. slowest: ' classifierStruct(ii).cutoff(jj).cutoff, '; featSpace =  ', num2str(ii*2)]);
            legend(['SVM (auc=' num2str(classifierStruct(ii).svm_500(jj).perf.auc(1), '%4.2f') ')'], ['Logit GLM (auc=' num2str(classifierStruct(ii).glm_500(jj).perf.auc(1), '%4.2f') ')'], ['Decision Tree (auc=' num2str(classifierStruct(ii).tree_500(jj).perf.auc(1), '%4.2f') ')'], 'location', 'best')
            set(gcf,'color','w')
            if xy==6
                xlabel('False Positive Rate')
                ylabel('True Positive Rate')
            end
            
        end
    end
    if cv
        suptitle(['ROC for Binary Classifiers: 500ms Trials; CV = ' cvMethod]);
    else
        suptitle('ROC for Binary Classifiers: 500ms Trials; ModelFit');
    end
end




%%
if saveon&&subploton
    disp('saving figures...')
    set(h3,'units','normalized','outerposition',[0 0 1 1])
    set(h4,'units','normalized','outerposition',[0 0 1 1])
    figure(h3);
    print([sprintf('%s_',subj) 'classifierROC_1500_' cvMethod '.eps'],'-depsc');
    print([sprintf('%s_',subj) 'classifierROC_1500_' cvMethod '.png'],'-dpng');
    figure(h4);
    print([sprintf('%s_',subj) 'classifierROC_500_' cvMethod '.eps'],'-depsc');
    print([sprintf('%s_',subj) 'classifierROC_500_' cvMethod '.png'],'-dpng');
end

if saveon
    disp('saving classifierStruct...');
    save([sprintf('%s',subj) '_classifierStruct_' cvMethod],'classifierStruct');
end