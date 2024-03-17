result=sim(network1,tr_in);
    for i=1:length(result)
        if result(1,i) > result(2,i)
            tr_r(:,i)=[1;-1];
        else
            tr_r(:,i)=[-1;1];
        end
    end
    tp=0;tn=0;fp=0;fn=0;
    for i=1:length(tr_tr)
        if (tr_r(1,i)==1) + (tr_tr(1,i)==1) ==2 
            tp=tp+1;
        elseif (tr_r(1,i)==-1) + (tr_tr(1,i)==1) ==2 
            fn=fn+1;
        elseif (tr_r(1,i)==-1) + (tr_tr(1,i)==-1) ==2
            tn=tn+1;
        elseif (tr_r(1,i)==1) + (tr_tr(1,i)==-1) ==2 
            fp=fp+1;
        end
    end
tr_mcc=((tp*tn)-(fp*fn))/((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))^.5
tr_sensitivity=tp/(tp+fn)*100; %tp/total positive
tr_precision=tp/(tp+fp); 
tr_recall=tp/(tp+fn);
tr_specificity=tn/(tn+fp)*100; %tn/total negative
tr_accu=(tp+tn)/length(tr_in)*100;
tr_r=[];

result=sim(network1,ts_in);
    for i=1:length(result)
        if result(1,i) > result(2,i)
            tr_r(:,i)=[1;-1];
        else
            tr_r(:,i)=[-1;1];
        end
    end
    tp=0;tn=0;fp=0;fn=0;
    for i=1:length(ts_tr)
        if (tr_r(1,i)==1) + (ts_tr(1,i)==1) ==2 
            tp=tp+1;
        elseif (tr_r(1,i)==-1) + (ts_tr(1,i)==1) ==2 
            fn=fn+1;
        elseif (tr_r(1,i)==-1) + (ts_tr(1,i)==-1) ==2
            tn=tn+1;
        elseif (tr_r(1,i)==1) + (ts_tr(1,i)==-1) ==2 
            fp=fp+1;
        end
    end
ts_mcc=((tp*tn)-(fp*fn))/((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))^.5
ts_sensitivity=tp/(tp+fn)*100; %tp/total positive
ts_precision=tp/(tp+fp); 
ts_recall=tp/(tp+fn);%=sensitivity
ts_specificity=tn/(tn+fp)*100; %tn/total negative
ts_accu=(tp+tn)/length(ts_in)*100;
tr_r=[];