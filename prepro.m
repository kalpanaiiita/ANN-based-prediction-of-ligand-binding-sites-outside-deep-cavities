% load('both.mat');
[sll sbb]=size(bi_all2);
all=[bi_all2;bi_n_all2];
for i=1:sbb
    all_max_min(i,:)=[max(all(:,i)),min(all(:,i))];
    all(:,i)=(all(:,i)-min(all(:,i)))/(max(all(:,i)-min(all(:,i))))*(1-0)+0;   
end
% for i=1:19
%     all(:,i)=(all(:,i)-mean(all(:,i)))/(sum(all(:,i)-mean(all(:,i)))/length(all));
% end

n_ac=all(1:length(bi_all2),:);
n_nac=all(length(bi_all2)+1:end,:);

t_ac=[ones(length(bi_all2),1),ones(length(bi_all2),1)*-1];
t_nac=[ones(length(bi_n_all2),1)*-1,ones(length(bi_n_all2),1)];
[v br]=sort(rand(1,length(bi_all2)));
[v nbr]=sort(rand(1,length(bi_n_all2)));
n_ac=n_ac(br,:);
n_nac=n_nac(nbr,:);

% ts_in=[n_ac(1:floor(30*length(n_ac)/100),:);n_nac(1:floor(30*length(n_nac)/100),:)]';
% tr_in=[n_ac(floor(30*length(n_ac)/100)+1:end,:);n_nac(floor(30*length(n_nac)/100)+1:end,:)]';
% ts_tr=[t_ac(1:floor(30*length(n_ac)/100),:);t_nac(1:floor(30*length(n_nac)/100),:)]';
% tr_tr=[t_ac(floor(30*length(n_ac)/100)+1:end,:);t_nac(floor(30*length(n_nac)/100)+1:end,:)]';

ts_in=[n_ac(1:floor(30*length(n_ac)/100),:);n_nac(1:floor(30*length(n_ac)/100),:)]';
tr_in=[n_ac(floor(30*length(n_ac)/100)+1:end,:);n_nac(floor(30*length(n_ac)/100)+1:length(n_ac),:)]';
ts_tr=[t_ac(1:floor(30*length(n_ac)/100),:);t_nac(1:floor(30*length(n_ac)/100),:)]';
tr_tr=[t_ac(floor(30*length(n_ac)/100)+1:end,:);t_nac(floor(30*length(n_ac)/100)+1:length(n_ac),:)]';

% ts_in=[n_ac(1:floor(30*length(n_nac)/100),:);n_nac(1:floor(30*length(n_nac)/100),:)]';
% tr_in=[n_ac(floor(30*length(n_nac)/100)+1:end,:);n_nac(floor(30*length(n_nac)/100)+1:length(n_nac),:)]';
% ts_tr=[t_ac(1:floor(30*length(n_nac)/100),:);t_nac(1:floor(30*length(n_nac)/100),:)]';
% tr_tr=[t_ac(floor(30*length(n_nac)/100)+1:end,:);t_nac(floor(30*length(n_nac)/100)+1:length(n_nac),:)]';


% wr=rand(length(tr_in),1);
% [v r_tr]=sort(wr);
% tr_in=tr_in(r_tr,:)';
% tr_tr=tr_tr(r_tr,:)';
% 
% wr=rand(length(ts_in),1);
% [v r_ts]=sort(wr);
% ts_in=ts_in(r_ts,:)';
% ts_tr=ts_tr(r_ts,:)';
