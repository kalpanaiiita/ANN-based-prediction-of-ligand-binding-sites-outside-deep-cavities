pdb2=pdbread('E:\28thjuly23\CRI\pdb_fd\all\26_pro\1mvp.pdb');
itw='E:\28thjuly23\CRI\pdb_fd\all\26_pro\d_1mvp.txt';
itw2='E:\28thjuly23\CRI\pdb_fd\all\26_pro\1mvp_c.txt';
ddi=dir('E:\28thjuly23\CRI\pdb_fd\all\26_pro\*.txt');
for i =1:length([pdb2.Model.Atom(:).AtomSerNo])
    pdb2.Model.Atom(i).AtomSerNo=i;
end           
                     
%             cg=[0 0 0];
            cg=[mean([pdb2.Model.Atom(:).X]'),mean([pdb2.Model.Atom(:).Y]'),mean([pdb2.Model.Atom(:).Z]')];
% 
            final_res_num=[pdb2.Model.Atom(:).AtomSerNo];
            c_n=19; 

            [cri,in2,cl_i,dd1]=kmeancl99(pdb2,final_res_num,c_n,cg);
            
            r=unique([cl_i{in2(1:2)}]);
                   
            dd=(([pdb2.Model.Atom(r(:)).X]-mean([pdb2.Model.Atom(r(:)).X])).^2 + ([pdb2.Model.Atom(r(:)).Y]-mean([pdb2.Model.Atom(r(:)).Y])).^2 + ([pdb2.Model.Atom(r(:)).Z]-mean([pdb2.Model.Atom(r(:)).Z])).^2).^0.5;
            mv_half=(max(dd)*(2/3));dd=[];           
            mv_all=find(dd1<=(mv_half));                 
                                  
            mv_all2=unique([mv_all,r]);            
             pdb=pdb2;
            [coordi1,ori_ind] = pdb_inv2(pdb);
            if pdb2.Model.Atom(end).AtomSerNo <= 1000
                pr=.5; cr=4;
                final_res_num2 = cyl_surf(coordi1,ori_ind,pr,cr);
            elseif pdb2.Model.Atom(end).AtomSerNo > 1000 && pdb2.Model.Atom(end).AtomSerNo <= 1500
                pr=.6; cr=4;
                final_res_num2 = cyl_surf(coordi1,ori_ind,pr,cr);
            elseif pdb2.Model.Atom(end).AtomSerNo > 1500 && pdb2.Model.Atom(end).AtomSerNo <= 3000  
                 pr=.8; cr=4;
                final_res_num2 = cyl_surf(coordi1,ori_ind,pr,cr);
            elseif pdb2.Model.Atom(end).AtomSerNo > 3000 %&& pdb2.Model.Atom(end).AtomSerNo <=5000
                 pr=.4; cr=5;%pr.3;cr=5
                final_res_num2 = cyl_surf(coordi1,ori_ind,pr,cr);    
%             elseif pdb2.Model.Atom(end).AtomSerNo > 5000
%                 pr=.4; cr=5;%pr.4;cr=5
%                 final_res_num2 = cyl_surf(coordi1,ori_ind,pr,cr);
            end
            in_mol=final_res_num2;
            
            clu_sur=[];
            for i=1:length(mv_all2)
                if length(find(in_mol==mv_all2(i)))>0
                    clu_sur=[clu_sur,mv_all2(i)];
                end
            end
         
%             [cl_br1 ws]=match_res2_whole2(pdb2,clu_sur);
            [mv_nb mv_br1 ws]=match_res2_whole222(pdb2,mv_all2);
            [mo_nb mo_br1 ws]=match_res2_whole222(pdb2,in_mol); 
            
            clu_sur_b=[];
            for i=1:length(mv_br1)            
                if sum(strcmp(mo_br1,mv_br1(i)))>0
                    clu_sur_b=[clu_sur_b,mv_br1(i)];
                end
            end
            AAindex ={'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR'};                   
            AAin={'A' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'K' 'L' 'M' 'N' 'P' 'Q' 'R' 'S' 'T' 'V' 'W' 'Y'};
            clu_sur_nb=[];
            for i=1:length(mv_nb)            
                if sum(strcmp(mo_nb,mv_nb(i)))>0  
                    clu_sur_nb=[clu_sur_nb,mv_nb(i)];
                end
            end     
            
            dssp = dssp_ff22(itw,clu_sur_b);
            con = con_surf2(itw2,clu_sur_b);
            [bi2_pdb all_amino_b]=feature_extraction2(pdb2,clu_sur_b,clu_sur_nb,con,dssp,final_res_num,cg);
            dssp = dssp_ff22(itw,clu_sur_nb);
            con = con_surf2(itw2,clu_sur_nb);
            [bi_n2_pdb all_amino_nb]=feature_extraction2(pdb2,clu_sur_nb,clu_sur_nb,con,dssp,final_res_num,cg);
            
       
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sll sbb]=size(bi2_pdb);
all_pdb=[bi2_pdb;bi_n2_pdb];
for i=1:sbb
    all_pdb(:,i)=(all_pdb(:,i)-all_max_min(i,2))/(all_max_min(i,1)-all_max_min(i,2))*(1-0)+0;
end
[bb bl]=size(bi2_pdb);
[nb nl]=size(bi_n2_pdb);
t_ac=[ones(bb,1),ones(bb,1)*-1];
t_nac=[ones(nb,1)*-1,ones(nb,1)];
all_pdb_tr=[t_ac;t_nac];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts_in=all_pdb';
ts_tr=all_pdb_tr';
result=sim(network1,ts_in);
    for i=1:length(result)
        if result(1,i) > result(2,i)
            tr_r(:,i)=[1;-1];
        else
            tr_r(:,i)=[-1;1];
        end
    end
    tp=0;tn=0;fp=0;fn=0;ts_bb=[];tsb=[];
    for i=1:length(ts_tr)
        if (tr_r(1,i)==1) + (ts_tr(1,i)==1) ==2 
            tp=tp+1;
            ts_bb=[ts_bb,i];
        elseif (tr_r(1,i)==-1) + (ts_tr(1,i)==1) ==2 
            fn=fn+1;
        elseif (tr_r(1,i)==-1) + (ts_tr(1,i)==-1) ==2
            tn=tn+1;
        elseif (tr_r(1,i)==1) + (ts_tr(1,i)==-1) ==2 
            fp=fp+1;
            tsb=[tsb,i];
        end
    end
ts_mcc=((tp*tn)-(fp*fn))/((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))^.5
ts_sensitivity=tp/(tp+fn)*100; %tp/total positive
ts_precision=tp/(tp+fp); 
ts_recall=tp/(tp+fn);%=sensitivity
ts_specificity=tn/(tn+fp)*100; %tn/total negative
ts_accu=(tp+tn)/length(ts_in)*100;
all_amo=[all_amino_b,all_amino_nb];
binding_site_res=all_amino_b
predicted_binding_site_res=all_amo(ts_bb)
all_amo(tsb)