function [all_val all_amino]=feature_extraction2(pdb2,all_cs,clu_sur_nb,con,dssp,final_res_num,cg)

            
            AAindex ={'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR'};                   
            atoma={'OD1','OD2','OE2','NE2','ND1','OE1','OG','OG1','OH'};
            atomd={'NE','NE1','NE2','ND1','ND2','NH1','NH2','OG','OG1','OH','NZ'};
            atomel={'C' 'N' 'O'};           
            load('aa_features.mat');
            aa=aa_features;
            all_val=[];all_amino=[];
            for i=1:length(all_cs)
                ax=all_cs{i};             
                cs=find([pdb2.Model.Atom(:).resSeq]==str2num(ax(1:end-5)));
                cs=cs(find(strcmp({pdb2.Model.Atom(cs).chainID},ax(end-1))==1));  
                if strcmp(ax(end),' ')
                else
                 cs=cs(find(strcmp({pdb2.Model.Atom(cs).iCode},ax(end))==1));
                end
                
                all_co=[mean([pdb2.Model.Atom(cs).X]),mean([pdb2.Model.Atom(cs).Y]),mean([pdb2.Model.Atom(cs).Z])];
                if cs(end)==final_res_num(end)
                    alla=[final_res_num(1):final_res_num(cs(1)-1)];
                elseif cs(1)==final_res_num(1)
                    alla=[final_res_num(cs(end)+1):final_res_num(end)];
                else                
                    alla=[final_res_num(1):final_res_num(cs(1)-1),final_res_num(cs(end)+1):final_res_num(end)];
                end                                                
                cen=sum((([pdb2.Model.Atom(alla).X]-all_co(1)).^2 + ([pdb2.Model.Atom(alla).Y]-all_co(2)).^2 + ([pdb2.Model.Atom(alla).Z]-all_co(3)).^2).^.5);
                cenr=(length([pdb2.Model.Atom(:).resSeq])-1)/cen;
                ds2=mean((([pdb2.Model.Atom(cs).X]-cg(1)).^2 + ([pdb2.Model.Atom(cs).Y]-cg(2)).^2 + ([pdb2.Model.Atom(cs).Z]-cg(3)).^2).^.5);
                bf=mean([pdb2.Model.Atom(cs).tempFactor]);
                
                
                dvv=[];
                for ui=1:length(cs)
                    dv=find((([pdb2.Model.Atom(:).X]-pdb2.Model.Atom(cs(ui)).X).^2 + ([pdb2.Model.Atom(:).Y]-pdb2.Model.Atom(cs(ui)).Y).^2 + ([pdb2.Model.Atom(:).Z]-pdb2.Model.Atom(cs(ui)).Z).^2).^.5<=1)';
                    dv2=(([pdb2.Model.Atom(:).X]-pdb2.Model.Atom(cs(ui)).X).^2 + ([pdb2.Model.Atom(:).Y]-pdb2.Model.Atom(cs(ui)).Y).^2 + ([pdb2.Model.Atom(:).Z]-pdb2.Model.Atom(cs(ui)).Z).^2).^.5';
                    dvv=[dvv;[dv,dv2(dv)]];                        
                end 
                [dvvv dis]=unique(dvv(:,1));
                [dss dsi]=sort(dvv(dis,2));
                dvv2=dvvv(dsi);
                
                dvv1=[];
                for gi=1:length(dvv2)
                    if sum(strcmp({'VAL','LYS','PRO','TRP','ARG','GLY','GLU','ILE','ASN','ALA','ASP','THR','SER','MET','GLN','PHE','HIS','LEU','CYS','TYR'},pdb2.Model.Atom(dvv2(gi)).resName))>0
                        if length(pdb2.Model.Atom(dvv2(gi)).iCode)==0
                            dvv1=[dvv1,{[num2str(pdb2.Model.Atom(dvv2(gi)).resSeq),pdb2.Model.Atom(dvv2(gi)).resName,pdb2.Model.Atom(dvv2(gi)).chainID,' ']}];
                        else
                            dvv1=[dvv1,{strcat(num2str(pdb2.Model.Atom(dvv2(gi)).resSeq),pdb2.Model.Atom(dvv2(gi)).resName,pdb2.Model.Atom(dvv2(gi)).chainID,pdb2.Model.Atom(dvv2(gi)).iCode)}];
                        end
                    end
                end   
                dvvv1=dvv1(1);
                for gi=2:length(dvv1)
                    if sum(strcmp(dvvv1,dvv1(gi)))<1
                        dvvv1(end+1)=dvv1(gi);
                    end
                end
                dvvx1=[];
                if length(dvvv1)>3
                    for hc=2:length(dvvv1)
                        for hcx=(hc+1):length(dvvv1)
                            if sum(strcmp(clu_sur_nb,dvvv1(1)))==1 
                                if sum(strcmp(clu_sur_nb,dvvv1(hc)))==1 && sum(strcmp(clu_sur_nb,dvvv1(hcx)))==1
                                    dvvx1=[dvvx1;[dvvv1(1),dvvv1(hc),dvvv1(hcx)]];
                                end
                            else
                                dvvx1=[dvvx1;[dvv1(1),dvvv1(hc),dvvv1(hcx)]];
                            end
                        end
                    end
                else
                    dvvx1=dvvv1;
                end
                
                [zl,zb]=size(dvvx1);   
                all_val2=[];all_amino2=[];
                for gxi=1:zl
                    dvl=dvvx1(gxi,:);
                    dvx=[];
                    for ggi=1:zb
                        ax=dvl{ggi};
                        ccs=find([pdb2.Model.Atom(:).resSeq]==str2num(ax(1:end-5)));
                        ccs=ccs(find(strcmp({pdb2.Model.Atom(ccs).chainID},ax(end-1))==1));
                        if strcmp(ax(end),' ')
                        else
                            ccs=ccs(find(strcmp({pdb2.Model.Atom(ccs).iCode},ax(end))==1));
                        end
                        dvx=[dvx,ccs];
                        cs2(ggi)=find(strcmp(AAindex,ax(end-4:end-2))==1);
                     end
                               
                    drC=sum(strcmp({pdb2.Model.Atom(dvx).element},atomel{1}));
                    drN=sum(strcmp({pdb2.Model.Atom(dvx).element},atomel{2}));
                    drO=sum(strcmp({pdb2.Model.Atom(dvx).element},atomel{3}));
                
                    rsn={pdb2.Model.Atom(dvx).AtomName};                    
                    drd=0;
                    for ir=1:length(atomd)
                        drd=drd+sum(strcmp(rsn,atomd{ir}));
                    end
                    dra=0;
                    for ir=1:length(atoma)
                        dra=dra+sum(strcmp(rsn,atoma{ir}));
                    end
                if length(cs2)==1
                    csx=aa(cs2,:);
                else
                    csx=sum(aa(cs2,:));
                end
                all_val2(gxi,:)= [csx,[ds2,bf,cenr],dssp(i,:),drC,drN,drO,drd,dra,con(i)];   
                all_amino2=[all_amino2,dvl(1)];
                drd=0;dra=0;drC=0;drN=0;drO=0;cs2=[];
                end
                all_val=[all_val;all_val2];
                all_amino=[all_amino,all_amino2];
                all_val2=[];dvvx1=[];dvv=[];dvv1=[];  dvvv1=[];dvx=[];cs=[];all_amino2=[];
            end