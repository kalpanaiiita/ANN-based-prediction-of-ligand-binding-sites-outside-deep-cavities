function[mv_nb cl_br1 ws]=match_res2_whole222(pdb2,mv_all2)

       Ain={'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'};
                 
            aw=[pdb2.Site.SiteDetail];            
            for lll=1:length(aw)
                a1{lll}=[aw(lll).ResDet(:).ResSeqNo];
                a2{lll}={aw(lll).ResDet(:).ChainID};
                a3{lll}={aw(lll).ResDet(:).ResName};
                a4{lll}={aw(lll).ResDet(:).InsCode};
            end 
            sr=[];a444=[];
            for ll=1:length(aw)    
                a11=a1{ll};
                a22=a2{ll};
                a33=a3{ll};
                a44=a4{ll};
                for i=1:length(a11)
                    if sum(strcmp({'VAL','LYS','PRO','TRP','ARG','GLY','GLU','ILE','ASN','ALA','ASP','THR','SER','MET','GLN','PHE','HIS','LEU','CYS','TYR'},a33(i)))>0
                        sr=[sr,strcat(num2str(a11(i)),(a33(i)),a22(i),a44(i))];
                        a444=[a444, a44{i}];
                    end
                end
            end
%             cee=[];
%             for ce=1:length(a444)
%                 if sum(strcmp(Ain,a444(ce)))==1
%                     cee=[cee,ce];
%                 end
%             end           
%                 
            sr=unique(sr);
            ws=length(sr);
            cl_br11=[];
            for i=1:length(mv_all2)
                if sum(strcmp({'VAL','LYS','PRO','TRP','ARG','GLY','GLU','ILE','ASN','ALA','ASP','THR','SER','MET','GLN','PHE','HIS','LEU','CYS','TYR'},pdb2.Model.Atom(mv_all2(i)).resName))>0
%                     if ~isempty(strcmp(sr,strcat(num2str(pdb2.Model.Atom(mv_all2(i)).resSeq),pdb2.Model.Atom(mv_all2(i)).chainID)))
                    if length(pdb2.Model.Atom(mv_all2(i)).iCode)==0
                        cl_br11=[cl_br11,{[num2str(pdb2.Model.Atom(mv_all2(i)).resSeq),pdb2.Model.Atom(mv_all2(i)).resName,pdb2.Model.Atom(mv_all2(i)).chainID,' ']}];
%                     end
                    else
                        cl_br11=[cl_br11,{strcat(num2str(pdb2.Model.Atom(mv_all2(i)).resSeq),pdb2.Model.Atom(mv_all2(i)).resName,pdb2.Model.Atom(mv_all2(i)).chainID,pdb2.Model.Atom(mv_all2(i)).iCode)}];
                    end
                end
            end
            cl_br11=unique(cl_br11);
            cl_br1=[];mv_nb=[];
            for i=1:length(cl_br11)
                if sum(strcmp(sr,cl_br11(i)))>0
                    cl_br1=[cl_br1,cl_br11(i)];
                else
                    mv_nb=[mv_nb,cl_br11(i)];
                end
            end
            
            
%             for ll=1:length(aw)    
%                 a11=a1{ll};
%                 a22=a2{ll};
%                 for li=1:length(pdb2.Model.Atom)                    
%                        for i=1:length(a11)
%                            if a11(i)==pdb2.Model.Atom(li).resSeq & strcmp(a22{i},pdb2.Model.Atom(li).chainID)==1
%                                wws=[wws;pdb2.Model.Atom(li).AtomSerNo];      
%                                wws2=[wws2;pdb2.Model.Atom(li).resSeq]; 
%                            end
%                        end
%                 end
%             end 
%             [x xx]=unique(wws);
%             wws=wws(xx);
%             wws2=wws2(xx);
%             ws=wws2(1); 
%             for i=1:length(wws2)
%                 if ws(end)==wws2(i)
%                 else
%                     ws=[ws,wws2(i)];
%                 end
%             end
%             ws=length(ws);       
%             
%             cl_br=[];
%             for j=1:length(mv_all2)   
%                 if length(find(wws==mv_all2(j)))>0 
%                     cl_br=[cl_br,mv_all2(j)];
%                 end
%             end
%              
%             cl_br=[pdb2.Model.Atom(cl_br).resSeq];
% %             cl_br0=[pdb2.Model.Atom(cl_br).AtomSerNo];
%             if length(cl_br)>0
%             cl_br2=cl_br(1);%cl_br3=cl_br0(1);
%             for i=1:length(cl_br)-1
%                 if cl_br2(end)==cl_br(i)
%                 else
%                     cl_br2=[cl_br2,cl_br(i)];
% %                     cl_br3=[cl_br3,cl_br0(i)];
%                 end
%             end
%             cl_br2=length(cl_br2);  
%             else
%                 cl_br2=0;
%             end
% %             cl_br1=pdb2.Model.Atom(cl_br3);
% 
% 
%             
%             
%             
%             