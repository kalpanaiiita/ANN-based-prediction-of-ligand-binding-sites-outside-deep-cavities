function [cri,in2,cl_i,dd1]=kmeancl99(pdb2,final_res_num,c_n,cg)
            dd1=(([pdb2.Model.Atom(:).X]-cg(1)).^2 + ([pdb2.Model.Atom(:).Y]-cg(2)).^2 + ([pdb2.Model.Atom(:).Z]-cg(3)).^2).^0.5;
            [mn md]=sort(dd1);
            c{1,1}=[pdb2.Model.Atom(md(1)).X,pdb2.Model.Atom(md(1)).Y,pdb2.Model.Atom(md(1)).Z];  c{2,1}=[pdb2.Model.Atom(md(end)).X,pdb2.Model.Atom(md(end)).Y,pdb2.Model.Atom(md(end)).Z];
            md_l=round(length(md)/(c_n-1));
            for i=1:c_n-2
                c{i+2,1}=[pdb2.Model.Atom(md(md_l*i)).X,pdb2.Model.Atom(md(md_l*i)).Y,pdb2.Model.Atom(md(md_l*i)).Z];  
            end
            for i=1:c_n
                cl{i,1}=[];      cl_i{i,1}=[];        cl_b{i,1}=[];  cl_br{i,1}=[];
            end
            while 1
                for i=1:length(final_res_num)
                    for j=1:c_n
                        d(j)=((pdb2.Model.Atom(i).X-c{j}(1))^2 + (pdb2.Model.Atom(i).Y-c{j}(2))^2 + (pdb2.Model.Atom(i).Z-c{j}(3))^2)^.5;
                    end
                    [m_d m_v]=min(d);
                    cl{m_v}=[cl{m_v};[pdb2.Model.Atom(i).X,pdb2.Model.Atom(i).Y,pdb2.Model.Atom(i).Z]]; cl_i{m_v}=[cl_i{m_v};i];
                end
                for j=1:c_n
                    new_c{j,1}=[mean(cl{j}(:,1)),mean(cl{j}(:,2)), mean(cl{j}(:,3))];
                end               
                if isequal(c,new_c)==1
                    break
                else
                    c=new_c;
                    for j=1:c_n
                        cl{j}=[];
                        cl_i{j}=[];
                    end
                end                
            end
            for j=1:c_n
                cl_i{j}=final_res_num(cl_i{j});
            end
            
            for j=1:c_n
                [lx,bx]=size(cl{j,1});
                cl_l(j)=lx;
                d_a=((cl{j}(:,1)-cg(1)).^2 + (cl{j}(:,2)-cg(2)).^2 + (cl{j}(:,3)-cg(3)).^2).^0.5;
                dis_a{j}=d_a;
                d_a=[];
            end
            for j=1:c_n
                sd{j}=dis_a{j}-mean(dis_a{j});
            end
            for j=1:c_n
                sdd(j)=sqrt(sum(sd{j}.^2)/length(sd{j}));
            end
            
            [v2 in2]=sort(sdd); 
            in2=fliplr(in2);
            cri=mean(sdd);

            
            
