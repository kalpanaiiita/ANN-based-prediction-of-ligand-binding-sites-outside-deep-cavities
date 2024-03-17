function final_res_num2 = cyl_surf(coordi1,ori_ind,pr,cr)
sur_res_num = [];
[stlPoints Apoints Bpoints Cpoints ]= sphereTriangulation(cr, 1);
% figure;
% c=ones(length(stlPoints),1);
% scatter3(stlPoints(:,1),stlPoints(:,2),stlPoints(:,3),10,c,'filled');       
xyz=stlPoints; 
% load('xyz.mat');
for j = 1:length(xyz)       % as angle with X axis is beta
        l = (xyz(j,1));     %     l = (xyz(j,1)) / sqrt((xyz(j,1))^2 + (xyz(j,2))^2 + (xyz(j,3))^2);
        m = (xyz(j,2));     %     m = (xyz(j,2)) / sqrt((xyz(j,1))^2 + (xyz(j,2))^2 + (xyz(j,3))^2);
        n = (xyz(j,3));     %     n= (xyz(j,3)) / sqrt((xyz(j,1))^2 + (xyz(j,2))^2 + (xyz(j,3))^2);
    for ii = 1:length(ori_ind)    
        r = ((l)*(coordi1(ii,1)) + (m)*(coordi1(ii,2)) + (n)*(coordi1(ii,3)))/ (l*l + m*m +n*n);
        x_line = l*r;
        y_line = m*r;
        z_line = n*r;
        dist_rr(ii) = sqrt((((coordi1(ii,1))-(x_line))^2) + (((coordi1(ii,2))-(y_line))^2) + (((coordi1(ii,3))-(z_line))^2));
    end
    index=find(dist_rr<=pr);           % index for position of res. in coordi1
    get_cord=coordi1(index,:);
    get_cord_i = ori_ind(index);    
    dist_c_r = [];
    for ii=1:length(index)
         dist_c_r(ii) = sqrt(((get_cord(ii,1))^2) + ((get_cord(ii,2))^2) + ((get_cord(ii,3))^2));
    end    
    [v ind11] = max(dist_c_r);       
    dist_c_r = [];
    for ii=1:length(index)
        dist_c_r(ii) = sqrt((get_cord(ind11,1)-get_cord(ii,1))^2+(get_cord(ind11,2)-get_cord(ii,2))^2+(get_cord(ind11,3)-get_cord(ii,3))^2);
    end
    [dist_max ind112] = max(dist_c_r);       
    sur_res_num=[sur_res_num get_cord_i(ind11)];
    sur_res_num=[sur_res_num get_cord_i(ind112)]; 
%     sur_res_cord = [coordi1(ind111,:);coordi1(ind112,:)];
end
final_res_num2=unique(sur_res_num);
sur_coordi1=coordi1(final_res_num2',:);

    
    
    
    
  