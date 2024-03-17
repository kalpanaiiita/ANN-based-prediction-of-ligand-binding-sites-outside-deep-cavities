function [coordi1,ori_ind] = pdb_inv2(pdb)
%%%% programme for converting the coordinate of pdb to invarient%%%%
%%%% co-ordinate%%%% system
ori_ind=[];
j=1;
for i=1:length(pdb.Model.Atom)
    xx(i,:)=[pdb.Model.Atom(1,i).X,pdb.Model.Atom(1,i).Y,pdb.Model.Atom(1,i).Z];%if (strcmp(pdb.Model.Atom(1,i).AtomName,'CA')==1)
    x(j)=pdb.Model.Atom(1,i).X;
    y(j)=pdb.Model.Atom(1,i).Y;
    z(j)=pdb.Model.Atom(1,i).Z;
    ori_ind=[ori_ind;pdb.Model.Atom(1,i).AtomSerNo];
    j = j+1;%     end
end
coord = [x;y;z];
coordi = coord';    %coord'is used for row to colomn interconvertion;
cx = mean(x); cy = mean(y); cz = mean(z);
cg = [cx cy cz]; %% assigning CG
x = x-cx; y = y - cy; z = z - cz;
d = (x.^2 + y.^2 + z.^2).^0.5;  %%finding out the distance from CG to each points
ind0 = find(d == max(d));
x_z = coordi((ind0),1); y_z = coordi((ind0),2); z_z = coordi((ind0),3);
Z=[x_z y_z z_z]; %%max point and Z axis
%%%finding equation of plane passing through cg and perpindicular to the
%%%normal from Z = maximum distance findind direction ratio
a=x_z-cx;
b=y_z-cy;
c=z_z-cz;
%%so equation of plane ax+by+cz=k having cg as a point in it page 51
k=((a*cx)+(b*cy)+(c*cz));
% %%perpendicular distance of points from given plane passing through cg
%%% direction cosine of the normal to the plane
l=a/(sqrt(a^2+b^2+c^2));
m=b/(sqrt(a^2+b^2+c^2));
n=c/(sqrt(a^2+b^2+c^2));
p=k/(sqrt(a^2+b^2+c^2));
%%finding distance from all points to the above plane
lenn=length(coordi);
for i=1:lenn    
    rr(i)=(l*coordi(i,1)) + (m*coordi(i,2)) + (n*coordi(i,3))-p;
    r(i)=abs(rr(i));
end
%%%finding maximum distance point among all points of distance 2
j=1;
for i=1:lenn
    if (r(i)<= 2)
        r1(j)=r(i);
        index(j)=i;
        j=j+1;
   end
end    
for j=1:length(index)    
    distance(j)=sqrt( (coordi(index(j),1)-cx)^2+(coordi(index(j),2)-cy)^2+(coordi(index(j),3)-cz)^2);
end
ind1 = find(distance == max(distance));
x_x = coordi(index(ind1),1);
y_x = coordi(index(ind1),2);
z_x = coordi(index(ind1),3);
%%%X_P=[coordi(index(ind1),1) coordi(index(ind1),2) coordi(index(ind1),3)];
X=[x_x y_x z_x];
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % getting the perpendicular from a given point to a given line
% % from page 115  and 90
% 
% % dist = sqrt((x_z-cx)^2+(y_z-cy)^2+(z_z-cz)^2);
% 
% % l = (x_z-cx)/dist3;
% % m = (y_z-cy)/dist3;
% % n = (z_z-cz)/dist3;
% 
% % equation of line through cg and z
% % ((x-x_z)/(cx-x_z)) = ((y-y_z)/(cy-y_z)) = ((z-z_z)/(cz-z_z));
l =  x_z - cx;
m =  y_z - cy;
n =  z_z - cz;
r = (l*(x_x-cx) + m*(y_x-cy) + n*(z_x-cz))/(l*l + m*m + n*n);
%%% so CG' is cx' cy' cz';
cx = cx + l*r;
cy = cy + m*r;
cz = cz + n*r;
% %%%xx for testing wheather it is cg-X and cg-Z are perpendicular or not
% %%%uncomment it and test
% xx = (cx-X(1))*(cx-Z(1)) + (cy-X(2))*(cy-Z(2)) + (cz-X(3))*(cz-Z(3));
CG = [cx cy cz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%program for getting Y from X.Z and CG
%%%% get help from web page..on dextop
mm = ((y_x-cy)*(z_z-cz)-(z_x-cz)*(y_z-cy));
nn =-((x_x-cx)*(z_z-cz)-(z_x-cz)*(x_z-cx));
pp = ((x_x-cx)*(y_z-cy)-(y_x-cy)*(x_z-cx));
dd =-(mm*cx + nn*cy + pp*cz);
% %%%so eqauation of the plane passes through X,Y and CG would be
%%%%%%%%%%%  mm*x + nn*y + pp*z =dd;
% %%%equation for perpedicular to the given plane
%%% we take CG for co-ordinates
r = 0.05;
x_y = cx - r*mm;
y_y = cy - r*nn;
z_y = cz - r*pp;

Y = [x_y y_y z_y];
%%%% for testing value of  r
%%%%test_r = -(mm*x_y + nn*y_y + pp*z_y + dd)/ (mm^2 + nn^2 + pp^2) ;
%%%%% page 235, 21
dist1 = sqrt((x_x-cx)^2+(y_x-cy)^2+(z_x-cz)^2);
l1 = (x_x-cx)/dist1;
m1 = (y_x-cy)/dist1;
n1 = (z_x-cz)/dist1;
dist2 = sqrt((x_y-cx)^2+(y_y-cy)^2+(z_y-cz)^2);
l2 = (x_y-cx)/dist2;
m2 = (y_y-cy)/dist2;
n2 = (z_y-cz)/dist2;
dist3 = sqrt((x_z-cx)^2+(y_z-cy)^2+(z_z-cz)^2);
l3 = (x_z-cx)/dist3;
m3 = (y_z-cy)/dist3;
n3 = (z_z-cz)/dist3;
% %changed co-ordinates are x', y' and z'
% %original co-ordinates are obtained from pdb in  x, y,z form
for i=1:lenn
coordi1(i,1) = (l1*(coordi(i,1)-cx) + m1*(coordi(i,2)-cy) + n1*(coordi(i,3)-cz));
coordi1(i,2) = (l2*(coordi(i,1)-cx) + m2*(coordi(i,2)-cy) + n2*(coordi(i,3)-cz));
coordi1(i,3) = (l3*(coordi(i,1)-cx) + m3*(coordi(i,2)-cy) + n3*(coordi(i,3)-cz));
end