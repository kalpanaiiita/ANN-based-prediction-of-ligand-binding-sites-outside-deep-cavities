function [stlPoints Apoints Bpoints Cpoints ]= sphereTriangulation(numIterations, radius)
A=[1 0 0]*radius;
B=[0 1 0]*radius;
C=[0 0 1]*radius;
triangles=[ A;  B;  C;...
            A;  B; -C;... 
           -A;  B;  C;...
           -A;  B; -C;...
           -A; -B;  C;...
           -A; -B; -C;...
            A; -B;  C;...
            A; -B; -C];
 selector         = 1:3:numel(triangles(:,1))-1 ;
 Apoints = triangles(selector  ,:)    ;
 Bpoints = triangles(selector+1,:)    ;
 Cpoints = triangles(selector+2,:)    ;
for iteration = 1:numIterations    
    AB_2   =  (Apoints+Bpoints)/2;
    AB_2   =  arsUnit(AB_2,radius);...
    AC_2   =  (Apoints+Cpoints)/2; AC_2   =  arsUnit(AC_2,radius);
    CB_2   =  (Cpoints+Bpoints)/2; CB_2   =  arsUnit(CB_2,radius);
    Apoints = [ Apoints;... 
                AB_2;...    
                AC_2;...    
                AC_2];...   
    Bpoints = [AB_2; Bpoints; AB_2; CB_2   ];
    Cpoints = [AC_2; CB_2   ; CB_2; Cpoints];
end

numPoints = numel(Apoints(:,1))                                 ; 
selector  = 1:numPoints                                         ;
selector  = selector(:)                                         ;    
selector  = [selector, selector+numPoints, selector+2*numPoints];    
selector  = selector'                                           ;
selector  = selector(:)                                         ;
stlPoints = [Apoints; Bpoints; Cpoints]                         ;
stlPoints = stlPoints(selector,:)                               ;
end

function [rez]=arsNorm(A)
    rez = A(:,1).*A(:,1) + A(:,2).*A(:,2) + A(:,3).*A(:,3);
	rez = sqrt(rez);
end 

function [rez]=arsUnit(A,radius)
	normOfA = arsNorm(A)                      ;
	rez       = A ./ [normOfA normOfA normOfA];
    rez       = rez * radius                  ;
end











