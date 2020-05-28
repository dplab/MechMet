function sx_vec_rstiff =  RotateStiffnessMatrix(sx_vec_stiff,myrot)

%  Rotates the elastic stiffness matrix in crystal bases to bases 
%  given by myrot 
%  Operates on 6x6 matrix form using method in Ting's Anisotropic
%  Elasticity Text
%


k_1 = zeros(3,3);
k_2 = zeros(3,3);
k_3 = zeros(3,3);
k_4 = zeros(3,3);

k_full = zeros(6,6);


k_1  =  [   myrot(1,1)^2       myrot(1,2)^2      myrot(1,3)^2  ; ...
            myrot(2,1)^2       myrot(2,2)^2      myrot(2,3)^2  ; ...
            myrot(3,1)^2       myrot(3,2)^2      myrot(3,3)^2  ] ;
        
        
k_2  =  [   myrot(1,2)*myrot(1,3)       myrot(1,3)*myrot(1,1)      myrot(1,1)*myrot(1,2)  ; ...
            myrot(2,2)*myrot(2,3)       myrot(2,3)*myrot(2,1)      myrot(2,1)*myrot(2,2)  ; ...
            myrot(3,2)*myrot(3,3)       myrot(3,3)*myrot(3,1)      myrot(3,1)*myrot(3,2)  ] ;
                
        
k_3  =  [   myrot(2,1)*myrot(3,1)       myrot(2,2)*myrot(3,2)      myrot(2,3)*myrot(3,3)  ; ...
            myrot(3,1)*myrot(1,1)       myrot(3,2)*myrot(1,2)      myrot(3,3)*myrot(1,3)  ; ...
            myrot(1,1)*myrot(2,1)       myrot(1,2)*myrot(2,2)      myrot(1,3)*myrot(2,3)  ] ;
        
k_4  =  [   myrot(2,2)*myrot(3,3)+myrot(2,3)*myrot(3,2)       myrot(2,3)*myrot(3,1)+myrot(2,1)*myrot(3,3)      myrot(2,1)*myrot(3,2)+myrot(2,2)*myrot(3,1)  ; ...
            myrot(3,2)*myrot(1,3)+myrot(3,3)*myrot(1,2)       myrot(3,3)*myrot(1,1)+myrot(3,1)*myrot(1,3)      myrot(3,1)*myrot(1,2)+myrot(3,2)*myrot(1,1)  ; ...
            myrot(1,2)*myrot(2,3)+myrot(1,3)*myrot(2,2)       myrot(1,3)*myrot(2,1)+myrot(1,1)*myrot(2,3)      myrot(1,1)*myrot(2,2)+myrot(1,2)*myrot(2,1)  ] ;
        
        
k_full = [ k_1 2*k_2; k_3 k_4];

sx_vec_rstiff = k_full*sx_vec_stiff*k_full';
