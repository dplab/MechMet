function shafac_structure = ShapeFunctionArrays()

%generates shape functions and shape function derivatives for 10-node tets

    qploc = zeros(3,15);
    qploc(1,1) =  0.333333333333333333d0; 
    qploc(1,2) =  0.333333333333333333d0; 
    qploc(1,3) =  0.333333333333333333d0; 
    qploc(1,4) =  0.0d0; 
    qploc(1,5) =  0.25d0; 
    qploc(1,6) =  0.909090909090909091d-1; 
    qploc(1,7) =  0.909090909090909091d-1;
    qploc(1,8) =  0.909090909090909091d-1; 
    qploc(1,9) =  0.727272727272727273d0; 
    qploc(1,10) = 0.665501535736642813d-1; 
    qploc(1,11) = 0.665501535736642813d-1; 
    qploc(1,12) = 0.665501535736642813d-1; 
    qploc(1,13) = 0.433449846426335728d0; 
    qploc(1,14) = 0.433449846426335728d0; 
    qploc(1,15) = 0.433449846426335728d0; 

    
    qploc(2,1) =  0.333333333333333333d0;      
    qploc(2,2) =  0.333333333333333333d0; 
    qploc(2,3) =  0.0d0; 
    qploc(2,4) =  0.333333333333333333d0; 
    qploc(2,5) =  0.25d0; 
    qploc(2,6) =  0.909090909090909091d-1; 
    qploc(2,7) =  0.909090909090909091d-1; 
    qploc(2,8) =  0.727272727272727273d0; 
    qploc(2,9) =  0.909090909090909091d-1; 
    qploc(2,10) = 0.665501535736642813d-1; 
    qploc(2,11) = 0.433449846426335728d0; 
    qploc(2,12) = 0.433449846426335728d0; 
    qploc(2,13) = 0.665501535736642813d-1; 
    qploc(2,14) = 0.665501535736642813d-1; 
    qploc(2,15) = 0.433449846426335728d0;

    qploc(3,1) =  0.333333333333333333d0;
    qploc(3,2) =  0.0d0;
    qploc(3,3) =  0.333333333333333333d0;
    qploc(3,4) =  0.333333333333333333d0;
    qploc(3,5) =  0.25d0;
    qploc(3,6) =  0.909090909090909091d-1;
    qploc(3,7) =  0.727272727272727273d0;
    qploc(3,8) =  0.909090909090909091d-1;
    qploc(3,9) =  0.909090909090909091d-1;
    qploc(3,10) = 0.433449846426335728d0; 
    qploc(3,11) = 0.665501535736642813d-1;
    qploc(3,12) = 0.433449846426335728d0;
    qploc(3,13) = 0.665501535736642813d-1;                    
    qploc(3,14) = 0.433449846426335728d0;                   
    qploc(3,15) = 0.665501535736642813d-1; 

    wtqp = zeros(1,15);
    
    wtqp(1:4)   = 0.602678571428571597d-2;
    wtqp(5)     = 0.302836780970891856d-1;
    wtqp(6:9)   = 0.116452490860289742d-1;
    wtqp(10:15) = 0.109491415613864534d-1;
  
      nqpt = 15;
      nnpe = 10;
      nnps =  6;
    
      sfac = zeros(nnpe,nqpt);
      dndxi = zeros(nnpe,nqpt);
      dndet = zeros(nnpe,nqpt);
      dndze = zeros(nnpe,nqpt);

   for i=1:1:nnpe
       for j=1:1:nqpt
           q1= qploc(1,j);
           q2= qploc(2,j);
           q3= qploc(3,j);
    sfac(i,j) =sfn3d(q1,q2,q3,i);
    dndxi(i,j) =fxi3d(q1,q2,q3,i);
    dndet(i,j) =fet3d(q1,q2,q3,i);
    dndze(i,j) =fze3d(q1,q2,q3,i); 
       end
   end
  
 
shafac_structure    = Data_Structure_Shafac(nnpe,nqpt,wtqp,sfac,dndxi,dndet,dndze);
end
function sf=sfn3d(q1,q2,q3,n)
% shape functions

switch n
    case 1
      sf=(1.0-2.0*q1-2.0*q2-2.0*q3)*(1.0-q1-q2-q3);	
    case 2
       sf=4.0*q1*(1.0-q1-q2-q3);
    case 3
      sf=(2.0*q1-1.0)*q1;
    case 4
       sf=4.0*q1*q2;
    case 5
      sf=(2.0*q2-1.0)*q2;
    case 6
       sf=4.0*q2*(1.0-q1-q2-q3);	
    case 7
       sf=4.0*q3*(1.0-q1-q2-q3);
    case 8
       sf=4.0*q1*q3;
    case 9
       sf=4.0*q2*q3;
    case 10
      sf=(2.0*q3-1.0)*q3;
end
end
function fi=fxi3d(q1,q2,q3,n) 
% shape function derivatives with respect to xi 
    
switch n
    case 1
      fi=((-2.0)*(1.0-q1-q2-q3))+((1.0-2.0*q1-2.0*q2-2.0*q3)*(-1.0));    
    case 2
       fi=(4.0*1.0*(1.0-q1-q2-q3))+(4.0*q1*(-1.0));
    case 3
      fi=((2.0)*q1)+((2.0*q1-1.0)*1.0);     
    case 4
       fi=4.0*(1.0)*q2;
    case 5
      fi=0.0;	
    case 6
       fi=4.0*q2*(-1.0);	
    case 7
       fi=4.0*q3*(-1.0);
    case 8
       fi=4.0*(1.0)*q3; 
    case 9
       fi=0.0;
    case 10
      fi=0.0;	 
    end
	  
end 
function ft=fet3d(q1,q2,q3,n) 
% shape function derivative wrt xi
switch n
    case 1
      ft=4.0*(q1+q2+q3) - 3.0;	
    case 2
       ft=4.0*q1*(-1.0);
    case 3
      ft=0.0;
    case 4
       ft=4.0*q1*(1.0);
    case 5
      ft=4.0*q2-1.0;		
    case 6
       ft=-4.0*(q1+2.0*q2+q3-1.0);	
    case 7
       ft=4.0*q3*(-1.0);
    case 8
       ft=0.0; 
    case 9
       ft=4.0*(1.0)*q3;
    case 10
      ft=0.0;      
end

end
function fz=fze3d(q1,q2,q3,n) 
% shape function derivative wrt xi
		  
switch n
    case 1
      fz=((-2.0)*(1.0-q1-q2-q3))+...
          ((1.0-2.0*q1-2.0*q2-2.0*q3)*(-1.0));	
    case 2
       fz=4.0*q1*(-1.0);
    case 3
      fz=0.0;
    case 4
       fz=0.0;
    case 5
      fz=0.0;	
    case 6
       fz=4.0*q2*(-1.0);	
    case 7
       fz=(4.0*q3*(-1.0))+(4.0*(1.0-q1-q2-q3));
    case 8
       fz=4.0*q1;
    case 9
       fz=4.0*q2;
    case 10
      fz=((2.0)*q3)+((2.0*q3-1.0)*1.0);	       
    end
end

  
