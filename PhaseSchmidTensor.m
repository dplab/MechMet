function phase_schmid_tensors = PhaseSchmidTensor(nphases,crystal_type_all,covera_all,maxss)
%
%  Calculates the Schmid tensor for either cubic (3) or hexagonal (6)
%  symmetry
%
%   return a 3x3 tensor for all ss of all grains, which is later convertedt
%   to a 6x1 tensor
%
phase_schmid_tensors = zeros(maxss,3,3,nphases);

for iphase = 1:nphases 
    
 crystal_type = crystal_type_all(iphase);
% 
%  for cubics, covera is the flag for fcc or bcc type ss
%  for hexagonals, covera is the c over a ratio
%
 covera = covera_all(iphase);

    num_ss = maxss;    
    m_ss = zeros(num_ss,3);
    s_ss = zeros(num_ss,3);
      
    if(crystal_type==3)
        
        if(covera==1)
        
    a=  1./sqrt(3.);
    b=  1./sqrt(2.);
  
    s_ss(1,1) = 0;
    s_ss(1,2) = b;
    s_ss(1,3) = -b;
    m_ss(1,1) = a;
    m_ss(1,2) = a;
    m_ss(1,3) = a;
    
    s_ss(2,1) = -b;
    s_ss(2,2) = 0;
    s_ss(2,3) = b;
    m_ss(2,1) = a;
    m_ss(2,2) = a;
    m_ss(2,3) = a;    
    
    s_ss(3,1) = b;
    s_ss(3,2) = -b;
    s_ss(3,3) = 0;
    m_ss(3,1) = a;
    m_ss(3,2) = a;
    m_ss(3,3) = a;
        
    
    s_ss(4,1) = 0;
    s_ss(4,2) = b;
    s_ss(4,3) = -b;
    m_ss(4,1) = -a;
    m_ss(4,2) = a;
    m_ss(4,3) = a;
        
    
    s_ss(5,1) = b;
    s_ss(5,2) = 0;
    s_ss(5,3) = b;
    m_ss(5,1) = -a;
    m_ss(5,2) = a;
    m_ss(5,3) = a;
        
    
    s_ss(6,1) = b;
    s_ss(6,2) = b;
    s_ss(6,3) = 0;
    m_ss(6,1) = -a;
    m_ss(6,2) = a;
    m_ss(6,3) = a;
        
    
    s_ss(7,1) = 0;
    s_ss(7,2) = b;
    s_ss(7,3) = b;
    m_ss(7,1) = -a;
    m_ss(7,2) = -a;
    m_ss(7,3) = a;
        
    
    s_ss(8,1) = b;
    s_ss(8,2) = 0;
    s_ss(8,3) = b;
    m_ss(8,1) = -a;
    m_ss(8,2) = -a;
    m_ss(8,3) = a;
        
    s_ss(9,1) = -b;
    s_ss(9,2) = b;
    s_ss(9,3) = 0;
    m_ss(9,1) = -a;
    m_ss(9,2) = -a;
    m_ss(9,3) = a;
        
    s_ss(10,1) = 0;
    s_ss(10,2) = b;
    s_ss(10,3) = b;
    m_ss(10,1) = a;
    m_ss(10,2) = -a;
    m_ss(10,3) = a;
        
    s_ss(11,1) = b;
    s_ss(11,2) = 0;
    s_ss(11,3) = -b;
    m_ss(11,1) = a;
    m_ss(11,2) = -a;
    m_ss(11,3) = a;
        
    s_ss(12,1) = b;
    s_ss(12,2) = b;
    s_ss(12,3) = 0;
    m_ss(12,1) = a;
    m_ss(12,2) = -a;
    m_ss(12,3) = a;
    
    elseif(covera==2)
            
    a=  1./sqrt(3.);
    b=  1./sqrt(2.);
    
    s_ss(1,1) = 0;
    s_ss(1,2) = b;
    s_ss(1,3) = -b;
    m_ss(1,1) = a;
    m_ss(1,2) = a;
    m_ss(1,3) = a;
    
    s_ss(2,1) = -b;
    s_ss(2,2) = 0;
    s_ss(2,3) = b;
    m_ss(2,1) = a;
    m_ss(2,2) = a;
    m_ss(2,3) = a;    
    
    s_ss(3,1) = b;
    s_ss(3,2) = -b;
    s_ss(3,3) = 0;
    m_ss(3,1) = a;
    m_ss(3,2) = a;
    m_ss(3,3) = a;
        
    
    s_ss(4,1) = 0;
    s_ss(4,2) = b;
    s_ss(4,3) = -b;
    m_ss(4,1) = -a;
    m_ss(4,2) = a;
    m_ss(4,3) = a;
        
    
    s_ss(5,1) = b;
    s_ss(5,2) = 0;
    s_ss(5,3) = b;
    m_ss(5,1) = -a;
    m_ss(5,2) = a;
    m_ss(5,3) = a;
        
    
    s_ss(6,1) = b;
    s_ss(6,2) = b;
    s_ss(6,3) = 0;
    m_ss(6,1) = -a;
    m_ss(6,2) = a;
    m_ss(6,3) = a;
        
    
    s_ss(7,1) = 0;
    s_ss(7,2) = b;
    s_ss(7,3) = b;
    m_ss(7,1) = -a;
    m_ss(7,2) = -a;
    m_ss(7,3) = a;
        
    
    s_ss(8,1) = b;
    s_ss(8,2) = 0;
    s_ss(8,3) = b;
    m_ss(8,1) = -a;
    m_ss(8,2) = -a;
    m_ss(8,3) = a;
        
    s_ss(9,1) = -b;
    s_ss(9,2) = b;
    s_ss(9,3) = 0;
    m_ss(9,1) = -a;
    m_ss(9,2) = -a;
    m_ss(9,3) = a;
        
    s_ss(10,1) = 0;
    s_ss(10,2) = b;
    s_ss(10,3) = b;
    m_ss(10,1) = a;
    m_ss(10,2) = -a;
    m_ss(10,3) = a;
        
    s_ss(11,1) = b;
    s_ss(11,2) = 0;
    s_ss(11,3) = -b;
    m_ss(11,1) = a;
    m_ss(11,2) = -a;
    m_ss(11,3) = a;
        
    s_ss(12,1) = b;
    s_ss(12,2) = b;
    s_ss(12,3) = 0;
    m_ss(12,1) = a;
    m_ss(12,2) = -a;
    m_ss(12,3) = a;
        end
            
  end  
    
    if(crystal_type==6)

    num_ss_bp = 6;
    num_ss_pyr = 12;
             
    a=  sqrt(3.)/2.;
    b=  1./2.;
    sqr3 = sqrt(3.);

      m_ss(1,1) = 0.0;
      m_ss(1,2) = 0.0;
      m_ss(1,3) = 1.0;
      m_ss(2,1) = 0.0;
      m_ss(2,2) = 0.0;
      m_ss(2,3) = 1.0;
      m_ss(3,1) = 0.0;
      m_ss(3,2) = 0.0;
      m_ss(3,3) = 1.0;
      m_ss(4,1) = 0.0;
      m_ss(4,2) = 1.0;
      m_ss(4,3) = 0.0;
      m_ss(5,1) = -a;
      m_ss(5,2) = -b;
      m_ss(5,3) = 0.0;
      m_ss(6,1) = a;
      m_ss(6,2) = -b;
      m_ss(6,3) = 0.0;
      
      s_ss(1,1) = 1.0;
      s_ss(1,2) = 0.0;
      s_ss(1,3) = 0.0;
      s_ss(2,1) = -b;
      s_ss(2,2) = a;
      s_ss(2,3) = 0.0;
      s_ss(3,1) = -b;
      s_ss(3,2) = -a;
      s_ss(3,3) = 0.0;
      s_ss(4,1) = 1.0;
      s_ss(4,2) = 0.0;
      s_ss(4,3) = 0.0;
      s_ss(5,1) = -b;
      s_ss(5,2) = a;
      s_ss(5,3) = 0.0;
      s_ss(6,1) = -b;
      s_ss(6,2) = -a;
      s_ss(6,3) = 0.0;

%       pyramidal 11-22 -1-123

%          c+a type 2

       ann(1,1)=-2.;
       ann(1,2)=1.;
       ann(1,3)=1.;
       ann(1,4)=2.;

       abb(1,1)=2.;
       abb(1,2)=-1.;
       abb(1,3)=-1.;
       abb(1,4)=3.;

       ann(2,1)=1.;
       ann(2,2)=-2.;
       ann(2,3)=1.;
       ann(2,4)=2.;

       abb(2,1)=-1.;
       abb(2,2)=2.;
       abb(2,3)=-1.;
       abb(2,4)=3.;

       ann(3,1)=1.;
       ann(3,2)=1.;
       ann(3,3)=-2.;
       ann(3,4)=2.;

       abb(3,1)=-1.;
       abb(3,2)=-1.;
       abb(3,3)=2.;
       abb(3,4)=3.;

       abb(3,1)=-1.;
       abb(3,2)=-1.;
       abb(3,3)=2.;
       abb(3,4)=3.;

       ann(4,1)=2.;
       ann(4,2)=-1.;
       ann(4,3)=-1.;
       ann(4,4)=2.;

       abb(4,1)=-2.;
       abb(4,2)=1.;
       abb(4,3)=1.;
       abb(4,4)=3.;

       ann(5,1)=-1.;
       ann(5,2)=2.;
       ann(5,3)=-1.;
       ann(5,4)=2.;

       abb(5,1)=1.;
       abb(5,2)=-2.;
       abb(5,3)=1.;
       abb(5,4)=3.;

       ann(6,1)=-1.;
       ann(6,2)=-1.;
       ann(6,3)=2.;
       ann(6,4)=2.;

       abb(6,1)=1.;
       abb(6,2)=1.;
       abb(6,3)=-2.;
       abb(6,4)=3.;
       
       ann(7,1)=-1.;
       ann(7,2)=0.;
       ann(7,3)=1.;
       ann(7,4)=1.;

       abb(7,1)=2.;
       abb(7,2)=-1.;
       abb(7,3)=-1.;
       abb(7,4)=3.;

       ann(8,1)=-1.;
       ann(8,2)=0.;
       ann(8,3)=1.;
       ann(8,4)=1.;

       abb(8,1)=1.;
       abb(8,2)=1.;
       abb(8,3)=-2.;
       abb(8,4)=3.;

       ann(9,1)=0.;
       ann(9,2)=-1.;
       ann(9,3)=1.;
       ann(9,4)=1.;

       abb(9,1)=1.;
       abb(9,2)=1.;
       abb(9,3)=-2.;
       abb(9,4)=3.;

       ann(10,1)=0.;
       ann(10,2)=-1.;
       ann(10,3)=1.;
       ann(10,4)=1.;

       abb(10,1)=-1.;
       abb(10,2)=2.;
       abb(10,3)=-1.;
       abb(10,4)=3.;

       ann(11,1)=1.;
       ann(11,2)=-1.;
       ann(11,3)=0.;
       ann(11,4)=1.;

       abb(11,1)=-1.;
       abb(11,2)=2.;
       abb(11,3)=-1.;
       abb(11,4)=3.;

       ann(12,1)=1.;
       ann(12,2)=-1.;
       ann(12,3)=0.;
       ann(12,4)=1.;

       abb(12,1)=-2.;
       abb(12,2)=1.;
       abb(12,3)=1.;
       abb(12,4)=3.;


      nslips = num_ss_bp + 1;
      num_ss_all = num_ss_bp+num_ss_pyr;

        for  i=nslips:num_ss_all

        jcount = i - num_ss_bp;

        m_ss(i,1) = ann(jcount,1);
        m_ss(i,2) = (2.*ann(jcount,2)+ann(jcount,1))/sqr3;
        m_ss(i,3) = ann(jcount,4)/covera;

        xnrm = m_ss(i,1)^2 + m_ss(i,2)^2 + m_ss(i,3)^2;
        xnrm = sqrt(xnrm);
        m_ss(i,1) = m_ss(i,1) / xnrm;
        m_ss(i,2) = m_ss(i,2) / xnrm;
        m_ss(i,3) = m_ss(i,3) / xnrm;

        s_ss(i,1) = 3.*abb(jcount,1)/2.;
        s_ss(i,2) = (abb(jcount,1)/2. + abb(jcount,2))*sqr3;
        s_ss(i,3) = abb(jcount,4) * covera;

        xnrs = s_ss(i,1)^2 + s_ss(i,2)^2 + s_ss(i,3)^2;
        xnrs = sqrt(xnrs);
        s_ss(i,1) = s_ss(i,1) / xnrs;
        s_ss(i,2) = s_ss(i,2) / xnrs;
        s_ss(i,3) = s_ss(i,3) / xnrs;
        
        end

    end
    

    
ss_schmid=zeros(3,3);

for i=1:num_ss
    
    for  j=1:3
        for k=1:3
         ss_schmid(j,k) = s_ss(i,j)*m_ss(i,k);
        end
    end
    
    sym_schmid =  (ss_schmid+ss_schmid')/2;
    
        for  j=1:3
        for k=1:3
             phase_schmid_tensor(i,j,k) = sym_schmid(j,k);
        end
    end

    
end

phase_schmid_tensors(:,:,:,iphase) = phase_schmid_tensor(:,:,:);

end



