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
    
    if(crystal_type==4)
       % if(covera==1)
          
    a=  1./sqrt(3.);
    b=  1./sqrt(2.);
    c=  1./sqrt(6.);
    d=  1./sqrt(14.);
    e=  1./sqrt(26.);

% {110}-type families    
    s_ss(1,1) = b;
    s_ss(1,2) = b;
    s_ss(1,3) = 0;
    m_ss(1,1) = -a;
    m_ss(1,2) = a;
    m_ss(1,3) = a;
    
    s_ss(2,1) = b;
    s_ss(2,2) = 0;
    s_ss(2,3) = b;
    m_ss(2,1) = -a;
    m_ss(2,2) = a;
    m_ss(2,3) = a;
    
    s_ss(3,1) = 0;
    s_ss(3,2) = -b;
    s_ss(3,3) = b;
    m_ss(3,1) = -a;
    m_ss(3,2) = a;
    m_ss(3,3) = a;
    
    s_ss(4,1) = b;
    s_ss(4,2) = -b;
    s_ss(4,3) = 0;
    m_ss(4,1) = a;
    m_ss(4,2) = a;
    m_ss(4,3) = a;
    
    s_ss(5,1) = b;
    s_ss(5,2) = 0;
    s_ss(5,3) = -b;
    m_ss(5,1) = a;
    m_ss(5,2) = a;
    m_ss(5,3) = a;
    
    s_ss(6,1) = 0;
    s_ss(6,2) = b;
    s_ss(6,3) = -b;
    m_ss(6,1) = a;
    m_ss(6,2) = a;
    m_ss(6,3) = a;
    
    s_ss(7,1) = b;
    s_ss(7,2) = -b;
    s_ss(7,3) = 0;
    m_ss(7,1) = a;
    m_ss(7,2) = a;
    m_ss(7,3) = -a;
    
    s_ss(8,1) = b;
    s_ss(8,2) = 0;
    s_ss(8,3) = b;
    m_ss(8,1) = a;
    m_ss(8,2) = a;
    m_ss(8,3) = -a;
    
    s_ss(9,1) = 0;
    s_ss(9,2) = b;
    s_ss(9,3) = b;
    m_ss(9,1) = a;
    m_ss(9,2) = a;
    m_ss(9,3) = -a;
    
    s_ss(10,1) = b;
    s_ss(10,2) = b;
    s_ss(10,3) = 0;
    m_ss(10,1) = a;
    m_ss(10,2) = -a;
    m_ss(10,3) = a;
    
    s_ss(11,1) = b;
    s_ss(11,2) = 0;
    s_ss(11,3) = -b;
    m_ss(11,1) = a;
    m_ss(11,2) = -a;
    m_ss(11,3) = a;
    
    s_ss(12,1) = 0;
    s_ss(12,2) = b;
    s_ss(12,3) = b;
    m_ss(12,1) = a;
    m_ss(12,2) = -a;
    m_ss(12,3) = a;
    
% {112}-type families    
    s_ss(13,1) = c;
    s_ss(13,2) = 2.*c;
    s_ss(13,3) = -c;
    m_ss(13,1) = -a;
    m_ss(13,2) = a;
    m_ss(13,3) = a;
    
    s_ss(14,1) = 2.*c;
    s_ss(14,2) = c;
    s_ss(14,3) = c;
    m_ss(14,1) = -a;
    m_ss(14,2) = a;
    m_ss(14,3) = a;
    
    s_ss(15,1) = c;
    s_ss(15,2) = -c;
    s_ss(15,3) = 2.*c;
    m_ss(15,1) = -a;
    m_ss(15,2) = a;
    m_ss(15,3) = a;
    
    s_ss(16,1) = c;
    s_ss(16,2) = -2.*c;
    s_ss(16,3) = c;
    m_ss(16,1) = a;
    m_ss(16,2) = a;
    m_ss(16,3) = a;
    
    s_ss(17,1) = -2.*c;
    s_ss(17,2) = c;
    s_ss(17,3) = c;
    m_ss(17,1) = a;
    m_ss(17,2) = a;
    m_ss(17,3) = a;
    
    s_ss(18,1) = c;
    s_ss(18,2) = c;
    s_ss(18,3) = -2.*c;
    m_ss(18,1) = a;
    m_ss(18,2) = a;
    m_ss(18,3) = a;
    
    s_ss(19,1) = -c;
    s_ss(19,2) = 2.*c;
    s_ss(19,3) = c;
    m_ss(19,1) = a;
    m_ss(19,2) = a;
    m_ss(19,3) = -a;
    
    s_ss(20,1) = -2.*c;
    s_ss(20,2) = c;
    s_ss(20,3) = -c;
    m_ss(20,1) = a;
    m_ss(20,2) = a;
    m_ss(20,3) = -a;
    
    s_ss(21,1) = c;
    s_ss(21,2) = c;
    s_ss(21,3) = 2.*c;
    m_ss(21,1) = a;
    m_ss(21,2) = a;
    m_ss(21,3) = -a;
    
    s_ss(22,1) = c;
    s_ss(22,2) = 2.*c;
    s_ss(22,3) = c;
    m_ss(22,1) = a;
    m_ss(22,2) = -a;
    m_ss(22,3) = a;
    
    s_ss(23,1) = 2.*c;
    s_ss(23,2) = c;
    s_ss(23,3) = -c;
    m_ss(23,1) = a;
    m_ss(23,2) = -a;
    m_ss(23,3) = a;
    
    s_ss(24,1) = -c;
    s_ss(24,2) = c;
    s_ss(24,3) = 2.*c;
    m_ss(24,1) = a;
    m_ss(24,2) = -a;
    m_ss(24,3) = a;
    
% {123}-type families    
    s_ss(25,1) = 3.*d;
    s_ss(25,2) = 2.*d;
    s_ss(25,3) = d;
    m_ss(25,1) = -a;
    m_ss(25,2) = a;
    m_ss(25,3) = a;
    
    s_ss(26,1) = 3.*d;
    s_ss(26,2) = 1.*d;
    s_ss(26,3) = 2.*d;
    m_ss(26,1) = -a;
    m_ss(26,2) = a;
    m_ss(26,3) = a;
    
    s_ss(27,1) = d;
    s_ss(27,2) = 3.*d;
    s_ss(27,3) = -2.*d;
    m_ss(27,1) = -a;
    m_ss(27,2) = a;
    m_ss(27,3) = a;
    
    s_ss(28,1) = 2.*d;
    s_ss(28,2) = 3.*d;
    s_ss(28,3) = -d;
    m_ss(28,1) = -a;
    m_ss(28,2) = a;
    m_ss(28,3) = a;
    
    s_ss(29,1) = d;
    s_ss(29,2) = -2.*d;
    s_ss(29,3) = 3.*d;
    m_ss(29,1) = -a;
    m_ss(29,2) = a;
    m_ss(29,3) = a;
    
    s_ss(30,1) = 2.*d;
    s_ss(30,2) = -d;
    s_ss(30,3) = 3.*d;
    m_ss(30,1) = -a;
    m_ss(30,2) = a;
    m_ss(30,3) = a;
    
    s_ss(31,1) = -3.*d;
    s_ss(31,2) = 2.*d;
    s_ss(31,3) = d;
    m_ss(31,1) = a;
    m_ss(31,2) = a;
    m_ss(31,3) = a;
    
    s_ss(32,1) = -3.*d;
    s_ss(32,2) = d;
    s_ss(32,3) = 2.*d;
    m_ss(32,1) = a;
    m_ss(32,2) = a;
    m_ss(32,3) = a;
    
    s_ss(33,1) = d;
    s_ss(33,2) = -3.*d;
    s_ss(33,3) = 2.*d;
    m_ss(33,1) = a;
    m_ss(33,2) = a;
    m_ss(33,3) = a;
    
    s_ss(34,1) = 2.*d;
    s_ss(34,2) = -3.*d;
    s_ss(34,3) = d;
    m_ss(34,1) = a;
    m_ss(34,2) = a;
    m_ss(34,3) = a;
    
    s_ss(35,1) = d;
    s_ss(35,2) = 2.*d;
    s_ss(35,3) = -3.*d;
    m_ss(35,1) = a;
    m_ss(35,2) = a;
    m_ss(35,3) = a;
    
    s_ss(36,1) = 2.*d;
    s_ss(36,2) = d;
    s_ss(36,3) = -3.*d;
    m_ss(36,1) = a;
    m_ss(36,2) = a;
    m_ss(36,3) = a;
    
    s_ss(37,1) = 3.*d;
    s_ss(37,2) = -2.*d;
    s_ss(37,3) = d;
    m_ss(37,1) = a;
    m_ss(37,2) = a;
    m_ss(37,3) = -a;
    
    s_ss(38,1) = 3.*d;
    s_ss(38,2) = -d;
    s_ss(38,3) = 2.*d;
    m_ss(38,1) = a;
    m_ss(38,2) = a;
    m_ss(38,3) = -a;
    
    s_ss(39,1) = -d;
    s_ss(39,2) = 3.*d;
    s_ss(39,3) = 2.*d;
    m_ss(39,1) = a;
    m_ss(39,2) = a;
    m_ss(39,3) = -a;
    
    s_ss(40,1) = -2.*d;
    s_ss(40,2) = 3.*d;
    s_ss(40,3) = d;
    m_ss(40,1) = a;
    m_ss(40,2) = a;
    m_ss(40,3) = -a;
    
    s_ss(41,1) = d;
    s_ss(41,2) = 2.*d;
    s_ss(41,3) = 3.*d;
    m_ss(41,1) = a;
    m_ss(41,2) = a;
    m_ss(41,3) = -a;
    
    s_ss(42,1) = 2.*d;
    s_ss(42,2) = d;
    s_ss(42,3) = 3.*d;
    m_ss(42,1) = a;
    m_ss(42,2) = a;
    m_ss(42,3) = -a;
    
    s_ss(43,1) = 3.*d;
    s_ss(43,2) = 2.*d;
    s_ss(43,3) = -d;
    m_ss(43,1) = a;
    m_ss(43,2) = -a;
    m_ss(43,3) = a;
    
    s_ss(44,1) = 3.*d;
    s_ss(44,2) = d;
    s_ss(44,3) = -2.*d;
    m_ss(44,1) = a;
    m_ss(44,2) = -a;
    m_ss(44,3) = a;
    
    s_ss(45,1) = d;
    s_ss(45,2) = 3.*d;
    s_ss(45,3) = 2.*d;
    m_ss(45,1) = a;
    m_ss(45,2) = -a;
    m_ss(45,3) = a;
    
    s_ss(46,1) = 2.*d;
    s_ss(46,2) = 3.*d;
    s_ss(46,3) = d;
    m_ss(46,1) = a;
    m_ss(46,2) = -a;
    m_ss(46,3) = a;
    
    s_ss(47,1) = -d;
    s_ss(47,2) = 2.*d;
    s_ss(47,3) = 3.*d;
    m_ss(47,1) = a;
    m_ss(47,2) = -a;
    m_ss(47,3) = a;
    
    s_ss(48,1) = -2.*d;
    s_ss(48,2) = d;
    s_ss(48,3) = 3.*d;
    m_ss(48,1) = a;
    m_ss(48,2) = -a;
    m_ss(48,3) = a;
    
% {134}-type families

    s_ss(49,1) = 4.*e;
    s_ss(49,2) = 3.*e;
    s_ss(49,3) = e;
    m_ss(49,1) = -a;
    m_ss(49,2) = a;
    m_ss(49,3) = a;
    
    s_ss(50,1) = 4.*e;
    s_ss(50,2) = e;
    s_ss(50,3) = 3.*e;
    m_ss(50,1) = -a;
    m_ss(50,2) = a;
    m_ss(50,3) = a;
    
    s_ss(51,1) = 3.*e;
    s_ss(51,2) = 4.*e;
    s_ss(51,3) = -e;
    m_ss(51,1) = -a;
    m_ss(51,2) = a;
    m_ss(51,3) = a;
    
    s_ss(52,1) = e;
    s_ss(52,2) = 4.*e;
    s_ss(52,3) = -3.*e;
    m_ss(52,1) = -a;
    m_ss(52,2) = a;
    m_ss(52,3) = a;
    
    s_ss(53,1) = 3.*e;
    s_ss(53,2) = -e;
    s_ss(53,3) = 4.*e;
    m_ss(53,1) = -a;
    m_ss(53,2) = a;
    m_ss(53,3) = a;
    
    s_ss(54,1) = e;
    s_ss(54,2) = -3.*e;
    s_ss(54,3) = 4.*e;
    m_ss(54,1) = -a;
    m_ss(54,2) = a;
    m_ss(54,3) = a;
    
    s_ss(55,1) = -4.*e;
    s_ss(55,2) = 3.*e;
    s_ss(55,3) = d;
    m_ss(55,1) = a;
    m_ss(55,2) = a;
    m_ss(55,3) = a;
    
    s_ss(56,1) = -4.*e;
    s_ss(56,2) = e;
    s_ss(56,3) = 3.*e;
    m_ss(56,1) = a;
    m_ss(56,2) = a;
    m_ss(56,3) = a;
    
    s_ss(57,1) = 3.*e;
    s_ss(57,2) = -4.*e;
    s_ss(57,3) = e;
    m_ss(57,1) = a;
    m_ss(57,2) = a;
    m_ss(57,3) = a;
    
    s_ss(58,1) = e;
    s_ss(58,2) = -4.*e;
    s_ss(58,3) = 3.*e;
    m_ss(58,1) = a;
    m_ss(58,2) = a;
    m_ss(58,3) = a;
    
    s_ss(59,1) = 3.*e;
    s_ss(59,2) = e;
    s_ss(59,3) = -4.*e;
    m_ss(59,1) = a;
    m_ss(59,2) = a;
    m_ss(59,3) = a;
    
    s_ss(60,1) = e;
    s_ss(60,2) = 3.*e;
    s_ss(60,3) = -4.*e;
    m_ss(60,1) = a;
    m_ss(60,2) = a;
    m_ss(60,3) = a;
    
    s_ss(61,1) = 4.*e;
    s_ss(61,2) = -3.*e;
    s_ss(61,3) = e;
    m_ss(61,1) = a;
    m_ss(61,2) = a;
    m_ss(61,3) = -a;
    
    s_ss(62,1) = 4.*e;
    s_ss(62,2) = -e;
    s_ss(62,3) = 3.*e;
    m_ss(62,1) = a;
    m_ss(62,2) = a;
    m_ss(62,3) = -a;
    
    s_ss(63,1) = -3.*e;
    s_ss(63,2) = 4.*e;
    s_ss(63,3) = e;
    m_ss(63,1) = a;
    m_ss(63,2) = a;
    m_ss(63,3) = -a;
    
    s_ss(64,1) = -e;
    s_ss(64,2) = 4.*e;
    s_ss(64,3) = 3.*e;
    m_ss(64,1) = a;
    m_ss(64,2) = a;
    m_ss(64,3) = -a;
    
    s_ss(65,1) = 3.*e;
    s_ss(65,2) = e;
    s_ss(65,3) = 4.*e;
    m_ss(65,1) = a;
    m_ss(65,2) = a;
    m_ss(65,3) = -a;
    
    s_ss(66,1) = e;
    s_ss(66,2) = 3.*e;
    s_ss(66,3) = 4.*e;
    m_ss(66,1) = a;
    m_ss(66,2) = a;
    m_ss(66,3) = -a;
    
    s_ss(67,1) = 4.*e;
    s_ss(67,2) = 3.*e;
    s_ss(67,3) = -e;
    m_ss(67,1) = a;
    m_ss(67,2) = -a;
    m_ss(67,3) = a;
   
    s_ss(68,1) = 4.*e;
    s_ss(68,2) = e;
    s_ss(68,3) = -3.*e;
    m_ss(68,1) = a;
    m_ss(68,2) = -a;
    m_ss(68,3) = a;
    
    s_ss(69,1) = 3.*e;
    s_ss(69,2) = 4.*e;
    s_ss(69,3) = e;
    m_ss(69,1) = a;
    m_ss(69,2) = -a;
    m_ss(69,3) = a;
    
    s_ss(70,1) = e;
    s_ss(70,2) = 4.*e;
    s_ss(70,3) = 3.*e;
    m_ss(70,1) = a;
    m_ss(70,2) = -a;
    m_ss(70,3) = a;
    
    s_ss(71,1) = -3.*e;
    s_ss(71,2) = e;
    s_ss(71,3) = 4.*e;
    m_ss(71,1) = a;
    m_ss(71,2) = -a;
    m_ss(71,3) = a;
    
    s_ss(72,1) = -e;
    s_ss(72,2) = 3.*e;
    s_ss(72,3) = 4.*e;
    m_ss(72,1) = a;
    m_ss(72,2) = -a;
    m_ss(72,3) = a;
    
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


%          c+a type 1 for consistency with FEpX

       ann(1,1)=1.;
       ann(1,2)=0.;
       ann(1,3)=-1.;
       ann(1,4)=1.;

       abb(1,1)=-2.;
       abb(1,2)=1.;
       abb(1,3)=1.;
       abb(1,4)=3.;

       ann(2,1)=1.;
       ann(2,2)=0.;
       ann(2,3)=-1.;
       ann(2,4)=1.;

       abb(2,1)=-1;
       abb(2,2)=-1;
       abb(2,3)=2;
       abb(2,4)=3;

       ann(3,1)=0.;
       ann(3,2)=1.;
       ann(3,3)=-1.;
       ann(3,4)=1.;

       abb(3,1)=-1.;
       abb(3,2)=-1.;
       abb(3,3)=2.;
       abb(3,4)=3.;

       ann(4,1)=0.;
       ann(4,2)=1.;
       ann(4,3)=-1.;
       ann(4,4)=1.;

       abb(4,1)=1.;
       abb(4,2)=-2.;
       abb(4,3)=1.;
       abb(4,4)=3.;

       ann(5,1)=-1.;
       ann(5,2)=1.;
       ann(5,3)=0.;
       ann(5,4)=1.;

       abb(5,1)=1.;
       abb(5,2)=-2.;
       abb(5,3)=1.;
       abb(5,4)=3.;

       ann(6,1)=-1.;
       ann(6,2)=0.;
       ann(6,3)=-1.;
       ann(6,4)=1.;

       abb(6,1)=2.;
       abb(6,2)=-1.;
       abb(6,3)=-1.;
       abb(6,4)=3.;
% 7 ok
       ann(7,1)=-1.;
       ann(7,2)=0.;
       ann(7,3)=1.;
       ann(7,4)=1.;
 
       abb(7,1)=2.;
       abb(7,2)=-1.;
       abb(7,3)=-1.;
       abb(7,4)=3.;
% 8 ok
       ann(8,1)=-1.;
       ann(8,2)=0.;
       ann(8,3)=1.;
       ann(8,4)=1.;

       abb(8,1)=1.;
       abb(8,2)=1.;
       abb(8,3)=-2.;
       abb(8,4)=3.;
% 9 ok
       ann(9,1)=0.;
       ann(9,2)=-1.;
       ann(9,3)=1.;
       ann(9,4)=1.;

       abb(9,1)=1.;
       abb(9,2)=1.;
       abb(9,3)=-2.;
       abb(9,4)=3.;
%  10  ok
       ann(10,1)=0.;
       ann(10,2)=-1.;
       ann(10,3)=1.;
       ann(10,4)=1.;

       abb(10,1)=-1.;
       abb(10,2)=2.;
       abb(10,3)=-1.;
       abb(10,4)=3.;
%  11  ok
       ann(11,1)=1.;
       ann(11,2)=-1.;
       ann(11,3)=0.;
       ann(11,4)=1.;

       abb(11,1)=-1.;
       abb(11,2)=2.;
       abb(11,3)=-1.;
       abb(11,4)=3.;
%  12  ok
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



