function schmid_tensors = CalcSchmidTensors(rotations,phase_schmid_tensors,phaseofgrain,maxss)
%
%  Calculates the Schmid tensor for each grain in global coordinates
% 
%  converts to vector (6x1) form for tensor
%
  
  
  [mo,no,numgrains] = size(rotations);
  myrot = zeros(3,3);
  schmid_tensors = zeros(maxss,6,numgrains);
  
  for ngr = 1:numgrains

 myrot = rotations(:,:,ngr);
 
 phasenum = phaseofgrain(ngr);
 
 schmid_tensor(:,:,:) = phase_schmid_tensors(:,:,:,phasenum);
 
 sx_temp1 = zeros(maxss,3,3);
 sx_temp2 = zeros(maxss,3,3);
% 
% 
 for l=1:3
     for m=1:3
        sx_temp1(:,:,l) =  sx_temp1(:,:,l) + myrot(m,l)*schmid_tensor(:,:,m);
     end
 end

 for l=1:3
     for m=1:3
        sx_temp2(:,l,:) =  sx_temp2(:,l,:) + myrot(m,l)*sx_temp1(:,m,:);
     end
 end

schmid_tensors(:,1,ngr) = sx_temp2(:,1,1);
schmid_tensors(:,2,ngr) = sx_temp2(:,2,2);
schmid_tensors(:,3,ngr) = sx_temp2(:,3,3);
schmid_tensors(:,4,ngr) = 2.0*sx_temp2(:,1,2);
schmid_tensors(:,5,ngr) = 2.0*sx_temp2(:,1,3);
schmid_tensors(:,6,ngr) = 2.0*sx_temp2(:,2,3);

    
  end
  
end

