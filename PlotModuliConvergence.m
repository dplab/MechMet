function status = PlotModuliConvergence(residual_history,sx_moduli_history,sx_compliance_history)
%  Program for plotting the changes to sx moduli with iterations  
%
%    Input:   (1) history of residual with iterations
%             (2) history of the sx moduli with iterations
%
%     Output: (1) plots of input variables
%
  
    i=10;
    
 %    [numiterations,num_sx] = size(residual_history);
     [maxiterations,num_sx] = size(sx_moduli_history);
     
     numiterations = maxiterations;
     for iter = maxiterations-1:-1:1
      if(residual_history(iter)<=0)
          numiterations =iter-1;
      end
      end
     
    for iter = 1:1:numiterations
    abs_vals(iter) = iter;
    end
    
    maxresidual = 1.2*max(residual_history(:));
    
    maxmoduli = 1.5*max(sx_moduli_history(:,1));
   
   ab_label = 'Iteration';
   ord_label = 'Residual';
   plot_title = 'Strain Residual ';

   linetype =[ '-b' '--r' '--g' '--b' ':r' ':g' ':b' '--r' '--g' '--b' ];
   
    i=i+1;
    figure(i)    
    clf;    

    hold 'on';
 
    p = plot(abs_vals([1:numiterations]),residual_history([1:numiterations]),linetype(1));

    set(gca,'Box','on','LineWidth',2,'FontName','Helvetica','FontSize',14);
    set(p,'LineWidth',2)

    
    axis([0 numiterations 0 maxresidual]);
    xlabel(ab_label); 
    ylabel(ord_label);   
    title(plot_title); 
        
    hold 'off'
    
   i=i+1;       
   figure(i)      
   clf; 
        
   ab_label = 'Iteration';
   ord_label = 'Moduli Value';
   plot_title = 'Moduli Changes with Iteration ';

    hold 'on';
  
    for isx = 1:1:num_sx
        
        sx_vals([1:numiterations]) = sx_moduli_history([1:numiterations],isx);
        
        p = plot(abs_vals,sx_vals,linetype(isx));

        set(gca,'Box','on','LineWidth',2,'FontName','Helvetica','FontSize',14);
        set(p,'LineWidth',2)
    end

    
    axis([0 numiterations 0 maxmoduli]);
    xlabel(ab_label); 
    ylabel(ord_label);   
    title(plot_title); 
    legend('c11', 'c12', 'c33' );
    
    maxcompliance = 1.5*max(sx_compliance_history(:,1));
    
    mincompliance = 1.5*min(sx_compliance_history(:,2));
        
    hold 'off'
    
   i=i+1;       
   figure(i)      
   clf; 
        
   ab_label = 'Iteration';
   ord_label = 'Compliance Value';
   plot_title = 'Compliance Changes with Iteration ';

    hold 'on';
  
    for isx = 1:1:num_sx
        
        sx_vals([1:numiterations]) = sx_compliance_history([1:numiterations],isx);
        
        p = plot(abs_vals,sx_vals,linetype(isx));

        set(gca,'Box','on','LineWidth',2,'FontName','Helvetica','FontSize',14);
        set(p,'LineWidth',2)
    end

    
    axis([0 numiterations mincompliance maxcompliance]);
    xlabel(ab_label); 
    ylabel(ord_label);   
    title(plot_title); 
    legend('s11', 's12', 's33' );
        
    hold 'off'
    
    status = 'done with plot';
    
    


