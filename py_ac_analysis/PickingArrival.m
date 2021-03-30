
picking = 'AIC';

B = fir1(256,0.1,'low');

for i = 8:8
    
    waveform = p4976_run1_ACdata_test(2048*(i-1)+1:2048*i);
    
    output = filtfilt(B,1,waveform);
    
    output(1:120) = 0;

    plot(waveform + 3000*(i-1)); 
    hold on; plot(output + 3000*(i-1));


 dcmObj = datacursormode;
 set(dcmObj,'UpdateFcn',@GoodCursor);
 
 AIC = zeros(1,length(waveform));
 
 switch picking
     
     case {'AIC'}
     
         for j=1:(length(output)-1)
             
                 AIC(j) = j.*log10(var(output(1:j)))+(length(output)-j-1).*log10(var(output(j+1:length(output))));
         
         end
         
         AIC(AIC == -Inf) = max(AIC); 
         AIC(length(output)) = max(AIC);

         [~, AIC_idx]= min(AIC);

         plot(AIC + 3000*(i-1));  
         plot(AIC_idx, output(AIC_idx)+3000*(i-1), 'k*');

     end

 
end