function [sx, diffs] = spike_filter(sx)
        
      threshold = 1e18*0.05;
      r = size(sx, 1);
      lim = 50;
      cnt = 0;
      diffs = zeros(1, size(sx, 2));
      for c = 1:size(sx, 2)
              if immse(sum(sum(abs(sx(1:r/2 - lim,c)))), sum(sum(abs(sx(end - r/2 + lim : end, c))))) < threshold
                       sx(1:r/2 - lim, c) = 0;
                       sx(end - r/2+ lim : end , c) = 0;
                       cnt  = cnt + 1;
              end
              diffs(c) = immse(sum(sum(abs(sx(1:r/2 - lim,c)))), sum(sum(abs(sx(end - r/2 + lim : end, c)))));
      end              
end