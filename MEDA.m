function MEDA


for w = 54:1:56
    for x = 9.9:0.1:10.1
        for y = 29.5:0.5:30.5            
            
            h_max = w-(x+y);
            
            h = h_max - 0.5;
            h_min = h_max - 1;
            

            
            if(w-x-y-h_min)<=0.7
                h
            end
            
        end
    end
end


end