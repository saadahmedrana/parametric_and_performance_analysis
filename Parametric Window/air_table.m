function s=air_table(value)
%the function takes temperature in Kelvin as input
f=importfile('air.xlsx');
[r c]=size(f);%this gives us the size of the matrix containg the data imported from the sheet where r=no.of rows and c=no.of columns of the matrix
for i=1:r %the for loop compares the given value with those in the steam table
    if value>f(end,1)%if the value is greater than the greatest value in the table it gives an error
        error('input is out of bound');
    elseif value<f(1,1)%similarly if it is smaller than the smallest value again an error is generated
        error('input is out of bound');
    elseif value==f(i,1)%if value equals some value from the steam table, all the correspoding values are stored in the different varaible below
        h=f(i,2);
        u=f(i,3);
        s=f(i,4);
        pr=f(i,5);
        vr=f(i,6);
    elseif f(i,1)<value && f(i+1,1)>value %if the value is inbetween any 2 values from the table, the function interpolates all the data (interpolation is done from line 46 to 58)
        h=(((value-f(i,1))/(f(i+1,1)-f(i,1)))*(f(i+1,2)-f(i,2)))+f(i,2);
        u=(((value-f(i,1))/(f(i+1,1)-f(i,1)))*(f(i+1,3)-f(i,3)))+f(i,3);
        s=(((value-f(i,1))/(f(i+1,1)-f(i,1)))*(f(i+1,4)-f(i,4)))+f(i,4);
        pr=(((value-f(i,1))/(f(i+1,1)-f(i,1)))*(f(i+1,5)-f(i,5)))+f(i,5);
        vr=(((value-f(i,1))/(f(i+1,1)-f(i,1)))*(f(i+1,6)-f(i,6)))+f(i,6);
    end
end
    %fprintf('T = %f K\nh = %f kj/kg\nu = %f kj/kg\ns = %f kj/K\npr = %f\nvr = %f\n',value,h,u,s,pr,vr);
    %this last line prints all the properties and their values (all the values are in SI units)
end