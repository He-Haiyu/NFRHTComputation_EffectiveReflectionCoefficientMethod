function var=checksqrt(x)
var1 = sqrt(x);
% var1(imag(var1)>0)=-var1(imag(var1)>0);
for i=1:1:length(var1)
    if imag(var1(i))<0
     var1(i) = -var1(i);
    end
end
var = var1;
end