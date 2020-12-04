x = rand(1,1000);
a=0.1
M=4; 
for i = 1:length(x)
sum=0;
y1(i)=0;
for j = 1:M+1
    if (i-j>0)
     sum = sum + a*y1(i-j);
    end
end
y1(i)=sum +a*x(i);
end
t = 1:1:length(x);
plot(t,x);
figure
plot(t,y1);
