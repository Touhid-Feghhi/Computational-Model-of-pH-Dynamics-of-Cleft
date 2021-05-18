t=0:10:1.5*10^6;
t=transpose(t);
[xx1, ia, ic] = unique(x1); 
yy1 = y1(ia);

[xx2, ia, ic] = unique(x2); 
yy2 = y2(ia);


[xx3, ia, ic] = unique(x3); 
yy3 = y3(ia);


[xx4, ia, ic] = unique(x4); 
yy4 = y4(ia);


[xx5, ia, ic] = unique(x5); 
yy5 = y5(ia);


[xx6, ia, ic] = unique(x6); 
yy6 = y6(ia);


[xx7, ia, ic] = unique(x7); 
yy7 = y7(ia);


[xx8, ia, ic] = unique(x8); 
yy8 = y8(ia);

[xx9, ia, ic] = unique(x9); 
yy9 = y9(ia);

pH1=pchip(xx1,yy1,t);
pH2=pchip(xx2,yy2,t);
pH3=pchip(xx3,yy3,t);
pH4=pchip(xx4,yy4,t);
pH5=pchip(xx5,yy5,t);
pH6=pchip(xx6,yy6,t);
pH7=pchip(xx7,yy7,t);
pH8=pchip(xx8,yy8,t);
pH9=pchip(xx9,yy9,t);
pH=1/9*(pH1+pH2+pH3+pH4+pH5+pH6+pH7+pH8+pH9);
plot(t,pH)
