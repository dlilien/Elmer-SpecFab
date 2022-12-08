lc=250.00;
p1=newp; Point(p1)={-20000.000, 0.000, 0.0, lc};
p2=newp; Point(p2)={20000.000, 0.000, 0.0, lc};
p3=newp; Point(p3)={20000.000, 8000.000, 0.0, lc};
p4=newp; Point(p4)={0.000, 8078.000, 0.0, lc};
p5=newp; Point(p5)={-20000.000, 8000.000, 0.0, lc};

s1=newreg; Line(s1) = {p1, p2};
s2=newreg; Line(s2) = {p2, p3};
s3=newreg; Line(s3) = {p3, p4, p5};
s4=newreg; Line(s4) = {p5, p1};

pl1=newreg; Physical Line(pl1)={s1};
pl2=newreg; Physical Line(pl2)={s2};
pl3=newreg; Physical Line(pl3)={s3};
pl4=newreg; Physical Line(pl4)={s4};
ll1=newreg; Line Loop(ll1)={s1,s2,s3,s4};
ps1=newreg; Plane Surface(ps1)={ll1};
ps2=newreg; Physical Surface(ps2)={ps1};
