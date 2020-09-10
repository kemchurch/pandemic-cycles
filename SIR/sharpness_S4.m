function pass = sharpness_S4(p0,p1,para)
w = intval(para(6));
p = p0 + infsup(0,1)*(p1-p0);
if p(2)-w>p(1) & p(1)>w
    pass = 1;
else
    pass = 0;
end
end