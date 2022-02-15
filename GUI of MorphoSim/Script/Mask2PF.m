function phi = Mask2PF(Mask)
filterSize = 9;
h = ones(filterSize,filterSize,filterSize)/filterSize^3;
blur = imfilter(Mask,h);
vOpt = sum(Mask,'all');
shift = fzero(@(x)Optfun(blur,x,vOpt,filterSize),0);
phi = (1+tanh((filterSize/2*(blur-0.5)+shift)*sqrt(2)))/2;
end

function loss = Optfun(blur,shift,v,a)
phi = (1+tanh((a/2*(blur-0.5)+shift)*sqrt(2)))/2;
loss = sum(phi,'all')-v;
end