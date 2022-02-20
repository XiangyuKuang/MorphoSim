function shape = ShapeInit(EggshellPara,L,xx,yy,zz)
o = L/2;
shape = EggshellPara(1)*(1-sqrt((xx-o(1)).^2/EggshellPara(1).^2+(yy-o(2)).^2/EggshellPara(2)^2+(zz-o(3)).^2/EggshellPara(3)^2));
if EggshellPara(4)>0
  deform = EggshellPara(4)-abs(zz-o(3));
  shape(deform<shape) = deform(deform<shape);
end