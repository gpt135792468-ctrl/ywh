function p = NH3_psat(T)
  TK = T + 273.15;
  p = 10^(4.86886 - 1113.928/(TK - 10.409)) * 100;
end
