function L = safe_lmtd(dT1, dT2)
  if dT1 <= 0 || dT2 <= 0
    L = -1e3; return
  end
  if abs(dT1 - dT2) < 1e-8
    L = dT1;
  else
    L = (dT1 - dT2) / log(dT1/dT2);
  end
end
