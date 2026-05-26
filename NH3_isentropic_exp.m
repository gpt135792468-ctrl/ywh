function [h3s, x] = NH3_isentropic_exp(s_in, T_out)
  sf = NH3_sf(T_out);  sg = NH3_sg(T_out);
  hf = NH3_hf(T_out);  hg = NH3_hg(T_out);
  x  = (s_in - sf) / (sg - sf);
  x  = max(0, min(1, x));
  h3s = hf + x*(hg - hf);
end
