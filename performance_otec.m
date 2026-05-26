function [h1,h2,h3,h4,s1,s2,s3,s4,Wt,Wp] = performance_otec(Teva, Tcon, mr, eta_t, eta_p)
    h2 = NH3_hg(Teva);   s2 = NH3_sg(Teva);
    h4 = NH3_hf(Tcon);   s4 = NH3_sf(Tcon);

    [h3s, ~] = NH3_isentropic_exp(s2, Tcon);
    h3 = h2 - eta_t*(h2 - h3s);

    hf3 = NH3_hf(Tcon); hg3 = NH3_hg(Tcon);
    sf3 = NH3_sf(Tcon); sg3 = NH3_sg(Tcon);
    x3  = max(0, min(1, (h3-hf3)/(hg3-hf3)));
    s3  = sf3 + x3*(sg3 - sf3);

    pe = NH3_psat(Teva);  pc = NH3_psat(Tcon);
    vl = 1/NH3_rhol(Tcon);
    wp = vl*(pe - pc)*1e3 / eta_p;

    h1 = h4 + wp/1000;
    s1 = s4;

    Wt = mr*(h2 - h3)*1e3;
    Wp = mr*wp;
end
