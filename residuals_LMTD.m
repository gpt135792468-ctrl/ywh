function F = residuals_LMTD(x, Twwi,Tcwi,mww,mcw,mr,UAe,UAc,cpw,eta_t,eta_p)
    Teva = x(1);  Tcon = x(2);
    Twwo = x(3);  Tcwo = x(4);

    h2 = NH3_hg(Teva);
    s2 = NH3_sg(Teva);
    h4 = NH3_hf(Tcon);

    pe = NH3_psat(Teva);
    pc = NH3_psat(Tcon);
    vl = 1 / NH3_rhol(Tcon);
    wp = vl*(pe - pc)*1e3 / eta_p;
    h1 = h4 + wp/1000;

    [h3s, ~] = NH3_isentropic_exp(s2, Tcon);
    h3 = h2 - eta_t*(h2 - h3s);

    Qeva_r = mr*(h2 - h1)*1e3;
    Qeva_w = mww*cpw*(Twwi - Twwo);
    Qcon_r = mr*(h3 - h4)*1e3;
    Qcon_w = mcw*cpw*(Tcwo - Tcwi);

    LMTDe = safe_lmtd(Twwi - Teva, Twwo - Teva);
    LMTDc = safe_lmtd(Tcon - Tcwi, Tcon - Tcwo);

    F = zeros(4,1);
    F(1) = Qeva_r - Qeva_w;
    F(2) = Qeva_r - UAe*LMTDe;
    F(3) = Qcon_r - Qcon_w;
    F(4) = Qcon_r - UAc*LMTDc;
end
