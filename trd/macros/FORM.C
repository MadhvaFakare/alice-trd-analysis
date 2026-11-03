void FORM(){
    double pTmin = 1.0;
    double pTmax = 3.0;
    int layer = 0;
    TString a = Form("fPt>%g && fPt<%g", pTmax, pTmin);
    TString b = Form("fQ1%d",layer);
    printf("a = %s\n",a.Data());
    printf("b = %s\n",b.Data());
}