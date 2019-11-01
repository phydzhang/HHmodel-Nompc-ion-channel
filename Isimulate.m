function   output=Isimulate(mgig,hgig,m0,h0,k,n,tm,th,Imax,t)

m = mgig-(mgig-m0).*exp(-t/tm);
h = hgig-(hgig-h0).*exp(-t/th);
output = -Imax.*m.^k.*h.^n-0;

end