function output=calRmse(mgig,hgig,m0,h0,tm,th,k,n,Imax,t,current_Data,i)

m = mgig-(mgig-m0).*exp(-t/tm);
h = hgig-(hgig-h0).*exp(-t/th);
I = -Imax.*m.^k.*h.^n-1032;

Ireshape = reshape(I,2292,1);
output = sqrt(sum((Ireshape(:)-current_Data(200:2491,i)).^2)/2292);


end