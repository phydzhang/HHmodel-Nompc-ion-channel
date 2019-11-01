function   output=Rmse(x,m0,h0,k,n,Imax,t,current_Data,i)
mgig = x(1);
hgig = x(2);

m = mgig-(mgig-m0).*exp(-t/x(3));
h = hgig-(hgig-h0).*exp(-t/x(4));
I = -Imax.*m.^k.*h.^n-1032;

Ireshape = reshape(I,2292,1);
output = sqrt(sum((Ireshape(:)-current_Data(200:2491,i)).^2)/2292);


end