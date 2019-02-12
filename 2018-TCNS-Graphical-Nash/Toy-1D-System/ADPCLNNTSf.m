function fi=ADPCLNNTSf(i,xi)

if i==0
    fi=0;
elseif i==1
fi = xi^2;
elseif i==2
fi = 0.5*xi^2;
elseif i==3
fi = xi^2 + 0.1*xi;
elseif i==4
fi = xi^2 + 0.5*xi;
else 
fi = xi^2 + 0.2*xi;
end