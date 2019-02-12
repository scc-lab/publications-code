function D = OCGenD(OrderOfCollocation)

[tau w]=OClgwt(OrderOfCollocation,-1,1);

tau=[-1; tau];
tau=sort(tau);
D = OCcollocD(tau);
D = D((2:end),:);

