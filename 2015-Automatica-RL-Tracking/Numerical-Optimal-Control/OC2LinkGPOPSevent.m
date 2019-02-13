function events = OC2LinkGPOPSevent(solevents)
e0 = solevents.initial.state;
ef = solevents.terminal.state;
x0 = [1.8 1.6 0 0];
xd0 = [1/2 1/3 0 0];
ei1 = e0 - (x0-xd0);
ei2 = ef;
events = [ei1;ei2];