function z=value(x,y)
global tau r w b na pp a beta e i
global sigma vold
    z=lininter(a,vold(:,y),x);
end