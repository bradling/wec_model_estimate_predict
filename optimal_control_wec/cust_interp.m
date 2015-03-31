function y = cust_interp(ti, t, y, dt)

idx = find(floor(ti/dt)*dt == t);
if (ti-t(idx)) > 1e-8
    y = (y(idx+1) - y(idx)) * (ti-t(idx)) / (t(idx+1) - t(idx)) + y(idx);
else
    y = y(idx);
end