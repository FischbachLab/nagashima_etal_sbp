function xjits = AddJitterX(ythis,yres,xwidth)

xjits = zeros(length(ythis),1);

ymin = min(ythis);
ymax = max(ythis);
yarr = ymin:yres:(ceil((ymax-ymin)/yres)*yres+ymin);
for i = 1:(length(yarr)-1)
    if i>1
        inds_this = find(ythis>=yarr(i) & ythis<yarr(i+1));
    else
        inds_this = find(ythis>yarr(i) & ythis<yarr(i+1));
    end
    num_here = length(inds_this);
    if num_here>1
        total_spread = xwidth*0.9*(1-exp(log(0.9)*(num_here-1)));
        spread_bin = total_spread/(num_here-1);
        if mod(num_here,2)==0
            offset = spread_bin/2;
        else
            offset = eps;
        end
        for j = 1:num_here
            jthis = inds_this(j);
            xjits(jthis) = xjits(jthis) + offset;
            offset = offset - sign(offset)*spread_bin*j;
        end
    end
end

end