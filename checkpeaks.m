function [pvf,ipf] = checkpeaks (a,pv,ip)

if (ip(4) > max(ip(1:3)) || ip(4) < min(ip(1:3))) && (ip(5) > max(ip(1:3)) || ip(5) < min(ip(1:3)))
    pvf = pv;
    ipf = ip;
else
    [pvt,ipt] = findpeaks(a,'Npeaks',20,'MinPeakDistance',5000,'SortStr','descend');
    pvf = pv(1:3);
    ipf = ip(1:3);
    for ii = 4:length(pvt)
        if ipt(ii) > max(ipf) || ipt(ii) < min(ipf)
            pvf = [pvf;pvt(ii)];
            ipf = [ipf;ipt(ii)];
        end
        if length(ipf) == 5
            break;
        end
    end
end
