%%Find the ripple margin in passband or stopband of a lowpass filter
function [p_max_ripple,s_max_ripple] = find_max_ripple(H,w,wp,ws)
    p_max_ripple = 0;
    s_max_ripple = 0;
    i=1;
    %Parse through H from 0 to passband, calculate ripples, and save the max
    while(w(i) < wp)
        ripple = abs(1-abs(H(i)));
        if ripple > p_max_ripple
            p_max_ripple = ripple;
        end
        i = i + 1;
    end
    %Skip the transition band
    while(w(i) < ws)
        i = i + 1;
    end
    %Parse through H from stopband to pi, calculate ripples, and save the max
    while(i < length(w))
        ripple = abs(H(i));
        if ripple > s_max_ripple
            s_max_ripple = ripple;
        end
        i = i + 1;
    end
    
end

