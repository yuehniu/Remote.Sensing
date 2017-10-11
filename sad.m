function out = sad(in1, in2)

    t1 = in1' * in2;
    t2 = norm(in1) * norm(in2);
    
    out = acosd(t1 / t2);
end