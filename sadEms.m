function sadOut = sadEms(emTrue, emObj, emNum)

    sadOut = zeros(emNum, 1);

    for em_i = 1:emNum
            tmp_sad = inf;
            for em_j = 1:emNum
                cur_sad = sad(emTrue(em_j,:)', emObj(em_i,:)');
                if cur_sad<tmp_sad
                    tmp_sad = cur_sad;
                end     
            end
            if(tmp_sad == inf)
                tmp_sad = mean(sadOut(em_i)/emNum);
            end
            sadOut(em_i) = tmp_sad / 180 * pi;
    end
    
end