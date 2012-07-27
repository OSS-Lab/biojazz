function score = MatlabScoringScript(t,input,output) % A scoring function that assumes a step function stimulus

    delta_e = 0;
    xbase = input';
    ybase = output';
    ycount = 0;
    y = [0];
    last_i = -1;
    for i=1:length(xbase)
        if xbase(i)>5.0 && (xbase(i-250) >= 9.9 && xbase(i+1) >= 9.9)
            if (last_i ~= -1) && (last_i + 1 ~= i)
                disp('ERROR: something wrong with slice');
            end
            last_i = i;
            ycount=ycount+1;
            y(ycount) = ybase(i);
        end
    end
    delta_y = y(2:end) - y(1:end-1);
    turn_pts = [0];
    count = 0;
    for i = 2:length(delta_y)
        if (delta_y(i)<0) && (delta_y(i-1)>=0)
            count = count+1;
            turn_pts(count) = y(i);
        end
        if (delta_y(i)>=0) && (delta_y(i-1)<0)
            count = count+1;
            turn_pts(count) = y(i);
        end
    end
    if length(turn_pts)>1
        delta_turn_pts = abs(turn_pts(2:end)-turn_pts(1:end-1));
    else
        delta_turn_pts = [0];
    end
    
    delta_tp_avg = mean(delta_turn_pts);
    
    yaverage = mean(turn_pts);
    ydesired = zeros(length(turn_pts));
    turn_pts
%     delta_a = zeros(length(turn_pts));
%     peak1=turn_pts(1);
%     peak2=turn_pts(2);
%     for i=1:length(turn_pts)
%         if mod(i,2)==1
%             ydesired(i) = peak1;
%         else
%             ydesired(i) = peak2;
%         end
%         if (peak1>=peak2 && mod(i,2)==1) || (peak1 < peak2 && mod(i,2)==0)
%             delta_e = delta_e+(ydesired(i)-turn_pts(i));
%         else
%             delta_e = delta_e-(ydesired(i)-turn_pts(i));
%         end
%         delta_a(i)=abs(yaverage-turn_pts(i));
%     end
    
    %delta_a_avg = mean(delta_a);
    
    size(y)
    %y
    figure(800); plot(xbase);
    figure(801); plot(y)
    %delta_e_avg = delta_e/length(turn_pts);
    score = 1000*delta_tp_avg - var(delta_turn_pts)