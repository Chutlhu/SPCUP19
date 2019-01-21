function plot_structure()

    vertex_array = [...
            0.0615   -0.0615   -0.0615;
            0.0615    0.0615    0.0615;
            0.0615    0.0615   -0.0615;
           -0.0615    0.0615    0.0615;
           -0.0615    0.0615   -0.0615;
           -0.0615   -0.0615    0.0615;
           -0.0615   -0.0615   -0.0615;
            0.0615   -0.0615    0.0615];
    micPos = [...
            0.0420    0.0615   -0.0410;
           -0.0420    0.0615    0.0410;
           -0.0615    0.0420   -0.0410;
           -0.0615   -0.0420    0.0410;
           -0.0420   -0.0615   -0.0410;
            0.0420   -0.0615    0.0410;
            0.0615   -0.0420   -0.0410;
            0.0615    0.0420    0.0410];
  
    array_to_drone_base = [...
           0.0615, -0.05 ,0.0615;
           0.0615,  0.05 ,0.0615;
          -0.0385,  0.05 ,0.0615
          -0.0385, -0.05 ,0.0615;
           0.0615, -0.05 ,0.1915;
           0.0615,  0.05 ,0.1915;
          -0.0385,  0.05 ,0.1915
          -0.0385, -0.05 ,0.1915;];
        
    nMic = size(micPos,1);

    line = @(p1, p2, c) plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], c);
    
    %%% plot microphones
    pairs_connected = ...
        [ 1,3; 3,5; 5,7; 7,1;
          2,4; 4,6; 6,8; 8,2;
          1,8; 3,2; 5,4; 7,6];
    nEdges = size(pairs_connected,1);
    for i = 1:nMic
        scatter3(micPos(i,1), micPos(i,2), micPos(i,3));
        hold on
        text(micPos(i,1), micPos(i,2), micPos(i,3), num2str(i-1), 'FontSize', 15)
        scatter3(vertex_array(i,1), vertex_array(i,2), vertex_array(i,3));
    end
    
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    %%% plot microphone array structure
    for i = 1:nEdges
        line(vertex_array(pairs_connected(i,1),:), ...
             vertex_array(pairs_connected(i,2),:), 'b');
    end

    %%% plot drone structure
    % array support
    for i = 1:size(array_to_drone_base,1)
        scatter3(array_to_drone_base(i,1), array_to_drone_base(i,2), array_to_drone_base(i,3));
    end
    pairs_connected = ...
        [ 1,2; 2,3; 3,4; 4,1;
          5,6; 6,7; 7,8; 8,5;
          1,5; 2,6; 3,7; 4,8];
    for i = 1:nEdges
        line(array_to_drone_base(pairs_connected(i,1),:), ...
             array_to_drone_base(pairs_connected(i,2),:), 'r');
    end
    % rotors supports
    r = 0.485/2;
    rotors_pos = zeros(4,3);
    angle = 45:90:359;
    offset = [mean(array_to_drone_base(1:4,1),1), mean(array_to_drone_base(1:4,2),1), 0];
    for i = 1:size(rotors_pos,1)
       rotors_pos(i,:) = [r*cosd(angle(i) + 90), r*sind(angle(i) + 90), array_to_drone_base(end,3)] ...
                         + offset;
       scatter3(rotors_pos(i,1), rotors_pos(i,2), rotors_pos(i,3))
       text(rotors_pos(i,1), rotors_pos(i,2), rotors_pos(i,3), ...
            num2str(i), 'FontSize',15)
    end
    pairs_connected = [ 1,3; 2,4];
    for i = 1:size(pairs_connected,1)
        line(rotors_pos(pairs_connected(i,1),:), ...
             rotors_pos(pairs_connected(i,2),:), 'r');
    end
    
    hold off
end