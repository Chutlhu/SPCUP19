function plot_structure()

    % coordinates of the verteces of the microphone array
    vertex_array = [...
            0.0615   -0.0615   -0.0615;
            0.0615    0.0615    0.0615;
            0.0615    0.0615   -0.0615;
           -0.0615    0.0615    0.0615;
           -0.0615    0.0615   -0.0615;
           -0.0615   -0.0615    0.0615;
           -0.0615   -0.0615   -0.0615;
            0.0615   -0.0615    0.0615];
        
    % coordinates of the microphone
    micPos = [...
            0.0420    0.0615   -0.0410;
           -0.0420    0.0615    0.0410;
           -0.0615    0.0420   -0.0410;
           -0.0615   -0.0420    0.0410;
           -0.0420   -0.0615   -0.0410;
            0.0420   -0.0615    0.0410;
            0.0615   -0.0420   -0.0410;
            0.0615    0.0420    0.0410];
  
    % coordinates of the drone core
    drone_core = [...
           0.0615, -0.05 ,0.0615;
           0.0615,  0.05 ,0.0615;
          -0.0385,  0.05 ,0.0615
          -0.0385, -0.05 ,0.0615;
           0.0615, -0.05 ,0.1915;
           0.0615,  0.05 ,0.1915;
          -0.0385,  0.05 ,0.1915
          -0.0385, -0.05 ,0.1915;];
        
    nMic = size(micPos,1);

    % function handle for plotting lines between two points
    line = @(p1, p2, c) plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], c);
    
    figure();
    
    %%% plot microphones
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
    pairs_connected = ...       % pair connetected with an edge
        [ 1,3; 3,5; 5,7; 7,1;
          2,4; 4,6; 6,8; 8,2;
          1,8; 3,2; 5,4; 7,6];
    nEdges = size(pairs_connected,1);
    for i = 1:nEdges
        line(vertex_array(pairs_connected(i,1),:), ...
             vertex_array(pairs_connected(i,2),:), 'b');
    end

    %%% plot drone structure
    % array core
    for i = 1:size(drone_core,1)
        scatter3(drone_core(i,1), drone_core(i,2), drone_core(i,3));
    end
    pairs_connected = ...
        [ 1,2; 2,3; 3,4; 4,1;
          5,6; 6,7; 7,8; 8,5;
          1,5; 2,6; 3,7; 4,8];
    for i = 1:nEdges
        line(drone_core(pairs_connected(i,1),:), ...
             drone_core(pairs_connected(i,2),:), 'r');
    end
    
    % drone rotors
    r = 0.485/2;
    rotorsPos = zeros(4,3);
    angle = 45:90:359;
    offset = [mean(drone_core(1:4,1),1), mean(drone_core(1:4,2),1), 0];
    for i = 1:size(rotorsPos,1)
       rotorsPos(i,:) = [r*cosd(angle(i) + 90), r*sind(angle(i) + 90), drone_core(end,3)] ...
                         + offset;
       scatter3(rotorsPos(i,1), rotorsPos(i,2), rotorsPos(i,3))
       text(rotorsPos(i,1), rotorsPos(i,2), rotorsPos(i,3), ...
            num2str(i), 'FontSize',15)
    end
    pairs_connected = [ 1,3; 2,4];
    for i = 1:size(pairs_connected,1)
        line(rotorsPos(pairs_connected(i,1),:), ...
             rotorsPos(pairs_connected(i,2),:), 'r');
    end
    
    hold off
end