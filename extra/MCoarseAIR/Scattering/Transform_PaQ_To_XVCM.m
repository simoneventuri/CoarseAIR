function [X, V] = Transform_PaQ_To_XVCM(PaQ, Masses)

    X(1,1:3) = PaQ(7:9);
    X(2,1:3) = PaQ(10:12);
    X(3,1:3) = - (Masses(1) .* X(1,1:3) + Masses(2) .* X(2,1:3)) / Masses(3);  

    V(1,1:3) = PaQ(1:3);
    V(2,1:3) = PaQ(4:6);
    V(3,1:3) = - (Masses(1) .* V(1,1:3) + Masses(2) .* V(2,1:3)) / Masses(3);
    
end