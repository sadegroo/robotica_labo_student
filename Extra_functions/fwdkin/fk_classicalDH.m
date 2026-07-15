function T = fk_classicalDH(DH_table, theta)
% Forward kinematics with classical DH (helper for verification).
% theta: nx1 vector of joint angles (radians). Overrides DH_table(:,4).
    T = eye(4);
    for i = 1:size(DH_table,1)
        al = DH_table(i,1);
        ai = DH_table(i,2);
        di = DH_table(i,3);
        th = theta(i);
        ct = cos(th); st = sin(th);
        ca = cos(al); sa = sin(al);
        Ai = [ct, -st*ca,  st*sa,  ai*ct;
              st,  ct*ca, -ct*sa,  ai*st;
               0,     sa,     ca,     di;
               0,      0,      0,      1];
        T = T * Ai;
    end
end