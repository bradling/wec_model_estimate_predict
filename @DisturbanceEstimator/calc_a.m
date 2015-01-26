function [A, varargout] = calc_a(obj, estState, myCase)
% Note - estState vector is : [z dz f df (parms)]^T

switch myCase
    case '10'
        c1 = obj.k/obj.m;
        c3 = obj.feFreq^2;
        Ac = [ 1   0           0  0   0;
             -c1 -estState(5)  1  0 -estState(2);
               0   0           0  1   0;
               0   0         -c3  0   0; 
               0   0           0  0   0];
         A = expm(obj.dt .* Ac);
         Fc = [ 1   0           0  0  0;
             -c1 -estState(5)  1  0  0;
               0   0           0  1   0;
               0   0         -c3  0   0; 
               0   0           0  0   0];
          varargout{1} = expm(obj.dt .* Fc);
    case '01'
        c1 = obj.k/obj.m;
        c2 = obj.b/obj.m;
        Ac = [  1   0           0  0           0 ;
              -c1 -c2           1  0           0 ;
                0   0           0  1           0 ;
                0   0 -estState(5) 0 -estState(3); 
                0   0           0  0           0 ];
        A = expm(obj.dt .* Ac);
        Fc = [  1   0           0  0           0 ;
              -c1 -c2           1  0           0 ;
                0   0           0  1           0 ;
                0   0 -estState(5) 0 0; 
                0   0           0  0           0 ];
        varargout{1} = expm(obj.dt .* Fc);
    case '11'
        c1 = obj.k/obj.m;
        Ac = [ 1   0           0  0           0   0;
             -c1 -estState(6)  1  0           0 -estState(2);
               0   0           0  1           0   0;
               0   0 -estState(5) 0 -estState(3)  0; 
               0   0           0  0           0   0;
               0   0           0  0           0   0];
         A = expm(obj.dt .* Ac);
         Fc = [ 1   0           0  0           0   0;
             -c1 -estState(6)  1  0           0 0;
               0   0           0  1           0   0;
               0   0 -estState(5) 0 0  0; 
               0   0           0  0           0   0;
               0   0           0  0           0   0];
           varargout{1} = expm(obj.dt .* Fc);
    case '00'
        c1 = obj.k/obj.m;
        c2 = obj.b/obj.m;
        c3 = obj.feFreq^2;
        Ac = [  1   0   0 0 ;
              -c1 -c2   1 0 ;
                0   0   0 1 ;
                0   0 -c3 0 ];
        A = expm(obj.dt .* Ac);
end %switch

