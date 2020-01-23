#
## A beam
#
algebraic3d


solid left = plane (0, 0, 0; -1, 0, 0) -bc=1;
solid right = plane(1, 0, 0;  1, 0, 0) -bc=2;

solid beam = left and right 
         and plane (0, 0, 0; 0, -1, 0)
         and plane (0, 0, 0; 0, 0, -1)
         and plane (1, 0.1, 0.1; 0, 0, 1)
         and plane (1, 0.1, 0.1; 0, 1, 0) -bc=3;

tlo beam;

