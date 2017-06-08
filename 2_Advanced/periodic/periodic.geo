##
## Example with periodic boundary conditions
##    by Joachim Schoeberl
##
##

algebraic3d

solid bot = plane (0, 0, 0; 0, 0, -1) -bc=1;
solid top = plane (1, 1, 1; 0, 0, 1) -bc=1;

solid left = plane (0, 0, 0; 0, -1, 0) -bc=2;
solid right = plane (1, 1, 1; 0, 1, 0) -bc=2;

solid front = plane (0, 0, 0; -1, 0, 0) -bc=1;
solid back = plane (1, 1, 1; 1, 0, 0) -bc=1;
 
solid cube = bot and top and left and right and front and back;

solid ball = sphere(0.5, 0.7, 0.5; 0.2) -bc=3;
solid matrix = cube and not ball;

tlo matrix -transparent -material=matrix -maxh=0.2;
tlo ball -col=[1,0,0] -material=ball -maxh=0.2;

identify periodic left right;

