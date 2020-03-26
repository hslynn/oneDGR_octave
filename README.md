# An octave(matlab) program to solve Einstein Field Equations in 1 dimension (spherical reduction)

* Working parameters: 
    1. inB = 1.5, outB = 2.5;
    2. N = 15, meshNum = 1;
    3. paragamma1 = 10;
    4. bdry_type = FREEZING.

* Want to make working: 
    * outB >= 5;
    * bdry = FREEZING or DIRICHLET;
    * paragamma1?

* How to combine freezing incoming characteristic fields condition with filter?
    1. Adjusting rhs(vmapO) actually affect all the values of the outermost cell;
    2. If we apply freezing condition first, then applying the filter will introduce instabilities(known from test);
    3.  

