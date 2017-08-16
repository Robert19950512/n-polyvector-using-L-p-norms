# N polyvector using L^p norms
===
This repository cintains materials which is used for testing how different norms affect curl of the field and the polyvector field quality based on different p value. the project is carry out under the supervision by professor [Daniele Panozzo](http://cs.nyu.edu/~panozzo/), and is based on the code and research from [Olga Diamanti](http://web.stanford.edu/~diamanti/). The test is carried out on summer 2016. some code and method is based on [Libigl](https://github.com/libigl/libigl/)
## introduction
directional field is wiedly used in many different computer graphics and geometric processing. it could be used to encode direction, value information on a surface or a point. the gradient information of a position can be represent as a vector, this could be consider as a  unit in directional field. According to the way it's used, there are many diiferent kinds of directional fields, as well as ways to generate them, PolyVector is one of those. In computer graphics and geometric processing, we name the generating process of directional field "synthesize". their have been significant developments on directional fields synthsis these years, which is cased by the increasing needed for its usage in for daiifernet reason, for instance: surface parametrization, mesh generation, texture synthesis, flow simulation, architectural geometry, fabrication and illustration. <br>
N-polyvector is one of these vector fields and its smoothness usually measured and optimized as the L2 norm between the coefficients of the polynomials which represent the field. This is done purely for efficiency reasons, but it is unclear how it relates to the curl of the field, a property which has more interesting applications. The idea of this project is to explore how different norms affect curl of the field and its singularity structure. We minimize Lp norms to rely on IRLS (iteratively reweighted least squares) which is only slightly more expensive than minimizing an L2 norm.

## how to use it
1.please follow the instruction and install [libigl](https://github.com/libigl/libigl/). <br>
2.replace the n-polyvector.cpp and n-polyvector.h in libigl using the one in this repository.
3. put the testing model(.obj) file into the tutorial/model folder. and revise the following code:
```cpp
b << 4550, 2321, 5413, 5350;

	

	// Load a mesh in OBJ format
	igl::readOBJ(TUTORIAL_SHARED_PATH "/aircraft.obj", V, F);
  ```

where b is the surface you want to manually insert vector fields, and  replace the aricraft.obj with your own tesing model.
4. run the main.cpp and see output.
## testing results
We use IRLS to minimize the Lp norm, but this arouse a problem, the length  of the vectorfield will become shorter and shorter after each iteration. currently we normalize the vectorfield after we generate it for visualization. It seems like the method we apply simply just minimize the energy by decrease it to nearly zero, and only the direction remain. so it’s not good for parametrization before this problem being solved. <br>
The code is only used for test and demonstration for Lpnorm with 2-polyvector field, when n is not equals to 2, the cold may crash.
## Q/A


