# Short description

This is the program of a mathematical model, used in paper "Modeling the Emergence and Evolution of Predation". <br>
To run this program you could build this project under IDE (for example CodeBlocks) or write makefile by yourself ))


# Files description
## C++ program files
**inout.cpp, inout.h** - file with input/output functionality  
**args.cpp, args.h** - classes declaration, constructors, etc.  
**clustering.cpp, clustering.h** - functionality, connected with populations clustering  
**mfunctions.cpp, mfunctions.h** - **M**odel functions, general functions, used in the model  
**selective_grad.cpp, selective_grad.h** - functionality, connected with selective gradient dynamics
## CONFIG FILES
**model_params.txt** - file with model parameters, has to be located in the same folder as exe file (or you could modify the pass in main)  
**populationz/0.txt** - initial population  
**res_dir/0.txt** - initial resourse density among the domain (it is a fucntion z=f(x,y), (x,y) in [-3,3]x[-3,3] square)  
## PLOT FILES
**evolution_visualization.py** - visualization for populations, you will probably have to modify the pass to the resulting files inside the program  
**resource_dense_XD** - same but for resource density
## ADDITIONAL
**video_res...** - resource density dynamic example
**video_pop...** - population dynamic example
