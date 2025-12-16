#include <iostream>
#include <stdio.h>
#include <stdlib.h>


#include "ELEMENTS.h"


#include "mesh_inputs.h"
#include "mesh_io.h"
#include "state.h"


int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope
    std::cout<<"Hello, Average Example!"<<std::endl;
    
    
    
    
    
} // end MATAR scope
MATAR_FINALIZE();

return 0;
}