/**
@cond

   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014-2016 MoDeNa Consortium, All rights reserved.

License
    This file is part of Modena.

    Modena is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    Modena is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.

@endcond
@file       Missing description
@todo       Document
@author     Jomas Mairhofer
@copyright  2014-2016, MoDeNa Project. GNU Public License.
@ingroup    twoTank
*/

#include <stdio.h>
#include <iostream>
#include "modena.h"

using namespace std;

int
main(int argc, char *argv[])
{
    double T = 270;
    double Tend = 290.0;
    
    // Instantiate index set
    //modena_index_set_t *indexSet = modena_index_set_new("species");


    // Instantiate a model
    modena_model_t *model = modena_model_new("SurfaceTension[A=AIR,B=THF]");   //muss das FunctionModule genau so heißen??
    if(modena_error_occurred())
    {
        return modena_error();
    }

    // Allocate memory and fetch arg positions
    modena_inputs_t *inputs = modena_inputs_new(model);           //How many inputs and outputs is defined in the function module!! 
    modena_outputs_t *outputs = modena_outputs_new(model);

    
    size_t Tpos = modena_model_inputs_argPos(model, "T");

    modena_model_argPos_check(model);

    
    while(T < Tend)
    {
      // Set input vector
      modena_inputs_set(inputs, Tpos, T);

      // Call the model
     int ret = modena_model_call(model, inputs, outputs);

     // Terminate, if requested
     if(modena_error_occurred())
            {
                modena_inputs_destroy(inputs);
                modena_outputs_destroy(outputs);
                modena_model_destroy(model);

                return modena_error();
            }

     // Fetch result
     double rho = modena_outputs_get(outputs, 0);
     
     cout << "T = " << T;
     
     
     T = T + 10.0;
    }    
    
       
    modena_inputs_destroy(inputs);
    modena_outputs_destroy(outputs);
    modena_model_destroy(model);

    return 0;
}
