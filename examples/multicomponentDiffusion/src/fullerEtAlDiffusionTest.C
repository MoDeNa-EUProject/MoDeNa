/*

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

Description
    Solving the two tank problem the MoDeNa way.

    A prototypical macros-scopic code embeds a micro-scale model (flowRate)
    through the MoDeNa interface library.

Authors
    Henrik Rusche

Contributors
*/

#include <stdio.h>
#include <iostream>
#include "modena.h"

using namespace std;

int
main(int argc, char *argv[])
{
    {
        // Instantiate index set
        modena_index_set_t *indexSet = modena_index_set_new("species");

        // Find index
        cout << "index[N2] = "
             << modena_index_set_get_index(indexSet, "N2")
             << endl;

        // Iterate index set and print names
        size_t start = modena_index_set_iterator_start(indexSet);
        size_t end = modena_index_set_iterator_end(indexSet);
        for(int i = start; i < end; i++)
        {
            cout << "index[" << i << "] = "
                 << modena_index_set_get_name(indexSet, i)
                 << endl;
        }

        modena_index_set_destroy(indexSet);
    }

    {
        // Instantiate function
        modena_function_t *function =
            modena_function_new("fullerEtAlDiffusion");

        // Get index set from function
        modena_index_set_t *indexSet =
           modena_function_get_index_set(function, "A");

        // Find index
        cout << "index[N2] = "
             << modena_index_set_get_index(indexSet, "N2")
             << endl;

        // Iterate index set and print names
        size_t start = modena_index_set_iterator_start(indexSet);
        size_t end = modena_index_set_iterator_end(indexSet);
        for(int i = start; i < end; i++)
        {
            cout << "index[" << i << "] = "
                 << modena_index_set_get_name(indexSet, i)
                 << endl;
        }

        modena_index_set_destroy(indexSet);
        // Not working - BUG!
        //modena_function_destroy(function);
    }

    // Instantiate model
    modena_model_t *model = modena_model_new
    (
        "fullerEtAlDiffusion[A=H2O,B=N2]"
    );

    if(modena_error_occurred())
    {
        return modena_error();
    }

    // Allocate memory and fetch arg positions
    modena_inputs_t *inputs = modena_inputs_new(model);
    modena_outputs_t *outputs = modena_outputs_new(model);

    size_t ppos = modena_model_inputs_argPos(model, "p");
    size_t TPos = modena_model_inputs_argPos(model, "T");

    size_t DPos = modena_model_outputs_argPos(model, "D[A]");

    modena_model_argPos_check(model);

    const double p = 1e5;
    const double T = 290;

    // Set input vector
    modena_inputs_set(inputs, ppos, p);
    modena_inputs_set(inputs, TPos, T);

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
    double D = modena_outputs_get(outputs, 0);

    cout << "p = " << p << " T = " << T << " D = " << D << endl;

    modena_inputs_destroy(inputs);
    modena_outputs_destroy(outputs);
    modena_model_destroy(model);

    return 0;
}
