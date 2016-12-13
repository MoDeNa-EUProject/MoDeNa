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

    The Modena interface library is free software; you can redistribute it
    and/or modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    Modena is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.

    This file should help to compile models on systems, where MoDeNa is not
    installed.
*/
#include <stddef.h>
#include <stdbool.h>
int *modena_inputs_new(int *self);
int *modena_inputs_new(int *self)
{
    int *i;
    return i;
}
int *modena_outputs_new(int *self);
int *modena_outputs_new(int *self)
{
    int *i;
    return i;
}
size_t modena_model_inputs_argPos(int *self, const char *name);
size_t modena_model_inputs_argPos(int *self, const char *name)
{
    size_t argPos;
    return argPos;
}
size_t modena_model_outputs_argPos(int *self, const char *name);
size_t modena_model_outputs_argPos(int *self, const char *name)
{
    size_t argPos;
    return argPos;
}
void modena_model_argPos_check(int *self);
void modena_model_argPos_check(int *self)
{

}
void modena_inputs_set(int *self, const size_t i, double x);
void modena_inputs_set(int *self, const size_t i, double x)
{

}
void modena_inputs_destroy(int *inputs);
void modena_inputs_destroy(int   *self)
{

}
void modena_outputs_destroy(int *inputs);
void modena_outputs_destroy(int   *self)
{

}
void modena_model_destroy(int *inputs);
void modena_model_destroy(int   *self)
{

}
double modena_outputs_get(int *self, const size_t i);
double modena_outputs_get(int *self, const size_t i)
{
    double d;
    return d;
}
int *modena_model_new(const char *modelId);
int *modena_model_new(const char *modelId)
{
    int *i;
    return i;
}
int modena_model_call(int *model, int *inputs, int *outputs);
int modena_model_call(int *model, int *inputs, int *outputs)
{
    int i;
    return i;
}
bool modena_error_occurred();
bool modena_error_occurred()
{
    bool b;
    return b;
}
int modena_error();
int modena_error()
{
    int i;
    return i;
}
