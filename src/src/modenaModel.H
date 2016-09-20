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

@endcond
@file
Low-level interface library
@author    Henrik Rusche
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@defgroup  C++_interface_library
MoDeNa C++ interface library
*/

#ifndef __MODENAMODEL_H__
#define __MODENAMODEL_H__

#include <exception>
#include <string>
#include "modena.h"

using namespace std;

class modenaException
:
    public exception
{
    // Private data

    int errorCode_;

public:

    modenaException(int errorCode)
    :
        errorCode_(errorCode)
    {}

    int errorCode () const
    {
        return errorCode_;
    }

    virtual const char * what () const throw ()
    {
        return "Modena Exception";
    }
};


class modenaModel
{
    // Private data

        // Pointer to MoDeNa model
        modena_model_t *model_;

        modena_inputs_t *inputs_;
        modena_outputs_t *outputs_;

public:
 
    // Constructors

        //- Construct from model name
        modenaModel
        (
            const string& name
        )
        {
            // Instantiate a model
            model_ = modena_model_new(name.c_str());

            if(modena_error_occurred())
            {
                throw modenaException(modena_error());
            }

            // Allocate memory and fetch arg positions
            inputs_ = modena_inputs_new(model_);
            outputs_ = modena_outputs_new(model_);
        }

    // Destructor

        ~modenaModel()
        {
            modena_inputs_destroy(inputs_);
            modena_outputs_destroy(outputs_);
            modena_model_destroy(model_);
        }

    // Member Functions

        inline size_t inputs_argPos(const string& name) const
        {
            return modena_model_inputs_argPos(model_, name.c_str());
        }

        inline size_t outputs_argPos(const string& name) const
        {
            return modena_model_outputs_argPos(model_, name.c_str());
        }

        inline size_t inputs_set(const size_t i, double x) const
        {
            modena_inputs_set(inputs_, i, x);
        }

        double outputs_get(const size_t i) const
        {
           return modena_outputs_get(outputs_, i);
        }

        inline void argPos_check() const
        {
            modena_model_argPos_check(model_);
        }

        inline void call() const
        {
            modena_model_call(model_, inputs_, outputs_);

            if(modena_error_occurred())
            {
                throw modenaException(modena_error());
            }
        }

};

#endif /* __MODENAMODEL_H__ */
