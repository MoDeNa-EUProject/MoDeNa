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
#include <vector>

#include "modena.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Modena
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class modenaException
:
    public std::exception
{
    // Private data

    int errorCode_;

public:

    modenaException(int errorCode)
    :
        errorCode_(errorCode)
    {}

    int errorCode() const
    {
        return errorCode_;
    }

    virtual const char*  what() const throw()
    {
        return "Modena Exception";
    }
};

/**
 * @brief    C++ wrapper class for MoDeNa Models
 * @author   Henrik Rusche
 * @date     2016
 * @details
 *           The class
 *
*/
class modenaModel
{

private :
        // Pointer to MoDeNa model
        modena_model_t *model_;

        modena_inputs_t *inputs_;
        modena_outputs_t *outputs_;

public :

    // Constructors

        modenaModel(const std::string& name); //- Construct from model name

    // Destructor

        ~modenaModel();

    // Overloaded Operators

        // Subscript
        double& operator[](const size_t& idx);
        const double& operator [](const size_t& idx) const;

        void operator()(void);

    // Member Functions (inline is technically redundant in declaration)

        inline size_t inputs_size() const;
        inline size_t outputs_size() const;
        inline size_t parameters_size() const;

        inline size_t inputs_argPos(const std::string& name) const;
        inline size_t outputs_argPos(const std::string& name) const;

        inline void inputs_set(const size_t i, double x) const;

        inline double inputs_get(const size_t i) const;
        inline double outputs_get(const size_t i) const;

        inline void argPos_check() const;
        inline void call() const;

        inline std::vector<std::string> inputs_names() const;
        inline std::vector<std::string> outputs_names() const;
        inline std::vector<std::string> parameters_names() const;
};

//- Construct from model name
modenaModel::modenaModel(const std::string& name)
{
    // Instantiate a model
    model_ = modena_model_new(name.c_str());

    if(modena_error_occurred())
        throw modenaException(modena_error());

    // Allocate memory and fetch arg positions
    inputs_ = modena_inputs_new(model_);
    outputs_ = modena_outputs_new(model_);
}

/**
 * @brief  Destructor
 *
*/
modenaModel::~modenaModel()
{
    modena_inputs_destroy(inputs_);
    modena_outputs_destroy(outputs_);
    modena_model_destroy(model_);
}

// Operators
double& modenaModel::operator[](const size_t& idx)
{
  std::cout << "Setter" << std::endl;
  return this->inputs_->inputs[idx];
}

const double& modenaModel::operator[](const size_t& idx) const
{
    const double ret = outputs_get(idx);
    return ret;
}

void modenaModel::operator()(void)
{
    this->call();
}

// Member Functions

inline size_t modenaModel::inputs_size() const
{
    return modena_model_inputs_size(model_);
}

inline size_t modenaModel::outputs_size() const
{
    return modena_model_outputs_size(model_);
}

inline size_t modenaModel::parameters_size() const
{
    return modena_model_parameters_size(model_);
}

inline size_t modenaModel::inputs_argPos(const std::string& name) const
{
    return modena_model_inputs_argPos(model_, name.c_str());
}

inline size_t modenaModel::outputs_argPos(const std::string& name) const
{
    return modena_model_outputs_argPos(model_, name.c_str());
}

inline void modenaModel::inputs_set(const size_t i, double x) const
{
    modena_inputs_set(inputs_, i, x);
}

inline double modenaModel::inputs_get(const size_t i) const
{
    return modena_inputs_get(inputs_, i);
}

inline double modenaModel::outputs_get(const size_t i) const
{
    return modena_outputs_get(outputs_, i);
}

inline void modenaModel::argPos_check() const
{
    modena_model_argPos_check(model_);
}

inline std::vector<std::string> modenaModel::inputs_names() const
{
    std::vector<std::string> v;

    const char** iNames = modena_model_inputs_names(model_);
    //char** iNames = modena_model_inputs_names(model_);
    for(size_t i=0; i<modena_model_inputs_size(model_); i++)
        v.push_back(std::string(iNames[i]));

    return v;
}

inline std::vector<std::string> modenaModel::outputs_names() const
{
    std::vector<std::string> v;

    const char** oNames = modena_model_outputs_names(model_);
    // char** oNames = modena_model_outputs_names(model_);
    for(size_t i=0; i<modena_model_outputs_size(model_); i++)
        v.push_back(std::string(oNames[i]));

    return v;
}

inline std::vector<std::string> modenaModel::parameters_names() const
{
    std::vector<std::string> v;

    const char** pNames = modena_model_parameters_names(model_);
    // char** pNames = modena_model_parameters_names(model_);
    for(size_t i=0; i<modena_model_parameters_size(model_); i++)
        v.push_back(std::string(pNames[i]));

    return v;
}

inline void modenaModel::call() const
{
    modena_model_call(model_, inputs_, outputs_);

    if(modena_error_occurred())
        throw modenaException(modena_error());
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Modena

#endif /* __MODENAMODEL_H__ */

// ************************************************************************* //
