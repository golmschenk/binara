#ifdef IS_PYTHON_LIBRARY_BUILD
#include <Python.h>
#include "python_interrupt_handling.h"
#include <nanobind/nanobind.h>

extern "C" void check_for_and_handle_python_interrupt()
{
    if (PyErr_CheckSignals() != 0)
    {
        throw nanobind::python_error();
    }
}

#else
extern "C" void check_for_and_handle_python_interrupt()
{

}
#endif
