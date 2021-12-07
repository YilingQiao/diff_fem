/*
Copyright 2021 by Inria, Mickaël Ly, Jean Jouve, Florence Bertails-Descoubes and
    Laurence Boissieux

This file is part of ProjectiveFriction.

ProjectiveFriction is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

ProjectiveFriction is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
ProjectiveFriction. If not, see <https://www.gnu.org/licenses/>.
*/
#include <chrono>

/// @brief Typedef for the chrono
using time_point_t = std::chrono::time_point<std::chrono::system_clock>;


#define TIMER_START(name)                       \
  const time_point_t timer_start_##name =       \
    std::chrono::high_resolution_clock::now();


#define TIMER_DELTA(name)                                               \
  ( std::chrono::high_resolution_clock::now() - timer_start_##name )
  
#define TIMER_DURATION(name, precision)                                 \
  (std::chrono::duration_cast< std::chrono::precision >(                \
    TIMER_DELTA(name) ).count())

//#define TIMER_PRINT
