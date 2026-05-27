## Code style
- For *new* files, format the code using clang-format with the configuratio in opm-simulators/.clang-format.
- For *existing* files, follow the existing code style in that file. 
- Prefer using `auto`.
- Naming: Write out abbreviations in full, except for well-known ones (e.g., `std`, `ptr`, `iter`). If the type does not have an abbreviation, the variable should not either.
- Use std::format for formatting strings, and avoid using stringstream or printf.
- Use OPM_THROW(excp, msg) for throwing exceptions, where excp is a standard exception type (e.g., std::runtime_error) and msg is a string message.
- Run clang-tidy with the configuration in .clang-tidy. Try to fix all warnings. 
- Use `const` and `constexpr` where appropriate.
- Use C++20 except modules.
## Testing
- We use Boost test
- You can tests under opm-simulators/tests

## Copyright and licensing
- Keep TODO ADD YEAR AND NAME OF AUTHOR to the copyright notice in new files, this should be handled manually. 
- Add the following copyright header:
```cpp
/*
  Copyright TODO ADD YEAR AND NAME OF AUTHOR

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
```
