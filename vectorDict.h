#ifndef VECTORDICT_H
#define VECTORDICT_H

#include <vector>
#include "TObject.h"

// Declaration of the custom type for dictionary generation
#ifdef __CLING__
#pragma link C++ class std::vector<std::vector<UChar_t>>+;
#endif

#endif  // VECTORDICT_H