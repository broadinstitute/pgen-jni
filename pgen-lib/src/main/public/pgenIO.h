//

#ifndef PGEN_LIB_PGENIO_H
#define PGEN_LIB_PGENIO_H
#include "pgenMeta.h"

PGEN_META *openPgen (const char* cFilename, long numberOfVariants, long sampleCount);
void closePgen (PGEN_META * const pGenMeta);

#endif //PGEN_LIB_PGENIO_H
