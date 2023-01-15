//

#ifndef PGEN_LIB_PGENIO_H
#define PGEN_LIB_PGENIO_H
#include "pgenContext.h"

PgenContext *openPgen (const char* cFilename, const long numberOfVariants, const long sampleCount);
void closePgen (const PgenContext * const pGenMeta);

#endif //PGEN_LIB_PGENIO_H
