//
// Created by Chris Norman on 4/8/23.
//

#include <map>
#include <stdio.h>
#include <string>
#include "pgenException.h"
#include "plink2_base.h"
#include "pgenUtils.h"

namespace pgenlib {

    std::map<plink2::PglErr::ec, std::string> pglErrMap {
            {plink2::PglErr::ec::kPglRetSuccess, "kPglRetSuccess"},
            {plink2::PglErr::ec::kPglRetSkipped, "kPglRetSkipped"},
            {plink2::PglErr::ec::kPglRetNomem, "kPglRetNomem"},
            {plink2::PglErr::ec::kPglRetOpenFail,"kPglRetOpenFail"},
            {plink2::PglErr::ec::kPglRetReadFail, "kPglRetReadFail"},
            {plink2::PglErr::ec::kPglRetWriteFail, "kPglRetWriteFail"},
            {plink2::PglErr::ec::kPglRetMalformedInput, "kPglRetMalformedInput"},
            {plink2::PglErr::ec::kPglRetInconsistentInput, "kPglRetInconsistentInput"},
            {plink2::PglErr::ec::kPglRetInvalidCmdline, "kPglRetInvalidCmdline"},
            {plink2::PglErr::ec::kPglRetThreadCreateFail, "kPglRetThreadCreateFail"},
            {plink2::PglErr::ec::kPglRetNetworkFail, "kPglRetNetworkFail"},
            {plink2::PglErr::ec::kPglRetVarRecordTooLarge, "kPglRetVarRecordTooLarge"},
            {plink2::PglErr::ec::kPglRetUnsupportedInstructions, "kPglRetUnsupportedInstructions"},
            {plink2::PglErr::ec::kPglRetDegenerateData, "kPglRetDegenerateData"},
            {plink2::PglErr::ec::kPglRetDecompressFail, "kPglRetDecompressFail"}, // also distinguish this from MalformedInput
            {plink2::PglErr::ec::kPglRetRewindFail, "kPglRetRewindFail"},
            {plink2::PglErr::ec::kPglRetGpuFail, "kPglRetGpuFail"},
            {plink2::PglErr::ec::kPglRetSampleMajorBed, "kPglRetSampleMajorBed"},
            {plink2::PglErr::ec::kPglRetNomemCustomMsg, "kPglRetNomemCustomMsg"},
            {plink2::PglErr::ec::kPglRetInternalError, "kPglRetInternalError"},
            {plink2::PglErr::ec::kPglRetWarningErrcode, "kPglRetWarningErrcode"},
            {plink2::PglErr::ec::kPglRetImproperFunctionCall, "kPglRetImproperFunctionCall"},
            {plink2::PglErr::ec::kPglRetNotYetSupported, "kPglRetNotYetSupported"},
            // These are only for internal use.  If any of these reach the top level
            // instead of being handled or converted to another error code, that's a bug,
            // and plink2 prints a message to that effect.
            {plink2::PglErr::ec::kPglRetHelp, "kPglRetHelp"},
            {plink2::PglErr::ec::kPglRetLongLine, "kPglRetLongLine"},
            {plink2::PglErr::ec::kPglRetEof, "kPglRetEof"}
    };

    template<typename MAP>
    const typename MAP::mapped_type& getOrDefault(
            const MAP& m,
            const typename MAP::key_type& key,
            const typename MAP::mapped_type& defval)
    {
        typename MAP::const_iterator it = m.find(key);
        if (it == m.end())
            return defval;
        return it->second;
    }

    // Throw an exception for an error that originated in underlying *PLINK* code.
    void throwOnPglErr(const plink2::PglErr pglErr, const char* message) {
        if (pglErr) { // PglErr evaluates to true if there was an error
            const int MSG_BUF_SIZE = 4096;
            char msgbuf[MSG_BUF_SIZE];
            std::string s = getOrDefault(pglErrMap, pglErr, "Unrecognized PglErr");
            snprintf(msgbuf,
                    MSG_BUF_SIZE,
                    "%s (PglErr: %d %s)",
                    message,
                    int32_t(pglErr), // include the integer value just in case...
                    s.c_str());
            throw PgenException(msgbuf);
        }
    }

} // pgenlib