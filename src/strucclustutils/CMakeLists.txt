set(strucclustutils_source_files
        strucclustutils/tmalign.cpp
        strucclustutils/aln2tmscore.cpp
        strucclustutils/structcreatedb.cpp
        strucclustutils/pareunaligner.cpp
        strucclustutils/PareunAlign.cpp
        strucclustutils/PareunAlign.h
        PARENT_SCOPE
        )

add_library(gemmiwrapper OBJECT
    GemmiWrapper.cpp
    GemmiWrapper.h
    )

get_target_property(COMPILE_TMP mmseqs-framework COMPILE_FLAGS)
get_target_property(LINK_TMP mmseqs-framework LINK_FLAGS)
set_target_properties(gemmiwrapper PROPERTIES COMPILE_FLAGS "${COMPILE_TMP} -w" LINK_FLAGS "${LINK_TMP} -w")
restore_exceptions(gemmiwrapper)
