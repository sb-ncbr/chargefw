target_sources(chargefw_core_objects
        PRIVATE
        classification_result.cpp
        parameter_classification.cpp
        parameter_classifier.cpp
        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        FILES
        ${PROJECT_SOURCE_DIR}/include/chargefw/parameters/classification/classification_result.h
        ${PROJECT_SOURCE_DIR}/include/chargefw/parameters/classification/parameter_classifier.h
        ${PROJECT_SOURCE_DIR}/include/chargefw/parameters/classification/parameter_classification.h
)