target_sources(chargefw_core_objects
        PRIVATE
        atomic_charges.cpp
        charge_collection.cpp

        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        FILES
        ${PROJECT_SOURCE_DIR}/include/chargefw/charges/atomic_charges.h
        ${PROJECT_SOURCE_DIR}/include/chargefw/charges/charge_collection.h
)
